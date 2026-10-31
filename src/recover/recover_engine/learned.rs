use std::path::{Path, PathBuf};

use anyhow::Result;
use log::{debug, info, warn};
use ndarray::{Array2, Axis};

use crate::cli::RunArgs;
use crate::coverage::depth_file::{self, CONTIGS};
use crate::coverage::read_stats::{Placed, tie_share, tied_reads};
use crate::defaults::TRAINING_FLOOR;
use crate::learned::corpus::Inputs;
use crate::learned::labelled::{self, Labelled};
use crate::learned::{Learned, Space, gemm, kept_carriers, store, windows::Windows};
use crate::markers::ContigMarkers;
use crate::recover::inputs::{laid_out, row_lengths};
use crate::recover::recover_engine::RecoverEngine;
use crate::recover::settings::Settings;
use crate::tables::Tables;

struct Training {
    composition: Array2<f64>,
    coverage: Array2<f64>,
    names: Vec<String>,
    lengths: Vec<usize>,
    rows: Vec<usize>,
}

// A given attach floor loads its rows in name hash order, so those are picked by length. Rows
// still to load sit longest first, so the training set takes a run off their front.
fn training(tables: &Tables, cutoff: usize) -> Result<Training> {
    let floor = TRAINING_FLOOR.min(cutoff);
    let lengths = &tables.coverage.contig_lengths;
    let loaded = tables.tnf.kmer_table.nrows();
    let rows = (0..loaded)
        .filter(|row| lengths[*row] >= floor)
        .collect::<Vec<_>>();
    let unread = loaded..loaded + lengths[loaded..].partition_point(|length| *length >= floor);
    let mut composition = tables.tnf.kmer_table.select(Axis(0), &rows);
    let mut coverage = tables.coverage.table.select(Axis(0), &rows);
    if !unread.is_empty() {
        let (more, depth) = tables.read_rows(unread.clone())?;
        composition.append(Axis(0), more.view())?;
        coverage.append(Axis(0), depth.view())?;
    }
    let rows = rows.into_iter().chain(unread).collect::<Vec<_>>();
    Ok(Training {
        composition,
        coverage,
        names: crate::rows::reorder(&tables.coverage.contig_names, &rows),
        lengths: crate::rows::reorder(lengths, &rows),
        rows,
    })
}

// The space reads no labels, so it trains while the markers are still being annotated. Its own
// pool keeps its row work from queueing behind gene calling in the shared one.
pub(super) fn learn_space(
    tables: &mut Tables,
    settings: &Settings,
    dir: &Path,
    args: &RunArgs,
    given: bool,
    markers: Option<&ContigMarkers>,
) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(rayon::current_num_threads())
        .build()?
        .install(|| learn(tables, settings, dir, args, given, markers))
}

fn learn(
    tables: &mut Tables,
    settings: &Settings,
    dir: &Path,
    args: &RunArgs,
    given: bool,
    markers: Option<&ContigMarkers>,
) -> Result<()> {
    let _gpu = gemm::gpu(settings.gpu);
    let table = args.reads.coverage.coverage_file.as_deref();
    let entry = match &settings.space_cache {
        Some(cache) => {
            space_key(tables, settings, dir, table)?.map(|key| (store::path(cache, &key), key))
        }
        None => None,
    };
    let loaded = entry
        .as_ref()
        .and_then(|(path, key)| match Learned::load(path, key) {
            Ok(loaded) => loaded.inspect(|_| info!("Learned space: read from {}", path.display())),
            Err(error) => {
                warn!("Space cache: {error}");
                None
            }
        });
    let learned = match loaded {
        Some(learned) => learned,
        None => {
            // A space trained on given bins alone serves no later run, so a kept one is trained
            // on every contig recover would read.
            let learned = match (&entry, given) {
                (Some(_), true) => {
                    debug!("Learned space: training on the whole assembly to keep it");
                    let whole = laid_out(args, settings, row_lengths(settings), None)?;
                    train(&whole, settings, dir, None)?
                }
                _ => train(tables, settings, dir, markers)?,
            };
            let Some(learned) = learned else {
                return Ok(());
            };
            if let Some((path, key)) = &entry
                && let Err(error) = learned.save(path, key)
            {
                warn!("Space cache: {error} ({})", path.display());
            }
            learned
        }
    };
    let rows = {
        let _timer = crate::timing::scope("learned_embed");
        learned.embed(tables.tnf.kmer_table.view(), tables.coverage.table.view())
    };
    tables.tnf.kmer_table = rows;
    tables.learned = Some(Space::Windowed(Box::new(learned)));
    Ok(())
}

fn train(
    tables: &Tables,
    settings: &Settings,
    dir: &Path,
    markers: Option<&ContigMarkers>,
) -> Result<Option<Learned>> {
    let training = training(tables, settings.lengths.cutoff)?;
    let spans = {
        let _timer = crate::timing::scope("learned_windows");
        let windows = Windows::read(dir, &tables.coverage.sample_names)?;
        windows.spans(&training.names, &training.lengths, settings.seeds.seed)?
    };
    Learned::train(&Inputs {
        assembly: &tables.assembly,
        names: &training.names,
        lengths: &training.lengths,
        composition: training.composition.view(),
        coverage: training.coverage.view(),
        kmer_size: tables.tnf.kmer_size(),
        seed: settings.seeds.seed,
        epochs: settings.epochs,
        patience: settings.patience,
        stretches: settings.stretches,
        fast: settings.fast_train,
        spans: &spans,
        families: match (markers, settings.marker_steps) {
            (Some(markers), Some(percent)) => {
                let carriers = coded(markers, &training.rows, markers.carriers(&training.rows));
                tied(dir, tables, &training.names, &carriers, percent)?
            }
            _ => Vec::new(),
        },
    })
}

fn memberships(families: &[Vec<u32>]) -> usize {
    families.iter().map(Vec::len).sum()
}

fn coded(markers: &ContigMarkers, rows: &[usize], families: Vec<Vec<u32>>) -> Vec<Vec<u32>> {
    let mut carriers = families.iter().flatten().copied().collect::<Vec<_>>();
    carriers.sort_unstable();
    carriers.dedup();
    let contigs = carriers
        .iter()
        .map(|at| rows[*at as usize])
        .collect::<Vec<_>>();
    let prokaryotic = markers.prokaryotic(&contigs);
    let kept = kept_carriers(&families, |at| {
        carriers
            .binary_search(&at)
            .is_ok_and(|found| prokaryotic[found])
    });
    info!(
        "Marker carriers: {} of {} code densely enough for a prokaryote",
        memberships(&kept),
        memberships(&families)
    );
    kept
}

fn tied(
    dir: &Path,
    tables: &Tables,
    names: &[String],
    families: &[Vec<u32>],
    percent: f64,
) -> Result<Vec<Vec<u32>>> {
    let samples = &tables.coverage.sample_names;
    let placed = Placed::read(dir, samples, names)?;
    let wanted = placed.wanted(families.iter().flatten().map(|row| *row as usize));
    let tied = tied_reads(dir, samples, &wanted)?;
    let kept = kept_carriers(families, |row| {
        placed.at[row as usize]
            .is_some_and(|contig| 100.0 * tie_share(tied[contig], placed.reads[contig]) >= percent)
    });
    info!(
        "Marker carriers: {} of {} have {percent}% of their reads tied",
        memberships(&kept),
        memberships(families)
    );
    Ok(kept)
}

fn space_key(
    tables: &Tables,
    settings: &Settings,
    dir: &Path,
    table: Option<&str>,
) -> Result<Option<String>> {
    let Some(content) = tables.content else {
        return Ok(None);
    };
    let mut fields = vec![
        store::FORMAT.to_string(),
        env!("CARGO_PKG_VERSION").to_string(),
        env!("ROSELLA_TRAINING_DIGEST").to_string(),
        format!("{content:016x}"),
        TRAINING_FLOOR.min(settings.lengths.cutoff).to_string(),
        format!("{:?}", settings.lengths.attach_floor),
        tables.tnf.kmer_size().to_string(),
        settings.seeds.seed.to_string(),
        settings.epochs.to_string(),
        settings.patience.to_string(),
        settings.stretches.to_string(),
        settings.fast_train.to_string(),
        settings.gpu.to_string(),
    ];
    fields.extend(
        settings
            .marker_steps
            .map(|percent| format!("marker_steps:{percent}")),
    );
    let files = tables
        .coverage
        .sample_names
        .iter()
        .map(|sample| depth_file::path(dir, sample))
        .chain([dir.join(CONTIGS)])
        .chain(table.map(PathBuf::from));
    for file in files {
        let found = std::fs::metadata(&file)?;
        fields.push(format!(
            "{}:{}:{}",
            file.file_name().unwrap_or_default().to_string_lossy(),
            found.len(),
            crate::markers::cache::stamp(&found)
        ));
    }
    Ok(Some(fields.join("\t")))
}

impl RecoverEngine {
    // A pass in k-mer space lends its near complete bins as labels, and the run then reads
    // composition from the space they teach.
    pub(super) fn learn_labelled(&mut self, contigs: &[usize]) -> Result<()> {
        let lengths = self.lengths();
        let cutoff = self.settings.lengths.cutoff;
        let long = contigs.partition_point(|contig| lengths[*contig] >= cutoff);
        let _gpu = gemm::gpu(self.settings.gpu);
        let labels = self.first_labels(&contigs[..long])?;
        let tables = &self.tables;
        let training = training(tables, cutoff)?;
        let learned = Labelled::train(&labelled::corpus::Inputs {
            assembly: &tables.assembly,
            names: &training.names,
            lengths: &training.lengths,
            composition: training.composition.view(),
            coverage: training.coverage.view(),
            labels: &labels,
            halves: &tables.halves,
            halves_floor: 2 * cutoff,
            kmer_size: tables.tnf.kmer_size(),
            seed: self.settings.seeds.seed,
            epochs: self.settings.epochs,
            patience: self.settings.patience,
            fast: self.settings.fast_train,
        })?;
        let Some(learned) = learned else {
            return Ok(());
        };
        let rows = {
            let _timer = crate::timing::scope("learned_embed");
            learned.embed(tables.tnf.kmer_table.view(), tables.coverage.table.view())
        };
        self.tables.tnf.kmer_table = rows;
        self.tables.learned = Some(Space::Labelled(Box::new(learned)));
        Ok(())
    }
}
