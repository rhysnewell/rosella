use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use log::{debug, info, warn};

use crate::external::hmmer_engine::HmmerEngine;
use crate::quality::{Quality, orfs};

mod annotator;
pub mod cache;
pub mod checkm;
pub mod fragments;
mod hit;
pub mod hmm_table;
mod proteins;
mod readings;
pub mod replicon;
pub mod sets;
mod strain;
mod table;
mod walk;

pub use annotator::Annotator;
pub use hit::{Hit, Place};
pub use readings::{CHECKM, COPY_BINS, Coding, Duplicate, GTDB, Reading, SetReading};
pub use strain::{heterogeneity, identity, rebuild};
pub use table::MarkerSet;
pub use walk::Shed;

pub(crate) const HMM_GZ: &[u8] = include_bytes!("../../data/gtdb_markers.hmm.gz");

// hmmsearch refuses a target past this, and an ORF this long is an uncovered N span in a gold
// standard assembly or a scaffold gap, never a marker gene.
const MAX_SEARCH_RESIDUES: usize = 100_000;

/// A marker gene cut by a contig end is still that marker gene, but two halves of one gene on
/// two contigs are not two copies, so presence and duplication read different columns.
#[derive(Clone, Copy, Debug)]
pub struct MarkerRules {
    pub fragment_span: f64,
}

impl Default for MarkerRules {
    fn default() -> Self {
        Self {
            fragment_span: fragments::DEFAULT_SPAN,
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Tally {
    complete: u32,
    any: u32,
}

pub struct MarkerAnnotation {
    rows: cache::Rows,
    set: MarkerSet,
    proteins: Option<proteins::Proteins>,
    checkm: bool,
}

impl MarkerAnnotation {
    pub fn build(annotator: &Annotator, band: std::ops::Range<usize>) -> Result<Self> {
        let (assembly, rules) = (annotator.assembly.as_str(), annotator.rules);
        let set = MarkerSet::embedded();
        let cached = annotator
            .cache
            .as_deref()
            .map(|directory| {
                cache::key(assembly, band.start, rules.fragment_span).map(|key| (directory, key))
            })
            .transpose()?;
        let held = cached
            .as_ref()
            .and_then(|(directory, key)| cache::find(directory, key))
            .and_then(|(path, floor)| match cache::read(&path, &set) {
                Ok(rows) => {
                    info!("Read the marker annotation from {}", path.display());
                    Some((floor, (path, rows)))
                }
                Err(error) => {
                    warn!("Ignoring {}: {error}", path.display());
                    None
                }
            });
        let (ceiling, held) = match held {
            Some((floor, (_, rows))) if floor <= band.start => {
                return Ok(Self {
                    rows: rows.at_least(band.start),
                    set,
                    proteins: None,
                    checkm: annotator.checkm,
                });
            }
            Some((floor, held)) if floor <= band.end => (floor, Some(held)),
            _ => (band.end, None),
        };
        let (mut rows, proteins) = annotate(annotator, band.start..ceiling, &set)?;
        let whole = held.is_some() || band.end == usize::MAX;
        let superseded = match held {
            Some((path, held)) => {
                rows = rows.fill_from(held, ceiling)?;
                Some(path)
            }
            None => None,
        };
        if let Some((directory, key)) = cached.as_ref().filter(|_| whole) {
            let path = cache::write_path(directory, key);
            match cache::write(&path, key, &set, &rows) {
                Ok(()) => {
                    info!("Wrote the marker annotation to {}", path.display());
                    if let Some(old) = superseded.filter(|old| *old != path) {
                        let _ = std::fs::remove_file(old);
                    }
                }
                Err(error) => warn!("Could not write {}: {error}", path.display()),
            }
        }
        Ok(Self {
            rows,
            set,
            proteins,
            checkm: annotator.checkm,
        })
    }

    /// Contigs under the annotated band are filled later or not at all, so a missing contig
    /// reads as featureless rather than refusing the whole table.
    pub fn select_present(self, names: &[String]) -> Result<ContigMarkers> {
        self.selected(names, true)
    }

    pub fn select(self, names: &[String]) -> Result<ContigMarkers> {
        self.selected(names, false)
    }

    fn selected(mut self, names: &[String], tolerate_missing: bool) -> Result<ContigMarkers> {
        let index = self
            .rows
            .names
            .iter()
            .enumerate()
            .map(|(position, name)| (name.as_str(), position))
            .collect::<HashMap<_, _>>();
        let mut per_contig = Vec::with_capacity(names.len());
        let mut shapes = Vec::with_capacity(names.len());
        let mut checkm = Vec::with_capacity(names.len());
        let mut placed = vec![None; self.rows.names.len()];
        for (contig, name) in names.iter().enumerate() {
            match index.get(name.as_str()) {
                Some(position) => {
                    placed[*position] = Some(contig);
                    per_contig.push(std::mem::take(&mut self.rows.hits[*position]));
                    shapes.push(Some(self.rows.shapes[*position]));
                    checkm.push(self.rows.checkm[*position].take());
                }
                None if tolerate_missing => {
                    per_contig.push(Default::default());
                    shapes.push(None);
                    checkm.push(None);
                }
                None => anyhow::bail!("the marker table does not hold {name}"),
            }
        }
        let markers = ContigMarkers::new(per_contig, self.set).with_shapes(shapes);
        Ok(match self.checkm {
            true => markers
                .with_checkm(checkm)
                .with_proteins(self.proteins.map(|kept| kept.renumbered(&placed))),
            false => markers,
        })
    }
}

// Only contigs in `band` are called and searched. Those at or above its end come back unscored,
// for the caller to fill from a cached annotation.
fn annotate(
    annotator: &Annotator,
    band: std::ops::Range<usize>,
    set: &MarkerSet,
) -> Result<(cache::Rows, Option<proteins::Proteins>)> {
    // Annotating is most of a run and the temp directory is dropped on any failure, so the
    // search has to be known to work before the gene calling is paid for.
    HmmerEngine::check_installed()?;
    let engine = HmmerEngine::new(annotator.threads, annotator.shards);
    let rules = annotator.rules;

    let directory = tempfile::tempdir()?;
    let hmm = directory.path().join("markers.hmm");
    inflate(HMM_GZ, &hmm)?;

    let (walked, called, pieces, shapes) = {
        let _timer = crate::timing::scope("genes");
        let mut sink = engine.protein_shards(directory.path())?;
        let mut called: Vec<orfs::Orf> = Vec::new();
        let mut shapes: Vec<replicon::Shape> = Vec::new();
        let walked = orfs::call_over(&annotator.assembly, band.clone(), |batch| {
            for mut orf in batch {
                if shapes.len() <= orf.contig {
                    shapes.resize(orf.contig + 1, replicon::Shape::default());
                }
                shapes[orf.contig].add(orf.bases);
                if searchable(&orf.protein) {
                    sink.write(called.len(), &orf.protein)?;
                }
                if !orf.partial() {
                    orf.protein = String::new();
                }
                called.push(orf);
            }
            Ok(())
        })?;
        let pieces = sink.finish()?;
        shapes.resize(walked.names.len(), replicon::Shape::default());
        let banded = walked.lengths.iter().filter(|length| **length < band.end);
        info!(
            "Called {} genes over {} contigs",
            called.len(),
            banded.count()
        );
        (walked, called, pieces, shapes)
    };

    let bars = fragments::gathering(&hmm)?;
    let table = {
        let _timer = crate::timing::scope("search");
        // CheckM reads this table for the models both panels share, so it reports down to
        // CheckM's lowest bar as well.
        let floor = fragments::floor(&bars, rules.fragment_span).min(set.checkm.floor());
        engine.search(
            &hmm,
            &pieces,
            directory.path(),
            &format!("{floor:.2}"),
            "markers",
        )?
    };
    let mut hits = fragments::complete(&table, &bars);
    {
        let _timer = crate::timing::scope("fragments");
        let cut = |protein: usize| {
            called.get(protein).is_some_and(|orf| orf.partial()) && !hits.contains_key(&protein)
        };
        let rescued = fragments::accepted(&table, &bars, rules.fragment_span, cut);
        debug!("{} markers rescued from cut genes", rescued.len());
        hits.extend(rescued);
    }

    let mut per_contig = vec![Vec::new(); walked.names.len()];
    for (protein, best) in hits {
        let Some(marker) = set.id(&best.model) else {
            continue;
        };
        let Some(orf) = called.get(protein) else {
            continue;
        };
        per_contig[orf.contig].push(Hit::called(marker, orf, &best));
    }
    in_marker_order(&mut per_contig);
    let carriers = per_contig.iter().filter(|hits| !hits.is_empty()).count();
    debug!(
        "{carriers} of {} contigs carry a single copy marker",
        walked.names.len()
    );
    let proteins = annotator
        .checkm
        .then(|| {
            let covered = walked.lengths.iter().map(|length| *length < band.end);
            proteins::Proteins::keep(directory, pieces, &table, &called, covered, engine)
        })
        .transpose()?;
    Ok((
        cache::Rows {
            checkm: vec![None; walked.names.len()],
            names: walked.names,
            lengths: walked.lengths,
            hits: per_contig,
            shapes,
        },
        proteins,
    ))
}

/// The search returns its hits in hash order, and the cache and the report are written from
/// these, so two annotations of one assembly would not diff against each other.
fn in_marker_order(per_contig: &mut [Vec<Hit>]) {
    for hits in per_contig {
        hits.sort_unstable_by_key(|hit| {
            (
                hit.marker,
                hit.partial,
                hit.place.gene_begin,
                hit.place.reverse,
            )
        });
    }
}

pub struct ContigMarkers {
    per_contig: Vec<Vec<Hit>>,
    checkm: Vec<Option<Vec<checkm::Copies>>>,
    proteins: Vec<proteins::Proteins>,
    shapes: Vec<replicon::Shape>,
    annotated: Vec<bool>,
    lengths: Vec<usize>,
    set: MarkerSet,
}

impl crate::quality::Scorer for ContigMarkers {
    fn features(&self, contigs: &[usize]) -> std::collections::HashSet<u32> {
        self.counts(contigs)
            .iter()
            .enumerate()
            .filter(|(_, tally)| tally.any > 0)
            .map(|(marker, _)| marker as u32)
            .collect()
    }

    /// Read against whichever lineage the bin's pattern of absences fits, since a reduced
    /// genome is missing markers a whole one of another lineage would carry.
    fn score(&self, contigs: &[usize]) -> Quality {
        let Some((chosen, counts)) = self.chosen(contigs) else {
            return Quality::default();
        };
        match self.gtdb_on(&counts, chosen) {
            reading if reading.markers == 0 => Quality::default(),
            reading => reading.quality(chosen),
        }
    }
}

fn observed(counts: &[Tally]) -> Vec<u16> {
    counts
        .iter()
        .enumerate()
        .filter(|(_, tally)| tally.any > 0)
        .map(|(marker, _)| marker as u16)
        .collect()
}

impl ContigMarkers {
    pub fn new(mut per_contig: Vec<Vec<Hit>>, set: MarkerSet) -> Self {
        in_marker_order(&mut per_contig);
        Self {
            per_contig,
            checkm: Vec::new(),
            proteins: Vec::new(),
            shapes: Vec::new(),
            annotated: Vec::new(),
            lengths: Vec::new(),
            set,
        }
    }

    pub fn set_name(&self, set: u16) -> &str {
        self.set.sets.name(set as usize)
    }

    pub fn report<'a>(
        &self,
        placed: impl IntoIterator<Item = (&'a str, usize)>,
        names: &[String],
        path: &Path,
    ) -> Result<()> {
        let mut sink = BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            sink,
            "bin\tcontig\tmodel\tpartial\tmodel_from\tmodel_to\tmodel_length\tprotein_from\t\
             protein_to\tscore\tgene_begin\tgene_end\tstrand\tcut_left\tcut_right"
        )?;
        for (bin, contig) in placed {
            for hit in &self.per_contig[contig] {
                writeln!(
                    sink,
                    "{bin}\t{}\t{}\t{}\t{}",
                    names[contig],
                    self.set.name(hit.marker),
                    u8::from(hit.partial),
                    hit.fields('\t')
                )?;
            }
        }
        sink.flush()?;
        Ok(())
    }

    pub fn with_lengths(mut self, lengths: Vec<usize>) -> Self {
        self.lengths = lengths;
        self
    }

    pub fn with_shapes(mut self, shapes: Vec<Option<replicon::Shape>>) -> Self {
        self.annotated = shapes.iter().map(Option::is_some).collect();
        self.shapes = shapes.into_iter().map(Option::unwrap_or_default).collect();
        self
    }

    pub fn with_checkm(mut self, checkm: Vec<Option<Vec<checkm::Copies>>>) -> Self {
        self.checkm = checkm;
        self
    }

    fn with_proteins(mut self, kept: Option<proteins::Proteins>) -> Self {
        self.proteins.extend(kept);
        self
    }

    pub fn departures<'a>(
        &self,
        bins: impl IntoIterator<Item = &'a [usize]>,
        composition: ndarray::ArrayView2<f64>,
        depths: ndarray::ArrayView2<f64>,
    ) -> replicon::Departures {
        if self.shapes.len() != self.lengths.len() {
            return replicon::Departures::default();
        }
        let carries = self
            .per_contig
            .iter()
            .map(|hits| !hits.is_empty())
            .collect::<Vec<_>>();
        replicon::departures(
            &self.shapes,
            &carries,
            &self.lengths,
            bins,
            composition,
            depths,
        )
    }

    fn chosen(&self, contigs: &[usize]) -> Option<(usize, Vec<Tally>)> {
        let counts = self.counts(contigs);
        let chosen = self
            .set
            .sets
            .choose(&observed(&counts), self.bin_bp(contigs))?;
        Some((chosen, counts))
    }

    pub fn points(&self, contigs: &[usize], weight: f64) -> Vec<(usize, f64)> {
        let Some((chosen, counts)) = self.chosen(contigs) else {
            return Vec::new();
        };
        counts
            .iter()
            .enumerate()
            .filter(|(marker, _)| self.set.sets.holds(chosen, *marker))
            .map(|(marker, tally)| {
                let extra = f64::from(tally.complete.saturating_sub(1));
                let present = f64::from(u8::from(tally.any >= 1));
                (marker, 100.0 * (present - weight * extra))
            })
            .collect()
    }

    pub fn hit_count(&self, contigs: &[usize]) -> usize {
        contigs
            .iter()
            .map(|contig| self.per_contig[*contig].len())
            .sum()
    }

    fn bin_bp(&self, contigs: &[usize]) -> usize {
        contigs
            .iter()
            .filter_map(|contig| self.lengths.get(*contig))
            .sum()
    }

    /// Whether the contig carries a whole marker copy the rest of the bin does not, which is
    /// the difference between a home and a second copy of what is already there.
    pub fn completes(&self, contigs: &[usize], contig: usize) -> bool {
        self.repeats(contigs, contig) == Some(false)
    }

    /// None when the contig holds no whole marker of the bin's set, so it is no evidence either way.
    pub fn repeats(&self, contigs: &[usize], contig: usize) -> Option<bool> {
        self.repeated(contigs, contig, true, |_, _| true)
    }

    pub fn repeats_any(&self, contigs: &[usize], contig: usize) -> Option<bool> {
        self.repeated(contigs, contig, false, |_, _| true)
    }

    pub fn repeats_in_place(&self, contigs: &[usize], contig: usize) -> Option<bool> {
        self.repeated(contigs, contig, false, Hit::same_part)
    }

    fn repeated(
        &self,
        contigs: &[usize],
        contig: usize,
        whole_only: bool,
        same: impl Fn(&Hit, &Hit) -> bool,
    ) -> Option<bool> {
        let (chosen, _) = self.chosen(contigs)?;
        let counted = |hit: &&Hit| {
            !(whole_only && hit.partial) && self.set.sets.holds(chosen, hit.marker as usize)
        };
        let held = self.per_contig[contig]
            .iter()
            .filter(counted)
            .collect::<Vec<_>>();
        if held.is_empty() {
            return None;
        }
        let rest = contigs
            .iter()
            .filter(|member| **member != contig)
            .flat_map(|member| self.per_contig[*member].iter().filter(counted))
            .filter(|other| held.iter().any(|hit| hit.marker == other.marker))
            .collect::<Vec<_>>();
        Some(held.iter().all(|hit| {
            rest.iter()
                .any(|other| other.marker == hit.marker && same(hit, other))
        }))
    }

    fn whole(&self, contig: usize, chosen: usize) -> Vec<usize> {
        let mut held = self.per_contig[contig]
            .iter()
            .filter(|hit| !hit.partial && self.set.sets.holds(chosen, hit.marker as usize))
            .map(|hit| hit.marker as usize)
            .collect::<Vec<_>>();
        held.dedup();
        held
    }

    fn counts(&self, contigs: &[usize]) -> Vec<Tally> {
        let mut counts = vec![Tally::default(); self.set.len()];
        for contig in contigs {
            for hit in &self.per_contig[*contig] {
                let tally = &mut counts[hit.marker as usize];
                tally.any += 1;
                tally.complete += u32::from(!hit.partial);
            }
        }
        counts
    }
}

fn searchable(protein: &str) -> bool {
    !protein.is_empty() && protein.len() <= MAX_SEARCH_RESIDUES
}

fn inflate(compressed: &[u8], target: &Path) -> Result<()> {
    let mut decoder = flate2::read::GzDecoder::new(compressed);
    let mut sink = BufWriter::new(std::fs::File::create(target)?);
    std::io::copy(&mut decoder, &mut sink)?;
    sink.flush()?;
    Ok(())
}
