use std::{
    collections::{HashMap, HashSet},
    path,
    thread::JoinHandle,
};

use anyhow::Result;
use log::{debug, warn};

use crate::{
    cli::{EnsembleArgs, RecoverArgs, ReportPaths, RunArgs},
    clustering::clusterer::{Partitioning, conserved},
    embedding::features::ContigFeatures,
    markers::{Annotator, ContigMarkers},
    recover::inputs::{Inputs, read_inputs},
    recover::origin::{BinSet, Given, Origin, Placed},
    recover::settings::Settings,
    recover::stage_report::StageReport,
    tables::Tables,
};

mod attach;
mod bin_writer;
mod incoherent;
mod learned;
mod offered;
mod partition;
mod sources;
mod stages;
mod tables;
mod weight;

use tables::Scoring;
use weight::Settled;

pub const UNBINNED: &str = "unbinned";

pub fn run_recover(args: &RecoverArgs) -> Result<()> {
    RecoverEngine::new(&args.run, Origin::Partition(args))?.run()
}

pub fn run_ensemble(args: &EnsembleArgs) -> Result<()> {
    let names = set_names(&args.bin_directories);
    let sets = args
        .bin_directories
        .iter()
        .zip(names)
        .map(|(source, name)| {
            let given = match path::Path::new(source).is_dir() {
                true => Given::read(
                    &crate::bins::fasta_in(source)?,
                    String::new(),
                    HashMap::new(),
                    None,
                )?,
                false => Given::table(path::Path::new(source))?,
            };
            Ok((name, given))
        })
        .collect::<Result<Vec<_>>>()?;
    debug!(
        "Bin sets: {} ({} bins)",
        sets.len(),
        sets.iter().map(|(_, given)| given.len()).sum::<usize>()
    );
    let mut engine = RecoverEngine::new(&args.recover.run, Origin::Ensemble(&args.recover, sets))?;
    if args.recover_copies {
        engine.quality.count_whole();
    }
    engine.reshed = args.reshed;
    engine.silent_offers = args.silent_offers;
    engine.grouped_offers = args.grouped_offers;
    engine.run()
}

fn set_names(directories: &[String]) -> Vec<String> {
    let base = |source: &String| {
        let at = path::Path::new(source);
        match at.is_file() {
            true => crate::bins::stem(at),
            false => at.file_name().map_or_else(
                || source.clone(),
                |name| name.to_string_lossy().into_owned(),
            ),
        }
    };
    let names = directories.iter().map(base).collect::<Vec<_>>();
    match names.iter().collect::<HashSet<_>>().len() == names.len() {
        true => names,
        false => directories.to_vec(),
    }
}

pub(crate) struct RecoverEngine {
    tables: Tables,
    annotator: Annotator,
    quality: ContigMarkers,
    settings: Settings,
    reports: ReportPaths,
    stages: Option<StageReport>,
    given: Option<Placed>,
    candidates: Vec<BinSet>,
    reshed: bool,
    silent_offers: bool,
    grouped_offers: bool,
    prefix: String,
    output_directory: String,
}

impl RecoverEngine {
    pub fn new(args: &RunArgs, origin: Origin) -> Result<Self> {
        let settings = Settings::of(args, origin.partition_seeds(), origin.split_contamination());
        let Inputs {
            mut tables,
            annotator,
            marking,
            given,
            candidates,
        } = read_inputs(args, &settings, &origin)?;
        let contigs = (0..tables.coverage.contig_names.len()).collect::<Vec<usize>>();
        let windows = tables.windows.clone();
        let learn = |tables: &mut Tables, markers: Option<&ContigMarkers>| match windows.as_deref()
        {
            Some(dir) => {
                learned::learn_space(tables, &settings, dir, args, given.is_some(), markers)
            }
            None => Ok(()),
        };
        let mut quality = match settings.marker_steps {
            Some(_) => {
                let quality = joined(marking)?;
                learn(&mut tables, Some(&quality))?;
                quality
            }
            None => {
                learn(&mut tables, None)?;
                joined(marking)?
            }
        };
        if !candidates.is_empty() {
            quality.count_in_place();
        }
        quality.count_sparse_as_presence();
        let mut engine = Self {
            tables,
            annotator,
            quality,
            settings,
            reports: args.reports.clone(),
            stages: args
                .reports
                .stages
                .as_deref()
                .map(StageReport::create)
                .transpose()?,
            given,
            candidates,
            reshed: false,
            silent_offers: false,
            grouped_offers: false,
            prefix: origin.prefix(),
            output_directory: args.common.output_directory.clone(),
        };
        // A table gives every view one mean, so a k-mer pass's near complete bins teach the space
        // as it first shipped. Those bins need the markers.
        if windows.is_none() {
            engine.learn_labelled(&contigs)?;
        }
        Ok(engine)
    }

    pub fn run(mut self) -> Result<()> {
        let _gpu = crate::learned::gemm::gpu(self.settings.gpu);
        let all_contigs = (0..self.tables.coverage.contig_names.len()).collect::<Vec<usize>>();
        let Settled {
            graph,
            knn,
            mut partitioning,
            weight,
            parked,
            worth,
        } = match self.given.as_ref() {
            Some(given) => self.given_partition(&all_contigs, &given.bins)?,
            None => self.partitioned(&all_contigs)?,
        };
        self.settings.features.distance.aggregate_weight = weight;
        self.trace_partition("stage1_settled", &partitioning, &all_contigs);
        self.settings.bars.genome = self.smallest_genome(&partitioning, &all_contigs);
        debug!(
            "Smallest genome: {} bp (smallest piece {} bp)",
            self.settings.bars.genome,
            self.settings.bars.piece()
        );
        debug!("Partition score {:?}", partitioning.score);
        debug!(
            "Outlier percentage: {:.2}",
            outlier_percentage(&partitioning)
        );

        self.evaluate_outliers(&mut partitioning, &knn)?;
        debug!(
            "Outlier percentage: {:.2}",
            outlier_percentage(&partitioning)
        );
        self.trace_partition("stage2_outliers", &partitioning, &all_contigs);

        let (mut refined, finished) = self.refine_clusters(partitioning, &graph, &knn);
        if !self.candidates.is_empty() {
            self.revote(&mut refined, &all_contigs);
            if self.reshed {
                self.finish_offered(&mut refined, finished);
            }
        }
        let mut bins = refined.bins;
        let mut outliers = refined
            .unbinned
            .into_iter()
            .chain(parked.iter().copied())
            .collect::<HashSet<_>>();
        let floor = self.attach(&mut bins, &mut outliers, parked, worth, &knn)?;
        self.trace_bins("stage6_attached", &bins);

        conserved(
            bins.values()
                .flatten()
                .copied()
                .chain(outliers.iter().copied()),
            &all_contigs.iter().copied().collect(),
        )?;
        let published = self.publish(bins, outliers, floor, &knn)?;
        let scoring = Scoring::of(&published, &self.prefix);

        let written = {
            let _timer = crate::timing::scope("write");
            self.write_clusters(&published, floor)?
        };
        self.write_tables(&scoring, &written)?;
        if !self.candidates.is_empty()
            && let Err(error) = self.write_sources(&scoring)
        {
            warn!("Bin sources: {error}");
        }
        if let Err(error) = self.write_origins(&published) {
            warn!("Replicon origins: {error}");
        }
        if let Err(error) = self.write_strays(&published) {
            warn!("Bridge strays: {error}");
        }

        crate::timing::report(
            path::Path::new(&self.output_directory).join(crate::timing::TIMINGS_FILE),
        )?;

        Ok(())
    }

    fn features(&self) -> ContigFeatures<'_> {
        self.tables.features(&self.settings.features)
    }

    fn lengths(&self) -> &[usize] {
        &self.tables.coverage.contig_lengths
    }

    fn trace_partition(&self, name: &str, held: &Partitioning, contigs: &[usize]) {
        if let Some(report) = &self.stages {
            report.partitioning(name, held, contigs, &self.tables.coverage.contig_names);
        }
    }

    fn trace_bins(&self, name: &str, bins: &std::collections::BTreeMap<usize, Vec<usize>>) {
        if let Some(report) = &self.stages {
            report.bins(name, bins, &self.tables.coverage.contig_names);
        }
    }

    // Rows run long first, so the long contigs are a prefix of any row-ordered list.
    fn long_end(&self, contigs: &[usize]) -> usize {
        let cutoff = self.settings.lengths.cutoff;
        contigs.partition_point(|contig| self.lengths()[*contig] >= cutoff)
    }
}

fn joined(marking: JoinHandle<Result<ContigMarkers>>) -> Result<ContigMarkers> {
    marking
        .join()
        .map_err(|_| anyhow::anyhow!("the marker annotation panicked"))?
}

fn outlier_percentage(partitioning: &Partitioning) -> f64 {
    let outliers = partitioning.outliers.len();
    let binned = partitioning
        .cluster_map
        .values()
        .map(Vec::len)
        .sum::<usize>();
    100.0 * outliers as f64 / (binned + outliers).max(1) as f64
}
