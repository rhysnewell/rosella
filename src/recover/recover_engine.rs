use std::{
    collections::{BTreeMap, HashMap, HashSet},
    path,
};

use anyhow::Result;
use log::{debug, warn};

use crate::{
    cli::{EnsembleArgs, PartitionProbes, RecoverArgs, ReportPaths, RunArgs},
    clustering::clusterer::{Partitioning, conserved},
    embedding::features::ContigFeatures,
    markers::{Annotator, ContigMarkers},
    recover::census::{Census, STAGES_FILE},
    recover::inputs::{Inputs, read_inputs},
    recover::origin::{BinSet, Given, Origin, Placed},
    recover::settings::Settings,
    tables::Tables,
};

mod attach;
mod attach_report;
mod bin_writer;
mod cuts;
mod learned;
mod offered;
mod partition;
mod probes;
mod sources;
mod stages;
mod tables;
mod weight;

pub use stages::{SHIPPED_ORDER, Stage, parse_order, stage_label};
use tables::Scoring;
use weight::Settled;

pub const UNBINNED: &str = "unbinned";

pub fn run_recover(args: &RecoverArgs) -> Result<()> {
    if let Some(path) = &args.probes.floor_report {
        return crate::recover::floor_report::write(&args.run, path);
    }
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
    RecoverEngine::new(&args.recover.run, Origin::Ensemble(&args.recover, sets))?.run()
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
    probes: PartitionProbes,
    oracle: Vec<Vec<usize>>,
    given: Option<Placed>,
    candidates: Vec<BinSet>,
    prefix: String,
    output_directory: String,
}

impl RecoverEngine {
    pub fn new(args: &RunArgs, origin: Origin) -> Result<Self> {
        let settings = Settings::of(args, origin.partition_seeds(), origin.split_contamination())?;
        let Inputs {
            mut tables,
            annotator,
            marking,
            oracle,
            given,
            candidates,
        } = read_inputs(args, &settings, &origin)?;
        let contigs = (0..tables.coverage.contig_names.len()).collect::<Vec<usize>>();
        let windows = tables.windows.clone();
        if let Some(dir) = &windows {
            learned::learn_space(&mut tables, &settings, dir, args, given.is_some())?;
        }
        let mut quality = marking
            .join()
            .map_err(|_| anyhow::anyhow!("the marker annotation panicked"))??;
        if !candidates.is_empty() {
            quality.count_in_place();
        }
        let mut engine = Self {
            tables,
            annotator,
            quality,
            settings,
            reports: args.reports.clone(),
            probes: origin.probes(),
            oracle,
            given,
            candidates,
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
        self.settings.bars.genome = self.smallest_genome(&partitioning, &all_contigs);
        debug!(
            "Smallest genome: {} bp (smallest piece {} bp)",
            self.settings.bars.genome,
            self.settings.bars.piece()
        );
        if self.probed(&all_contigs, &knn)? {
            return Ok(());
        }

        debug!("Partition score {:?}", partitioning.score);
        debug!(
            "Outlier percentage: {:.2}",
            outlier_percentage(&partitioning)
        );

        let mut census = Census::default();
        let start = match self.given {
            Some(_) => "given",
            None => "partition",
        };
        self.census_of(&mut census, start, &partitioning);

        self.evaluate_outliers(&mut partitioning, &knn)?;
        debug!(
            "Outlier percentage: {:.2}",
            outlier_percentage(&partitioning)
        );
        self.census_of(&mut census, "outlier_pool", &partitioning);

        let mut refined = self.refine_clusters(partitioning, &graph, &knn, &mut census);
        if !self.candidates.is_empty() {
            self.revote(&mut refined, &all_contigs);
            self.census_bins(&mut census, "revote", &refined.bins, &refined.unbinned);
        }
        let mut bins = refined.bins;
        let mut outliers = refined
            .unbinned
            .into_iter()
            .chain(parked.iter().copied())
            .collect::<HashSet<_>>();
        let floor = self.attach(&mut bins, &mut outliers, parked, worth, &knn)?;

        conserved(
            bins.values()
                .flatten()
                .copied()
                .chain(outliers.iter().copied()),
            &all_contigs.iter().copied().collect(),
        )?;
        let published = self.publish_traced(bins, outliers, refined.cuts)?;
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
        if self.settings.cycle.replicons
            && let Err(error) = self.write_origins(&published)
        {
            warn!("Replicon origins: {error}");
        }

        crate::timing::report(
            path::Path::new(&self.output_directory).join(crate::timing::TIMINGS_FILE),
        )?;
        census.write(path::Path::new(&self.output_directory).join(STAGES_FILE))?;

        Ok(())
    }

    fn census_of(&self, census: &mut Census, stage: &str, result: &Partitioning) {
        self.census_bins(census, stage, &result.cluster_map, &result.outliers);
    }

    fn census_bins(
        &self,
        census: &mut Census,
        stage: &str,
        bins: &BTreeMap<usize, Vec<usize>>,
        unbinned: &[usize],
    ) {
        census.record(
            stage,
            bins.values().map(|contigs| contigs.iter().copied()),
            unbinned.iter().copied(),
            self.lengths(),
        );
    }

    fn features(&self) -> ContigFeatures<'_> {
        self.tables.features(&self.settings.features)
    }

    fn lengths(&self) -> &[usize] {
        &self.tables.coverage.contig_lengths
    }

    fn names(&self) -> &[String] {
        &self.tables.coverage.contig_names
    }
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
