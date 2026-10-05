use std::{
    collections::{BTreeMap, HashSet},
    path,
    sync::Arc,
};

use anyhow::Result;
use log::debug;

use crate::{
    cli::{RecoverArgs, ReportPaths},
    clustering::clusterer::{Partitioning, conserved},
    embedding::features::ContigFeatures,
    markers::{Annotator, ContigMarkers},
    recover::census::{Census, STAGES_FILE},
    recover::inputs::{Inputs, read_inputs},
    recover::settings::Settings,
    tables::Tables,
};

mod attach;
mod attach_report;
mod bin_writer;
mod cuts;
mod evidence;
mod partition;
mod probes;
mod stages;
mod tables;
mod weight;

pub use stages::{SHIPPED_ORDER, Stage, parse_order, stage_label};
use tables::Scoring;
use weight::Settled;

pub const UNBINNED: &str = "unbinned";

pub fn run_recover(args: &RecoverArgs) -> Result<()> {
    if let Some(path) = &args.reports.floor_report {
        return crate::recover::floor_report::write(args, path);
    }
    RecoverEngine::new(args)?.run()
}

pub(crate) struct RecoverEngine {
    tables: Tables,
    annotator: Annotator,
    quality: ContigMarkers,
    settings: Settings,
    reports: ReportPaths,
    oracle: Vec<Vec<usize>>,
    output_directory: String,
}

impl RecoverEngine {
    pub fn new(args: &RecoverArgs) -> Result<Self> {
        let settings = Settings::of(args)?;
        let Inputs {
            tables,
            annotator,
            quality,
            oracle,
        } = read_inputs(args, &settings)?;
        Ok(Self {
            tables,
            annotator,
            quality,
            settings,
            reports: args.reports.clone(),
            oracle,
            output_directory: args.common.output_directory.clone(),
        })
    }

    pub fn run(mut self) -> Result<()> {
        let all_contigs = (0..self.tables.coverage.contig_names.len()).collect::<Vec<usize>>();
        self.settings.features.distance.evidence = self.evidence(&all_contigs)?.map(Arc::new);

        debug!("Embedding.");
        let Settled {
            graph,
            knn,
            mut partitioning,
            weight,
            unit,
            parked,
            worth,
        } = self.partitioned(&all_contigs)?;
        self.settings.features.distance.aggregate_weight = weight;
        self.settings.features.distance.variance_unit = unit;
        if self.probed(&all_contigs, &knn)? {
            return Ok(());
        }

        debug!("Partition score {:?}", partitioning.score);
        debug!(
            "Outlier percentage: {:.2}",
            outlier_percentage(&partitioning)
        );

        let mut census = Census::default();
        self.census_of(&mut census, "partition", &partitioning);

        debug!("Rescuing unbinned.");
        self.evaluate_outliers(&mut partitioning, &knn)?;
        debug!(
            "Outlier percentage: {:.2}",
            outlier_percentage(&partitioning)
        );
        self.census_of(&mut census, "outlier_pool", &partitioning);

        if self.settings.cycle.refine_rounds > 0 {
            debug!("Refining bins.");
        }
        let refined = self.refine_clusters(partitioning, &graph, &knn, &mut census);
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
        let published = self.publish_traced(bins, outliers, refined.cuts);
        let scoring = Scoring::of(&published);

        debug!("Writing clusters.");
        let written = {
            let _timer = crate::timing::scope("write");
            self.write_clusters(&published, floor)?
        };
        self.write_tables(&scoring, &written)?;

        crate::timing::report(
            path::Path::new(&self.output_directory).join(crate::timing::TIMINGS_FILE),
        )?;
        census.write(path::Path::new(&self.output_directory).join(STAGES_FILE))?;

        Ok(())
    }

    fn census_of(&self, census: &mut Census, stage: &str, result: &Partitioning) {
        census.record(
            stage,
            result
                .cluster_map
                .values()
                .map(|contigs| contigs.iter().copied()),
            result.outliers.iter().copied(),
            &self.tables.coverage.contig_lengths,
        );
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
            &self.tables.coverage.contig_lengths,
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
