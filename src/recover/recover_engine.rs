use std::{
    collections::{BTreeMap, HashMap, HashSet},
    path,
};

use anyhow::Result;
use log::{debug, warn};

use crate::{
    cli::RecoverArgs,
    clustering::{
        clusterer::{Partitioning, conserved, find_partitions},
        graph_partition::Partition,
    },
    coverage::coverage_table::CoverageTable,
    embedding::{
        KNN_ASSEMBLY, KNN_POOL, features::ContigFeatures, knn::KnnGraph, metrics::DistanceSettings,
    },
    kmers::kmer_counting::KmerFrequencyTable,
    kmers::sketch::ContigSketches,
    recover::census::{Census, STAGES_FILE},
    recover::inputs::{Inputs, read_inputs},
    recover::ladder::{Judge, best_per_arm, combine},
    recover::partition_report::PartitionReport,
    recover::settings::seeds,
    refine::{
        dissolve::{MIN_RESCUE_CONTIGS, PoolView, RoundParams},
        splitter::{RefineSettings, Refiner},
    },
    seeds::Seeds,
};

mod attach;
mod attach_report;
mod bin_writer;
mod cuts;
mod stages;
mod tables;
mod weight;

pub use stages::{SHIPPED_ORDER, Stage, parse_order, stage_label};
use tables::Scoring;

pub const UNBINNED: &str = "unbinned";

pub fn run_recover(args: &RecoverArgs) -> Result<()> {
    if let Some(path) = &args.reports.floor_report {
        return crate::recover::floor_report::write(args, path);
    }
    RecoverEngine::new(args)?.run()
}

pub(crate) struct RecoverEngine {
    output_directory: String,
    assembly: String,
    coverage_table: CoverageTable,
    coverage_file: String,
    tnf_table: KmerFrequencyTable,
    n_neighbours: usize,
    knn_candidates: usize,
    seeds: Seeds,
    n_contigs: usize,
    min_bin_size: usize,
    min_contig_size: usize,
    cutoff: usize,
    attach_given: bool,
    parked: Vec<usize>,
    worth_spread: f64,
    kept: f64,
    annotator: crate::markers::Annotator,
    max_bin_size: usize,
    max_retries: usize,
    worth: f64,
    links: Option<Vec<crate::assembly_graph::Link>>,
    link_weight: f32,
    sketches: Option<ContigSketches>,
    distance: DistanceSettings,
    dissolve: bool,
    dissolve_rounds: usize,
    dissolve_passes: usize,
    partition_seeds: usize,
    join: bool,
    min_completeness: f64,
    contamination_bar: f64,
    quality: crate::markers::ContigMarkers,
    oracle: Vec<Vec<usize>>,
    partition: Partition,
    stage_order: Vec<Stage>,
    knn_report: Option<std::path::PathBuf>,
    marker_report: Option<std::path::PathBuf>,
    attach_report: Option<std::path::PathBuf>,
    partition_report: Option<std::path::PathBuf>,
    reach_report: Option<Vec<std::path::PathBuf>>,
    audit_report: Option<std::path::PathBuf>,
    shed_report: Option<std::path::PathBuf>,
    cut_report: Option<std::path::PathBuf>,
    pool_report: Option<std::path::PathBuf>,
    combine_report: Option<std::path::PathBuf>,
    halves: crate::kmers::halves::Halves,
}

impl RecoverEngine {
    pub fn new(args: &RecoverArgs) -> Result<Self> {
        let Inputs {
            output_directory,
            assembly,
            min_contig_size,
            coverage_table,
            coverage_file,
            tnf_table,
            sketches,
            halves,
            links,
            quality,
            oracle,
            distance,
            partition,
            dissolve,
            cutoff,
            attach_given,
            annotator,
        } = read_inputs(args)?;

        let n_neighbours = args.graph.n_neighbours;
        let knn_candidates = args.graph.candidates();
        let seeds = seeds(&args.seeds);
        let min_bin_size = args.binning.min_bin_size;

        let n_contigs = coverage_table.contig_names.len();
        let max_bin_size = args.binning.max_bin_size;
        let max_retries = usize::from(!args.no_refine);

        Ok(Self {
            output_directory,
            assembly,
            coverage_table,
            coverage_file,
            tnf_table,
            n_neighbours,
            knn_candidates,
            seeds,
            n_contigs,
            min_bin_size,
            min_contig_size,
            cutoff,
            attach_given,
            parked: Vec::new(),
            worth_spread: 0.0,
            kept: 0.0,
            annotator,
            max_bin_size,
            max_retries,
            worth: args.rescue.worth_contamination,
            links,
            link_weight: args.graph.assembly_graph_weight as f32,
            sketches,
            distance,
            dissolve,
            dissolve_rounds: args.rescue.dissolve_rounds as usize,
            dissolve_passes: args.rescue.dissolve_passes as usize,
            partition_seeds: args.rescue.partition_seeds as usize,
            join: args.join,
            min_completeness: args.rescue.min_completeness,
            contamination_bar: args.rescue.max_contamination,
            quality,
            oracle,
            partition,
            stage_order: parse_order(&args.rescue.stage_order)?,
            knn_report: args.reports.knn_report.clone(),
            marker_report: args.reports.marker_report.clone(),
            attach_report: args.reports.attach_report.clone(),
            partition_report: args.reports.partition_report.clone(),
            reach_report: args.reports.reach_report.clone(),
            audit_report: args.reports.audit_report.clone(),
            shed_report: args.reports.shed_report.clone(),
            cut_report: args.reports.cut_report.clone(),
            pool_report: args.reports.pool_report.clone(),
            combine_report: args.reports.combine_report.clone(),
            halves,
        })
    }

    pub fn run(mut self) -> Result<()> {
        let all_contigs = (0..self.n_contigs).collect::<Vec<usize>>();

        debug!("Embedding.");
        let (graph, knn, mut partitioning) = self.partitioned(&all_contigs)?;
        let induced = &knn;
        if self.knn_report.is_some() {
            self.load(self.n_contigs)?;
        }

        if let Some(path) = &self.knn_report {
            self.write_knn_report(&all_contigs, path)?;
            debug!("Wrote the kNN report to {}.", path.display());
            return Ok(());
        }

        if let Some(paths) = &self.reach_report {
            let groups = crate::refine::oracle::read_groups(
                &paths[0].to_string_lossy(),
                &self.coverage_table.contig_names,
            )?;
            crate::embedding::reach::write(
                &paths[1],
                &groups,
                &self.features(),
                &knn,
                &self.coverage_table.contig_lengths,
                &self.coverage_table.contig_names,
            )?;
            debug!("Wrote the reach report to {}.", paths[1].display());
            return Ok(());
        }

        if self.partition_report.is_some() {
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
        self.evaluate_outliers(&mut partitioning, induced)?;
        debug!(
            "Outlier percentage: {:.2}",
            outlier_percentage(&partitioning)
        );
        self.census_of(&mut census, "outlier_pool", &partitioning);

        if self.max_retries > 0 {
            debug!("Refining bins.");
        }
        let mut cuts = None;
        let (mut bins, mut outliers) =
            self.refine_clusters(partitioning, &graph, induced, &mut census, &mut cuts);
        outliers.extend(self.parked.iter().copied());
        self.attach(&mut bins, &mut outliers, induced)?;

        conserved(
            bins.values()
                .flatten()
                .copied()
                .chain(outliers.iter().copied()),
            &all_contigs.iter().copied().collect(),
        )?;
        let published = self.publish_traced(bins, outliers, cuts);
        let scoring = Scoring::of(&published);

        debug!("Writing clusters.");
        let written = {
            let _timer = crate::timing::scope("write");
            self.write_clusters(&published)?
        };
        self.write_tables(&scoring, &written)?;

        crate::timing::report(
            path::Path::new(&self.output_directory).join(crate::timing::TIMINGS_FILE),
        )?;
        census.write(path::Path::new(&self.output_directory).join(STAGES_FILE))?;

        Ok(())
    }

    fn partition_all(
        &self,
        graph: &crate::embedding::Graph,
        contigs: &[usize],
        report: Option<&mut PartitionReport>,
    ) -> Result<Partitioning> {
        let mut ladder = Vec::new();
        for step in 0..self.partition_seeds {
            ladder.extend(self.partition_of(
                graph,
                contigs,
                self.partition,
                true,
                self.seeds.partition.wrapping_add(step as u64),
            )?);
        }
        Ok(self.pick_partition(ladder, contigs, report))
    }

    fn census_of(&self, census: &mut Census, stage: &str, result: &Partitioning) {
        census.record(
            stage,
            result
                .cluster_map
                .values()
                .map(|contigs| contigs.iter().copied()),
            result.outliers.iter().copied(),
            &self.coverage_table.contig_lengths,
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
            &self.coverage_table.contig_lengths,
        );
    }

    fn evaluate_outliers(&self, partitioning: &mut Partitioning, induced: &KnnGraph) -> Result<()> {
        if partitioning.outliers.len() < MIN_RESCUE_CONTIGS {
            return Ok(());
        }
        let outliers = std::mem::take(&mut partitioning.outliers);
        let (knn, order) =
            self.pool_neighbours(outliers, self.n_neighbours, PoolView::Combined, induced);
        let partitioning_of_filtered_contigs = self
            .evaluate_subset(
                &knn,
                &order,
                RoundParams {
                    n_neighbours: self.n_neighbours,
                    ladder: false,
                },
            )?
            .swap_remove(0);

        debug!(
            "New Partition score {:?}",
            partitioning_of_filtered_contigs.score
        );
        debug!(
            "Number of clusters: {}",
            partitioning_of_filtered_contigs.cluster_map.len()
        );
        partitioning.merge(partitioning_of_filtered_contigs);

        Ok(())
    }

    fn refine_clusters(
        &self,
        partitioning: Partitioning,
        assembly: &crate::embedding::Graph,
        induced: &KnnGraph,
        census: &mut Census,
        cuts: &mut Option<crate::refine::cut_report::CutLog>,
    ) -> (BTreeMap<usize, Vec<usize>>, HashSet<usize>) {
        let settings = RefineSettings {
            min_bin_size: self.min_bin_size,
            max_bin_size: self.max_bin_size,
            n_neighbours: self.n_neighbours,
            knn_candidates: self.knn_candidates,
            max_retries: self.max_retries,
            seeds: self.seeds,
            max_contamination: None,
            partition: self.partition,
        };
        let mut refiner = Refiner::new(
            self.features(),
            settings,
            partitioning.cluster_map,
            partitioning.outliers,
        )
        .with_assembly(assembly)
        .with_quality(&self.quality)
        .with_cuts(self.cut_report.is_some());
        refiner.run();
        self.census_bins(census, "refine", &refiner.bins, &refiner.unbinned);

        let bars = self.bars();

        let mut seen: HashMap<Stage, usize> = HashMap::new();
        let mut cycle = stages::Cycle {
            refiner: &mut refiner,
            census,
            induced,
            bars,
            finished: crate::refine::finished::Finished::default(),
        };
        for stage in &self.stage_order {
            let pass = seen.entry(*stage).or_insert(0);
            self.run_stage(*stage, &mut cycle, *pass);
            *pass += 1;
        }

        *cuts = refiner.cuts.take();
        let mut bins = std::mem::take(&mut refiner.bins);
        for contigs in bins.values_mut() {
            contigs.sort_unstable();
        }
        (bins, refiner.unbinned.iter().copied().collect())
    }

    fn bars(&self) -> crate::refine::rung::Bars {
        crate::refine::rung::Bars {
            min_bin_size: self.min_bin_size,
            completeness: self.min_completeness,
            contamination: self.contamination_bar,
            worth: self.worth,
        }
    }

    fn pick_partition(
        &self,
        ladder: Vec<Partitioning>,
        contigs: &[usize],
        mut partitions: Option<&mut PartitionReport>,
    ) -> Partitioning {
        let judge = Judge {
            quality: &self.quality,
            contigs,
            bars: self.bars(),
        };
        let report = self.combine_report.as_ref().and_then(|path| {
            crate::recover::combine_report::CombineReport::create(
                path,
                &self.coverage_table.contig_names,
                &self.coverage_table.contig_lengths,
            )
            .map_err(|error| warn!("Could not write the combine report: {error}"))
            .ok()
        });
        if let Some(partitions) = partitions.as_deref_mut() {
            partitions.ladder(&ladder);
        }
        let arms = best_per_arm(ladder, &judge);
        let chosen = combine(&arms, &judge, report.as_ref());
        if let Some(report) = &report {
            report.flush();
        }
        if let Some(partitions) = partitions {
            partitions.chosen(&arms);
            partitions.add("combined", &chosen);
        }
        chosen
    }

    // `contigs` index the contig list as it stands after the initial length filter.
    fn partition_of(
        &self,
        graph: &crate::embedding::Graph,
        contigs: &[usize],
        kind: Partition,
        rank_rungs: bool,
        partition_seed: u64,
    ) -> Result<Vec<Partitioning>> {
        find_partitions(
            graph,
            &self.features().contig_lengths(contigs),
            partition_seed,
            kind,
            rank_rungs,
        )
    }

    fn embed(&self, contigs: &[usize]) -> (crate::embedding::Graph, KnnGraph) {
        let features = self.features();
        let built = features.knn_of(
            contigs,
            self.n_neighbours,
            self.knn_candidates,
            self.seeds.knn,
            KNN_ASSEMBLY,
        );
        let graph = features.graph_from_knn(contigs, &built);
        (graph, built)
    }

    fn write_knn_report(&self, contigs: &[usize], path: &path::Path) -> Result<()> {
        let combined = self.features().knn_of(
            contigs,
            self.n_neighbours,
            self.knn_candidates,
            self.seeds.knn,
            KNN_ASSEMBLY,
        );
        let rho = self
            .features()
            .with_distance(self.distance.composition_only())
            .knn_of(
                contigs,
                self.n_neighbours,
                self.knn_candidates,
                self.seeds.knn,
                KNN_ASSEMBLY,
            );
        let names = contigs
            .iter()
            .map(|index| self.coverage_table.contig_names[*index].as_str())
            .collect::<Vec<_>>();
        crate::embedding::knn::write_report(&[("combined", combined), ("rho", rho)], &names, path)
    }

    fn pool_neighbours(
        &self,
        order: Vec<usize>,
        n_neighbours: usize,
        view: PoolView,
        induced: &KnnGraph,
    ) -> (KnnGraph, Vec<usize>) {
        if view == PoolView::Combined
            && let Some(built) = induced.induced(&order)
        {
            return (built, order);
        }
        let features = match view {
            PoolView::Combined => self.features(),
            PoolView::Composition => self
                .features()
                .with_distance(self.distance.composition_only()),
        };
        let knn = features.knn_of(
            &order,
            n_neighbours,
            self.knn_candidates,
            self.seeds.knn,
            KNN_POOL,
        );
        (knn, order)
    }

    fn evaluate_subset(
        &self,
        knn: &KnnGraph,
        ordered_indices: &[usize],
        round: RoundParams,
    ) -> Result<Vec<Partitioning>> {
        let subset_graph = self.features().graph_from_knn(ordered_indices, knn);
        // Label propagation returns one labelling, and the rungs need a ladder to walk.
        let kind = match round.ladder && !self.partition.runs_leiden() {
            true => {
                debug!("Pool rungs need a ladder, so this round runs Leiden");
                Partition::Leiden
            }
            false => self.partition,
        };
        let mut results = self.partition_of(
            &subset_graph,
            ordered_indices,
            kind,
            !round.ladder,
            self.seeds.partition,
        )?;
        if !round.ladder {
            results.truncate(1);
        }
        debug!("Partition score {:?}", results[0].score);

        let expected = ordered_indices.iter().copied().collect();
        for result in results.iter_mut() {
            result.reindex_clusters(ordered_indices);
            conserved(
                result
                    .cluster_map
                    .values()
                    .flatten()
                    .copied()
                    .chain(result.outliers.iter().copied()),
                &expected,
            )?;
        }

        Ok(results)
    }

    fn features(&self) -> ContigFeatures<'_> {
        ContigFeatures::new(
            &self.coverage_table.table,
            &self.tnf_table.kmer_table,
            &self.coverage_table.contig_lengths,
        )
        .with_distance(self.distance)
        .with_links(self.links.as_deref(), self.link_weight)
        .with_sketches(self.sketches.as_ref())
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
