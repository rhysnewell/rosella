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
    recover::bin_writer::{Published, REPLICON_PREFIX},
    recover::census::{Census, STAGES_FILE},
    recover::inputs::{Inputs, read_inputs},
    recover::ladder::{Judge, best_per_arm, combine},
    recover::partition_report::PartitionReport,
    recover::settings::seeds,
    refine::{
        dissolve::{PoolView, RoundParams},
        splitter::{RefineSettings, Refiner},
    },
    seeds::Seeds,
};

mod attach;
mod attach_report;
mod cuts;
mod stages;
mod weight;

pub use stages::{SHIPPED_ORDER, Stage, parse_order, stage_label};

pub const UNBINNED: &str = "unbinned";

/// The fuzzy set needs two neighbours, and a subset of three is the smallest that has them.
pub(crate) const MIN_RESCUE_CONTIGS: usize = 3;

pub fn run_recover(args: &RecoverArgs) -> Result<()> {
    if let Some(path) = &args.reports.floor_report {
        return crate::recover::floor_report::write(args, path);
    }
    RecoverEngine::new(args)?.run()
}

pub(crate) struct RecoverEngine {
    pub(crate) output_directory: String,
    pub(crate) assembly: String,
    pub(crate) coverage_table: CoverageTable,
    coverage_file: String,
    pub(crate) tnf_table: KmerFrequencyTable,
    pub(crate) n_neighbours: usize,
    pub(crate) knn_candidates: usize,
    pub(crate) seeds: Seeds,
    pub(crate) n_contigs: usize,
    pub(crate) min_bin_size: usize,
    pub(crate) min_contig_size: usize,
    cutoff: usize,
    attach_given: bool,
    parked: Vec<usize>,
    worth_spread: f64,
    kept: f64,
    annotator: crate::markers::Annotator,
    pub(crate) max_bin_size: usize,
    pub(crate) max_retries: usize,
    worth: f64,
    links: Option<Vec<crate::assembly_graph::Link>>,
    link_weight: f32,
    sketches: Option<ContigSketches>,
    pub(crate) distance: DistanceSettings,
    dissolve: bool,
    dissolve_rounds: usize,
    dissolve_passes: usize,
    partition_seeds: usize,
    join: bool,
    min_completeness: f64,
    contamination_bar: f64,
    pub(crate) quality: crate::markers::ContigMarkers,
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
            marker_report: args
                .reports
                .marker_report
                .as_ref()
                .map(std::path::PathBuf::from),
            attach_report: args.reports.attach_report.clone(),
            partition_report: args.reports.partition_report.clone(),
            reach_report: args.reports.reach_report.clone(),
            audit_report: args.reports.audit_report.clone(),
            shed_report: args.reports.shed_report.clone(),
            cut_report: args.reports.cut_report.clone(),
            pool_report: args
                .reports
                .pool_report
                .as_ref()
                .map(std::path::PathBuf::from),
            combine_report: args.reports.combine_report.clone(),
        })
    }

    pub fn run(mut self) -> Result<()> {
        let all_contigs = (0..self.n_contigs).collect::<Vec<usize>>();

        debug!("Embedding.");
        let (graph, knn, mut partitioning) = self.partitioned(&all_contigs)?;
        let induced = &knn;
        if self.knn_report.is_some() || self.reach_report.is_some() {
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
            "Outlier percentage: {}",
            partitioning.outliers.len() as f64 / self.n_contigs as f64
        );

        let mut census = Census::default();
        self.census_of(&mut census, "partition", &partitioning);

        debug!("Rescuing unbinned.");
        self.evaluate_outliers(&mut partitioning, induced)?;
        debug!(
            "Outlier percentage: {}",
            partitioning.outliers.len() as f64 / self.n_contigs as f64
        );
        self.census_of(&mut census, "outlier_pool", &partitioning);

        if self.max_retries > 0 {
            debug!("Refining bins.");
        }
        let mut cuts = None;
        let (mut cluster_map, mut outliers) =
            self.refine_clusters(partitioning, &graph, induced, &mut census, &mut cuts);
        outliers.extend(self.parked.iter().copied());
        self.attach(&mut cluster_map, &mut outliers, induced)?;
        if let Some(path) = &self.marker_report {
            self.quality
                .report(&self.coverage_table.contig_names, path)?;
        }

        conserved(
            cluster_map
                .values()
                .flatten()
                .copied()
                .chain(outliers.iter().copied()),
            &all_contigs.iter().copied().collect(),
        )?;
        let published = self.publish_traced(cluster_map, outliers, cuts);
        self.write_quality(&published);
        let cluster_results = self.get_cluster_result(published);
        debug!("Length of cluster results: {}", cluster_results.len());

        debug!("Writing clusters.");
        {
            let _timer = crate::timing::scope("write");
            self.write_clusters(&cluster_results)?;
        }

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
                self.seeds.partition + step as u64,
            )?);
        }
        Ok(self.pick_partition(ladder, contigs, report))
    }

    /// Written from the bins that are written out, not from the refiner's last pass, so the
    /// table and the assignments never describe different partitions.
    fn write_quality(&mut self, published: &Published) {
        let binned = published
            .bins
            .values()
            .flatten()
            .copied()
            .collect::<Vec<_>>();
        let names = &self.coverage_table.contig_names;
        if let Err(error) = self
            .annotator
            .complete_checkm(&mut self.quality, &binned, names)
        {
            warn!("Could not search the CheckM models: {error}");
        }
        let mut sorted = published
            .bins
            .iter()
            .map(|(bin, contigs)| {
                let mut contigs = contigs.iter().copied().collect::<Vec<_>>();
                contigs.sort_unstable();
                (*bin, contigs)
            })
            .collect::<Vec<_>>();
        sorted.sort_unstable_by_key(|(bin, _)| *bin);
        let report = crate::quality::write_report(
            &self.quality,
            sorted
                .iter()
                .map(|(bin, contigs)| (format!("rosella_bin_{bin}"), contigs.as_slice()))
                .chain(published.replicons.iter().enumerate().map(|(at, contig)| {
                    (
                        format!("rosella_bin_{REPLICON_PREFIX}{}", at + 1),
                        std::slice::from_ref(contig),
                    )
                })),
            &self.coverage_table.contig_lengths,
            &path::Path::new(&self.output_directory).join(crate::defaults::QUALITY_FILE),
        );
        if let Err(error) = report {
            warn!("Could not write the quality table: {error}");
        }
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
        let outliers = std::mem::take(&mut partitioning.outliers);
        if outliers.len() < MIN_RESCUE_CONTIGS {
            partitioning.outliers = outliers;
            return Ok(());
        }
        let (knn, order) =
            self.pool_neighbours(&outliers, self.n_neighbours, PoolView::Combined, induced)?;
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
    ) -> (HashMap<usize, HashSet<usize>>, HashSet<usize>) {
        let bins = partitioning
            .cluster_map
            .into_iter()
            .map(|(bin_id, contigs)| {
                let mut contigs = contigs.into_iter().collect::<Vec<_>>();
                contigs.sort_unstable();
                (bin_id, contigs)
            })
            .collect::<BTreeMap<_, _>>();
        let mut unbinned = partitioning.outliers.into_iter().collect::<Vec<_>>();
        unbinned.sort_unstable();

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
        let mut refiner = Refiner::new(self.features(), settings, bins, unbinned)
            .with_assembly(assembly)
            .with_quality(&self.quality)
            .with_cuts(self.cut_report.is_some());
        refiner.run();
        self.census_bins(census, "refine", &refiner.bins, &refiner.unbinned);

        let bars = self.bars();

        let mut seen: std::collections::HashMap<Stage, usize> = std::collections::HashMap::new();
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
        let cluster_map = refiner
            .bins
            .iter()
            .map(|(bin_id, contigs)| (*bin_id, contigs.iter().copied().collect::<HashSet<_>>()))
            .collect::<HashMap<_, _>>();
        (cluster_map, refiner.unbinned.iter().copied().collect())
    }

    fn bars(&self) -> crate::refine::rung::Bars {
        crate::refine::rung::Bars {
            min_bin_size: self.min_bin_size,
            completeness: self.min_completeness,
            contamination: self.contamination_bar,
            worth: self.worth,
            rung_floor: crate::refine::rung::DEFAULT_RUNG_FLOOR,
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

    /// Partition a subset of contigs. `contigs` are indices into the contig list as it
    /// stands after the initial length filter.
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
            self.seeds,
            self.knn_candidates,
            KNN_ASSEMBLY,
        );
        let graph = features.graph_from_knn(contigs, &built);
        (graph, built)
    }

    fn write_knn_report(&self, contigs: &[usize], path: &path::Path) -> Result<()> {
        let combined = self.features().knn_of(
            contigs,
            self.n_neighbours,
            self.seeds,
            self.knn_candidates,
            KNN_ASSEMBLY,
        );
        let rho = self
            .features()
            .with_distance(self.distance.composition_only())
            .knn_of(
                contigs,
                self.n_neighbours,
                self.seeds,
                self.knn_candidates,
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
        contig_indices: &HashSet<usize>,
        n_neighbours: usize,
        view: PoolView,
        induced: &KnnGraph,
    ) -> Result<(KnnGraph, Vec<usize>)> {
        let mut order = contig_indices.iter().copied().collect::<Vec<_>>();
        order.sort_unstable();
        if view == PoolView::Combined
            && let Some(built) = induced.induced(&order)
        {
            return Ok((built, order));
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
            self.seeds,
            self.knn_candidates,
            KNN_POOL,
        );
        Ok((knn, order))
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
