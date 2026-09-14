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
    quality::Scorer,
    recover::census::{Census, STAGES_FILE},
    recover::inputs::{Inputs, read_inputs},
    recover::ladder::{Judge, best_per_arm, combine},
    recover::settings::seeds,
    refine::{
        dissolve::{PoolView, RoundParams},
        splitter::{RefineSettings, Refiner},
    },
    seeds::Seeds,
};

pub const UNBINNED: &str = "unbinned";

/// umap-rs asks for two neighbours, and a subset of three is the smallest that has them.
const MIN_RESCUE_CONTIGS: usize = 3;

pub fn run_recover(args: RecoverArgs) -> Result<()> {
    RecoverEngine::new(&args)?.run()
}

pub(crate) struct RecoverEngine {
    pub(crate) output_directory: String,
    pub(crate) assembly: String,
    pub(crate) coverage_table: CoverageTable,
    pub(crate) tnf_table: KmerFrequencyTable,
    pub(crate) n_neighbours: usize,
    pub(crate) seeds: Seeds,
    pub(crate) n_contigs: usize,
    pub(crate) min_bin_size: usize,
    pub(crate) min_contig_size: usize,
    pub(crate) max_bin_size: usize,
    pub(crate) max_retries: usize,
    worth: crate::quality::Worth,
    rung_floor: f64,
    links: Option<Vec<(usize, usize)>>,
    link_weight: f32,
    pub(crate) knn_candidates: usize,
    pub(crate) distance: DistanceSettings,
    bisect: bool,
    dissolve: bool,
    dissolve_rounds: usize,
    dissolve_passes: usize,
    recruit_near_bar: Option<f64>,
    partition_seeds: usize,
    join: bool,
    min_completeness: f64,
    contamination_bar: f64,
    quality: crate::markers::ContigMarkers,
    oracle: Vec<Vec<usize>>,
    level_quantile: f64,
    partition: Partition,
    partition_resolution: Option<f64>,
    partition_theta: Option<f64>,
    knn_report: Option<std::path::PathBuf>,
    pool_report: Option<std::path::PathBuf>,
}

impl RecoverEngine {
    pub fn new(args: &RecoverArgs) -> Result<Self> {
        let Inputs {
            output_directory,
            assembly,
            min_contig_size,
            coverage_table,
            tnf_table,
            links,
            quality,
            oracle,
            distance,
            partition,
            dissolve,
        } = read_inputs(args)?;

        let n_neighbours = args.graph.n_neighbours;
        let seeds = seeds(&args.seeds);
        let min_bin_size = args.binning.min_bin_size;

        let n_contigs = coverage_table.table.nrows();
        let max_bin_size = args.binning.max_bin_size;
        let max_retries = if args.no_refine {
            0
        } else {
            args.refine.max_retries
        };

        Ok(Self {
            output_directory,
            assembly,
            coverage_table,
            tnf_table,
            n_neighbours,
            seeds,
            n_contigs,
            min_bin_size,
            min_contig_size,
            max_bin_size,
            max_retries,
            worth: crate::quality::Worth {
                contamination: args.rescue.worth_contamination,
                allowance: args.rescue.worth_allowance,
            },
            rung_floor: args.rescue.rung_floor,
            links,
            link_weight: args.assembly_graph_weight as f32,
            knn_candidates: args.graph.knn_candidates,
            distance,
            bisect: args.refine.bisect,
            dissolve,
            dissolve_rounds: args.rescue.dissolve_rounds as usize,
            dissolve_passes: args.rescue.dissolve_passes as usize,
            recruit_near_bar: args.rescue.recruit_near_bar,
            partition_seeds: args.rescue.partition_seeds as usize,
            join: !args.no_join,
            min_completeness: args.rescue.min_completeness,
            contamination_bar: args.rescue.max_contamination,
            quality,
            oracle,
            partition,
            partition_resolution: args.binning.partition_resolution,
            partition_theta: args.binning.partition_theta,
            knn_report: args.reports.knn_report.clone(),
            pool_report: args.reports.pool_report.as_ref().map(std::path::PathBuf::from),
            level_quantile: args.refine.split_level_quantile,
        })
    }

    /// Runs through the rosella bin recovery pipeline
    pub fn run(self) -> Result<()> {
        let all_contigs = (0..self.n_contigs).collect::<Vec<usize>>();

        if let Some(path) = &self.knn_report {
            self.write_knn_report(&all_contigs, path)?;
            debug!("Wrote the kNN report to {}.", path.display());
            return Ok(());
        }

        debug!("Embedding.");
        let (graph, knn) = self.embed(&all_contigs);
        let induced = &knn;

        debug!("Clustering.");
        let mut ladder = Vec::new();
        for step in 0..self.partition_seeds {
            ladder.extend(self.partition_of(
                &graph,
                &all_contigs,
                self.partition,
                true,
                self.seeds.partition + step as u64,
            )?);
        }
        let mut partitioning = self.pick_partition(ladder, &all_contigs);
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
        let (cluster_map, outliers) =
            self.refine_clusters(partitioning, &graph, induced, &mut census);

        conserved(
            cluster_map
                .values()
                .flatten()
                .copied()
                .chain(outliers.iter().copied()),
            &all_contigs.iter().copied().collect(),
        )?;
        let cluster_results = self.get_cluster_result(cluster_map, outliers);
        debug!("Length of cluster results: {}", cluster_results.len());

        debug!("Writing clusters.");
        {
            let _timer = crate::timing::scope("write");
            self.write_clusters(cluster_results)?;
        }

        crate::timing::report(
            path::Path::new(&self.output_directory).join(crate::timing::TIMINGS_FILE),
        )?;
        census.write(path::Path::new(&self.output_directory).join(STAGES_FILE))?;

        Ok(())
    }

    fn census_of(&self, census: &mut Census, stage: &'static str, result: &Partitioning) {
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
        stage: &'static str,
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
            max_retries: self.max_retries,
            seeds: self.seeds,
            max_contamination: None,
            bisect: self.bisect,
            level_quantile: self.level_quantile,
            partition: self.partition,
            partition_resolution: self.partition_resolution,
            partition_theta: self.partition_theta,
            knn_candidates: self.knn_candidates,
        };
        let mut refiner = Refiner::new(self.features(), settings, bins, unbinned)
            .with_assembly(assembly);
        refiner.run();
        self.census_bins(census, "refine", &refiner.bins, &refiner.unbinned);

        let bars = self.bars();

        if self.dissolve {
            // Stale by a round, since merge and both eject arms move the bins it was
            // measured on. Recomputing it here was measured and lost bins.
            let settings = crate::refine::dissolve::DissolveSettings {
                bars,
                genome_floor: refiner.genome_floor,
                min_contigs: MIN_RESCUE_CONTIGS,
                rounds: self.dissolve_rounds,
                passes: self.dissolve_passes,
                n_neighbours: self.n_neighbours,
                max_bin_size: self.max_bin_size,
            };
            let report = self.pool_report.as_ref().and_then(|path| {
                crate::refine::pool_report::PoolReport::create(
                    path,
                    &self.coverage_table.contig_names,
                )
                .map_err(|error| warn!("No pool report at {}: {error}", path.display()))
                .ok()
            });
            let ledger = crate::refine::dissolve::dissolve(
                &self.features(),
                &self.quality,
                &mut refiner.bins,
                &mut refiner.unbinned,
                settings,
                &self.oracle,
                report.as_ref(),
                |pool, n_neighbours, view| self.pool_neighbours(pool, n_neighbours, view, induced),
                |knn, order, round| self.evaluate_subset(knn, order, round),
            );
            if let Some(report) = report.as_ref() {
                report.flush();
            }
            debug!("Dissolve pool: {ledger}");
            self.census_bins(census, "dissolve", &refiner.bins, &refiner.unbinned);
        }

        if self.join {
            let _timer = crate::timing::scope("join");
            let ledger = crate::refine::join::join(
                &self.features(),
                &self.quality,
                &mut refiner.bins,
                crate::refine::join::JoinSettings {
                    completeness: bars.completeness,
                    contamination: self.contamination_bar,
                    max_bin_size: self.max_bin_size,
                },
            );
            debug!("Join: {ledger}");
            self.census_bins(census, "join", &refiner.bins, &refiner.unbinned);
        }

        if let Some(margin) = self.recruit_near_bar {
            let _timer = crate::timing::scope("recruit");
            let ledger = crate::refine::recruit::recruit(
                &self.features(),
                &self.quality,
                &mut refiner.bins,
                crate::refine::recruit::RecruitSettings {
                    completeness: bars.completeness,
                    contamination: bars.contamination,
                    margin,
                    min_bin_size: self.min_bin_size,
                    worth: self.worth,
                },
            );
            debug!("Recruit: {ledger}");
            self.census_bins(census, "recruit", &refiner.bins, &refiner.unbinned);
        }

        {
            let report = crate::quality::write_report(
                &self.quality,
                refiner
                    .bins
                    .iter()
                    .map(|(bin, contigs)| (format!("rosella_bin_{bin}"), contigs.as_slice())),
                &self.coverage_table.contig_lengths,
                &path::Path::new(&self.output_directory).join(crate::defaults::QUALITY_FILE),
            );
            if let Err(error) = report {
                warn!("Could not write the quality table: {error}");
            }
        }

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
            completeness: self.quality.completeness_bar(self.min_completeness),
            contamination: self.contamination_bar,
            worth: self.worth,
            rung_floor: self.rung_floor,
        }
    }

    fn pick_partition(&self, ladder: Vec<Partitioning>, contigs: &[usize]) -> Partitioning {
        let judge = Judge {
            quality: &self.quality,
            contigs,
            bars: self.bars(),
        };
        combine(best_per_arm(ladder, &judge), &judge)
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
            self.partition_resolution,
            self.partition_theta,
            rank_rungs,
        )
    }

    fn embed(&self, contigs: &[usize]) -> (crate::embedding::Graph, KnnGraph) {
        let features = self.features();
        let knn = features.knn_of(
            contigs,
            self.n_neighbours,
            self.seeds,
            self.knn_candidates,
            KNN_ASSEMBLY,
        );
        let graph = features.graph_from_knn(contigs, &knn);
        (graph, knn)
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
        let knn = features.knn_of(&order, n_neighbours, self.seeds, self.knn_candidates, KNN_POOL);
        Ok((knn, order))
    }

    fn evaluate_subset(
        &self,
        knn: &KnnGraph,
        ordered_indices: &[usize],
        round: RoundParams,
    ) -> Result<Vec<Partitioning>> {
        let contig_id_map = ordered_indices
            .iter()
            .enumerate()
            .map(|(position, index)| (position, *index))
            .collect::<HashMap<_, _>>();

        let subset_graph = self
            .features()
            .graph_from_knn(ordered_indices, knn);
        let kind = match round.ladder && !self.partition.reads_ladder() {
            true => Partition::Leiden,
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
            result.reindex_clusters(contig_id_map.clone());
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
    }
}
