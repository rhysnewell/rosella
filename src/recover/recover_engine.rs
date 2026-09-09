use std::{
    collections::{BTreeMap, HashMap, HashSet},
    path,
};

use anyhow::Result;
use log::{debug, info, warn};

use crate::{
    cli::RecoverArgs,
    clustering::{
        clusterer::{Partitioning, conserved, find_partitions},
        graph_partition::{NodeSize, Partition},
        objective::{Objective, ObjectiveChoice},
    },
    coverage::coverage_table::CoverageTable,
    embedding::{
        features::ContigFeatures,
        metrics::DistanceSettings,
        umap::EmbedOverrides,
    },
    kmers::kmer_counting::KmerFrequencyTable,
    kmers::sketch::ContigSketches,
    recover::census::{Census, STAGES_FILE},
    recover::inputs::{Inputs, read_inputs},
    recover::settings::{embed_overrides, seeds},
    refine::{
        bin_stats::LevelSource,
        dissolve::RoundParams,
        duplication::DuplicationSettings,
        splitter::{RefineSettings, Refiner},
    },
    seeds::Seeds,
};

pub const RECOVER_FASTA_EXTENSION: &str = ".fna";
pub const UNBINNED: &str = "unbinned";
pub(crate) const REFINING_BIN_SIZE: usize = 1000000;

/// umap-rs asks for two neighbours, and a subset of three is the smallest that has them.
const QUALITY_FILE: &str = "quality.tsv";
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
    objective: ObjectiveChoice,
    pub(crate) min_contig_size: usize,
    pub(crate) max_bin_size: usize,
    pub(crate) max_retries: usize,
    eject_duplicated: bool,
    duplication: DuplicationSettings,
    sketches: Option<ContigSketches>,
    pub(crate) overrides: EmbedOverrides,
    pub(crate) distance: DistanceSettings,
    bisect: bool,
    dissolve: bool,
    dissolve_rounds: usize,
    dissolve_passes: usize,
    min_completeness: f64,
    max_completeness_contamination: f64,
    quality: Option<crate::quality::ContigQuality>,
    oracle: Vec<Vec<usize>>,
    levels: LevelSource,
    level_quantile: f64,
    partition: Partition,
    node_size: NodeSize,
    partition_resolution: Option<f64>,
    partition_theta: Option<f64>,
    knn_report: Option<std::path::PathBuf>,
}

impl RecoverEngine {
    pub fn new(args: &RecoverArgs) -> Result<Self> {
        let Inputs {
            output_directory,
            assembly,
            min_contig_size,
            coverage_table,
            tnf_table,
            sketches,
            quality,
            oracle,
            distance,
            partition,
        } = read_inputs(args)?;

        let n_neighbours = args.binning.n_neighbours;
        let seeds = seeds(args.common.seed, &args.seeds);
        let min_bin_size = args.binning.min_bin_size;

        let n_contigs = coverage_table.table.nrows();
        let max_bin_size = args.binning.max_bin_size;
        let max_retries = if args.no_refine {
            0
        } else {
            args.binning.max_retries
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
            eject_duplicated: !args.no_eject_duplicated,
            duplication: DuplicationSettings {
                bar: args.duplication_bar,
                link: args.duplication_link,
                min_hashes: args.duplication_min_hashes,
            },
            sketches,
            overrides: embed_overrides(&args.overrides),
            distance,
            bisect: args.binning.bisect,
            dissolve: !args.no_dissolve,
            dissolve_rounds: args.dissolve_rounds as usize,
            dissolve_passes: args.dissolve_passes as usize,
            min_completeness: args.min_completeness,
            max_completeness_contamination: args.max_contamination,
            quality,
            oracle,
            partition,
            node_size: NodeSize::parse(&args.binning.node_size).expect("clap restricts the value"),
            partition_resolution: args.binning.partition_resolution,
            partition_theta: args.binning.partition_theta,
            knn_report: args.binning.knn_report.clone(),
            levels: LevelSource::parse(&args.binning.split_levels)
                .expect("clap restricts the value"),
            level_quantile: args.binning.split_level_quantile,
            objective: ObjectiveChoice::parse(&args.binning.objective)
                .ok_or_else(|| anyhow::anyhow!("unknown objective {}", args.binning.objective))?,
        })
    }

    /// Runs through the rosella bin recovery pipeline
    pub fn run(self) -> Result<()> {
        let all_contigs = (0..self.n_contigs).collect::<Vec<usize>>();

        if let Some(path) = &self.knn_report {
            self.write_knn_report(&all_contigs, path)?;
            info!("Wrote the kNN report to {}.", path.display());
            return Ok(());
        }

        info!("Embedding.");
        let graph = self.embed(&all_contigs);

        info!("Clustering.");
        let mut partitioning = self
            .partition_of(&graph, &all_contigs, self.partition)?
            .swap_remove(0);
        debug!("Partition score {}", partitioning.score);
        debug!(
            "Outlier percentage: {}",
            partitioning.outliers.len() as f64 / self.n_contigs as f64
        );

        let mut census = Census::default();
        self.census_of(&mut census, "partition", &partitioning);

        info!("Rescuing unbinned.");
        self.evaluate_outliers(&mut partitioning)?;
        info!(
            "Outlier percentage: {}",
            partitioning.outliers.len() as f64 / self.n_contigs as f64
        );
        self.census_of(&mut census, "outlier_pool", &partitioning);

        if self.max_retries > 0 {
            info!("Refining bins.");
        }
        let (cluster_map, outliers) =
            self.refine_clusters(partitioning, &mut census);

        conserved(
            cluster_map
                .values()
                .flatten()
                .copied()
                .chain(outliers.iter().copied()),
            &all_contigs.iter().copied().collect(),
        )?;
        let cluster_results = self.get_cluster_result(cluster_map, outliers, None);
        info!("Length of cluster results: {}", cluster_results.len());

        info!("Writing clusters.");
        {
            let _timer = crate::timing::scope("write");
            self.write_clusters(cluster_results, false)?;
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

    fn evaluate_outliers(&self, partitioning: &mut Partitioning) -> Result<()> {
        let outliers = std::mem::take(&mut partitioning.outliers);
        if outliers.len() < MIN_RESCUE_CONTIGS {
            partitioning.outliers = outliers;
            return Ok(());
        }
        let partitioning_of_filtered_contigs = self
            .evaluate_subset(
                &outliers,
                RoundParams {
                    n_neighbours: self.n_neighbours,
                    ladder: false,
                },
            )?
            .swap_remove(0);

        debug!(
            "New Partition score {}",
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
            levels: self.levels,
            level_quantile: self.level_quantile,
            partition: self.partition,
            node_size: self.node_size,
            partition_resolution: self.partition_resolution,
            partition_theta: self.partition_theta,
            overrides: self.overrides,
        };
        let scorer = self.scorer();
        let mut refiner = Refiner::new(
            self.features(),
            &scorer,
            settings,
            bins,
            unbinned,
        );
        refiner.run();
        self.census_bins(census, "refine", &refiner.bins, &refiner.unbinned);

        if self.eject_duplicated {
            let ejected = crate::refine::duplication::eject_duplicated(
                &self.features(),
                &mut refiner.bins,
                self.duplication,
                self.min_bin_size,
            );
            info!("Ejected {} contigs their bin holds twice.", ejected.len());
            refiner.unbinned.extend(ejected);
            self.census_bins(census, "eject_duplicated", &refiner.bins, &refiner.unbinned);
        }

        if self.dissolve {
            // Stale by a round, since merge and both eject arms move the bins it was
            // measured on. Recomputing it here was measured and lost bins.
            let settings = crate::refine::dissolve::DissolveSettings {
                bars: crate::refine::rung::Bars {
                    min_bin_size: self.min_bin_size,
                    duplication_bar: self.duplication.bar,
                    completeness: self.min_completeness,
                    contamination: self.max_completeness_contamination,
                },
                genome_floor: refiner.genome_floor,
                min_contigs: MIN_RESCUE_CONTIGS,
                rounds: self.dissolve_rounds,
                passes: self.dissolve_passes,
                n_neighbours: self.n_neighbours,
            };
            let ledger = crate::refine::dissolve::dissolve(
                &self.features(),
                self.quality.as_ref(),
                &mut refiner.bins,
                &mut refiner.unbinned,
                settings,
                &self.oracle,
                |pool, round| self.evaluate_subset(pool, round),
            );
            info!("Dissolve pool: {ledger}");
            self.census_bins(census, "dissolve", &refiner.bins, &refiner.unbinned);

        }

        if let Some(quality) = self.quality.as_ref() {
            let ledger = crate::refine::join::join(
                &self.features(),
                quality,
                &mut refiner.bins,
                crate::refine::join::JoinSettings {
                    completeness: self.min_completeness,
                    contamination: self.max_completeness_contamination,
                    max_bin_size: self.max_bin_size,
                },
            );
            info!("Join: {ledger}");
            self.census_bins(census, "join", &refiner.bins, &refiner.unbinned);
        }

        if let Some(quality) = self.quality.as_ref() {
            let report = quality.write_report(
                &refiner.bins,
                &self.coverage_table.contig_lengths,
                &path::Path::new(&self.output_directory).join(QUALITY_FILE),
            );
            if let Err(error) = report {
                warn!("Could not write {QUALITY_FILE}: {error}");
            }
        }

        let cluster_map = refiner
            .bins
            .iter()
            .map(|(bin_id, contigs)| (*bin_id, contigs.iter().copied().collect::<HashSet<_>>()))
            .collect::<HashMap<_, _>>();
        (cluster_map, refiner.unbinned.iter().copied().collect())
    }

    /// Partition a subset of contigs. `contigs` are indices into the contig list as it
    /// stands after the initial length filter.
    fn partition_of(
        &self,
        graph: &crate::embedding::Graph,
        contigs: &[usize],
        kind: Partition,
    ) -> Result<Vec<Partitioning>> {
        find_partitions(
            graph,
            &self.features().contig_lengths(contigs),
            self.node_size,
            &self.scorer(),
            self.seeds.partition,
            kind,
            self.partition_resolution,
            self.partition_theta,
        )
    }

    fn embed(&self, contigs: &[usize]) -> crate::embedding::Graph {
        self.embed_with(contigs, self.n_neighbours)
    }

    fn embed_with(&self, contigs: &[usize], n_neighbours: usize) -> crate::embedding::Graph {
        self.features()
            .graph_of(contigs, n_neighbours, self.seeds, &self.overrides)
    }

    fn write_knn_report(&self, contigs: &[usize], path: &path::Path) -> Result<()> {
        let knn = self
            .features()
            .knn_of(contigs, self.n_neighbours, self.seeds, &self.overrides);
        let labels = vec!["combined"];
        let names = contigs
            .iter()
            .map(|index| self.coverage_table.contig_names[*index].as_str())
            .collect::<Vec<_>>();
        crate::embedding::knn::write_report(
            &labels.into_iter().zip(knn).collect::<Vec<_>>(),
            &names,
            path,
        )
    }

    fn evaluate_subset(
        &self,
        contig_indices: &HashSet<usize>,
        round: RoundParams,
    ) -> Result<Vec<Partitioning>> {
        let mut ordered_indices = contig_indices.iter().copied().collect::<Vec<_>>();
        ordered_indices.sort_unstable();
        let contig_id_map = ordered_indices
            .iter()
            .enumerate()
            .map(|(position, index)| (position, *index))
            .collect::<HashMap<_, _>>();

        let subset_graph = self.embed_with(&ordered_indices, round.n_neighbours);
        // The size rule that sends a large assembly to label propagation is about the assembly,
        // and the pool is a fraction of it, so a ladder is available here either way.
        let kind = match round.ladder && !self.partition.reads_ladder() {
            true => Partition::Leiden,
            false => self.partition,
        };
        let mut results = self.partition_of(&subset_graph, &ordered_indices, kind)?;
        if !round.ladder {
            results.truncate(1);
        }
        debug!("Partition score {}", results[0].score);

        for result in results.iter_mut() {
            result.reindex_clusters(contig_id_map.clone());
            conserved(
                result
                    .cluster_map
                    .values()
                    .flatten()
                    .copied()
                    .chain(result.outliers.iter().copied()),
                contig_indices,
            )?;
        }

        Ok(results)
    }

    fn scorer(&self) -> Objective {
        self.objective.build()
    }

    fn features(&self) -> ContigFeatures<'_> {
        ContigFeatures::new(
            &self.coverage_table.table,
            &self.tnf_table.kmer_table,
            &self.coverage_table.contig_lengths,
        )
        .with_distance(self.distance)
        .with_sketches(self.sketches.as_ref())
    }
}
