use std::{
    collections::{BTreeMap, HashMap, HashSet},
    path,
};

use anyhow::Result;
use log::{debug, info, warn};
use ndarray::Array2;

use crate::{
    cli::RecoverArgs,
    clustering::{
        clusterer::{HDBSCANResult, conserved, find_best_clusters, find_best_partition},
        contract::Contraction,
        graph_partition::{NodeSize, Partition},
        objective::{ClusterObjective, Objective, ObjectiveChoice},
    },
    coverage::coverage_table::CoverageTable,
    embedding::{
        features::ContigFeatures,
        metrics::{DistanceSettings, View},
        umap::EmbedOverrides,
    },
    homology::Homology,
    kmers::kmer_counting::KmerFrequencyTable,
    kmers::sketch::ContigSketches,
    markers::ContigMarkers,
    recover::census::{Census, STAGES_FILE},
    recover::inputs::{Inputs, read_inputs},
    recover::settings::{embed_overrides, seeds},
    refine::{
        bin_stats::LevelSource,
        dissolve::RoundParams,
        duplication::DuplicationSettings,
        gates::SplitGate,
        merger::{MergeBar, MergeSettings},
        solo::SoloPool,
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
    merge: bool,
    merge_singles: bool,
    merge_settings: MergeSettings,
    recruit: bool,
    eject: bool,
    eject_factor: f64,
    eject_duplicated: bool,
    duplication: DuplicationSettings,
    sketches: Option<ContigSketches>,
    pub(crate) overrides: EmbedOverrides,
    pub(crate) distance: DistanceSettings,
    pub(crate) largest_cluster: usize,
    gate: SplitGate,
    bisect: bool,
    solo: bool,
    solo_scatter: bool,
    solo_pool: SoloPool,
    homology_trigger: bool,
    fusion_bar: f64,
    dissolve: bool,
    dissolve_scope: crate::refine::dissolve::DissolveScope,
    dissolve_rounds: usize,
    dissolve_ladder: bool,
    dissolve_improve: bool,
    dissolve_select: crate::refine::select::Selection,
    min_completeness: f64,
    max_completeness_contamination: f64,
    quality: Option<crate::quality::ContigQuality>,
    recruit_rescued: bool,
    levels: LevelSource,
    level_quantile: f64,
    partition: Partition,
    node_size: NodeSize,
    partition_resolution: Option<f64>,
    partition_theta: Option<f64>,
    knn_report: Option<std::path::PathBuf>,
    homology: Option<Homology>,
    components: Option<Vec<usize>>,
    markers: Option<ContigMarkers>,
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
            homology,
            components,
            markers,
            quality,
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
            merge: !args.no_merge,
            merge_singles: !args.binning.no_merge_singles,
            merge_settings: MergeSettings {
                bar: MergeBar::parse(&args.binning.merge_bar).expect("clap restricts the value"),
                mutual: args.binning.merge_mutual,
                short_side: args.binning.merge_short_side,
                ..MergeSettings::default()
            },
            recruit: !args.no_recruit,
            eject: !args.no_eject,
            eject_factor: args.eject_factor,
            eject_duplicated: !args.no_eject_duplicated,
            duplication: DuplicationSettings {
                bar: args.duplication_bar,
                link: args.duplication_link,
                min_hashes: args.duplication_min_hashes,
            },
            sketches,
            overrides: embed_overrides(&args.overrides),
            distance,
            largest_cluster: args.binning.max_cluster_size,
            gate: SplitGate::parse(&args.binning.split_gate).expect("clap restricts the value"),
            bisect: args.binning.bisect,
            solo: !args.binning.no_solo,
            solo_scatter: !args.binning.no_solo_scatter,
            solo_pool: SoloPool::parse(&args.binning.solo_pool).expect("clap restricts the value"),
            homology_trigger: args.binning.homology_trigger,
            fusion_bar: args.fusion_bar,
            dissolve: !args.no_dissolve,
            dissolve_scope: crate::refine::dissolve::DissolveScope::parse(&args.dissolve_scope)
                .ok_or_else(|| anyhow!("unknown dissolve scope {}", args.dissolve_scope))?,
            dissolve_rounds: args.dissolve_rounds as usize,
            dissolve_ladder: args.dissolve_ladder,
            dissolve_improve: args.dissolve_improve,
            dissolve_select: crate::refine::select::Selection::parse(&args.dissolve_select)
                .expect("clap restricts the value"),
            min_completeness: args.min_completeness,
            max_completeness_contamination: args.max_contamination,
            quality,
            recruit_rescued: args.recruit_rescued,
            partition,
            node_size: NodeSize::parse(&args.binning.node_size).expect("clap restricts the value"),
            partition_resolution: args.binning.partition_resolution,
            partition_theta: args.binning.partition_theta,
            knn_report: args.binning.knn_report.clone(),
            homology,
            components,
            markers,
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
        let (embeddings, graph) = self.embed(&all_contigs)?;

        info!("Clustering.");
        let mut hdbscan_result = self.partition_of(&graph, embeddings.as_ref(), &all_contigs)?;
        debug!("HDBSCAN score {}", hdbscan_result.score);
        debug!(
            "HDBSCAN outlier percentage: {}",
            hdbscan_result.outliers.len() as f64 / self.n_contigs as f64
        );

        let mut census = Census::default();
        self.census_of(&mut census, "partition", &hdbscan_result);

        if self.recruit {
            let outliers = std::mem::take(&mut hdbscan_result.outliers)
                .into_iter()
                .collect::<Vec<_>>();
            let mut bins = hdbscan_result
                .cluster_map
                .iter()
                .map(|(id, contigs)| {
                    let mut contigs = contigs.iter().copied().collect::<Vec<_>>();
                    contigs.sort_unstable();
                    (*id, contigs)
                })
                .collect::<BTreeMap<_, _>>();
            let (left_over, recruited) = crate::refine::recruit::recruit(
                &self.features(),
                &mut bins,
                outliers,
                self.seeds.sample,
            );
            info!("Recruited {recruited} outliers into existing bins.");
            hdbscan_result.cluster_map = bins
                .into_iter()
                .map(|(id, contigs)| (id, contigs.into_iter().collect()))
                .collect();
            hdbscan_result.outliers = left_over.into_iter().collect();
            self.census_of(&mut census, "recruit", &hdbscan_result);
        }

        info!("Rescuing unbinned.");
        self.evaluate_outliers(&mut hdbscan_result)?;
        info!(
            "HDBSCAN outlier percentage: {}",
            hdbscan_result.outliers.len() as f64 / self.n_contigs as f64
        );
        self.census_of(&mut census, "outlier_pool", &hdbscan_result);

        if self.max_retries > 0 {
            info!("Refining bins.");
        }
        let (cluster_map, outliers) =
            self.refine_clusters(hdbscan_result, embeddings.as_ref(), &mut census);

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

    fn census_of(&self, census: &mut Census, stage: &'static str, result: &HDBSCANResult) {
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

    fn evaluate_outliers(&self, hdbscan_result: &mut HDBSCANResult) -> Result<()> {
        let outliers = std::mem::take(&mut hdbscan_result.outliers);
        if outliers.len() < MIN_RESCUE_CONTIGS {
            hdbscan_result.outliers = outliers;
            return Ok(());
        }
        let hdbscan_result_of_filtered_contigs = self.evaluate_subset(
            &outliers,
            RoundParams {
                n_neighbours: self.n_neighbours,
            },
        )?;

        debug!(
            "New HDBSCAN score {}",
            hdbscan_result_of_filtered_contigs.score
        );
        debug!(
            "Number of clusters: {}",
            hdbscan_result_of_filtered_contigs.cluster_map.len()
        );
        hdbscan_result.merge(hdbscan_result_of_filtered_contigs);

        Ok(())
    }

    fn refine_clusters(
        &self,
        hdbscan_result: HDBSCANResult,
        embeddings: Option<&Array2<f64>>,
        census: &mut Census,
    ) -> (HashMap<usize, HashSet<usize>>, HashSet<usize>) {
        let bins = hdbscan_result
            .cluster_map
            .into_iter()
            .map(|(bin_id, contigs)| {
                let mut contigs = contigs.into_iter().collect::<Vec<_>>();
                contigs.sort_unstable();
                (bin_id, contigs)
            })
            .collect::<BTreeMap<_, _>>();
        let mut unbinned = hdbscan_result.outliers.into_iter().collect::<Vec<_>>();
        unbinned.sort_unstable();

        let settings = RefineSettings {
            min_bin_size: self.min_bin_size,
            max_bin_size: self.max_bin_size,
            n_neighbours: self.n_neighbours,
            max_retries: self.max_retries,
            seeds: self.seeds,
            max_contamination: None,
            gate: self.gate,
            bisect: self.bisect,
            solo: self.solo,
            solo_scatter: self.solo_scatter,
            solo_pool: self.solo_pool,
            homology_trigger: self.homology_trigger,
            fusion_bar: self.fusion_bar,
            levels: self.levels,
            level_quantile: self.level_quantile,
            partition: self.partition,
            node_size: self.node_size,
            partition_resolution: self.partition_resolution,
            partition_theta: self.partition_theta,
            overrides: self.overrides,
            largest_cluster: self.largest_cluster,
        };
        let scorer = self.scorer();
        let mut refiner = Refiner::new(
            self.features(),
            embeddings,
            &scorer,
            settings,
            bins,
            unbinned,
        );
        refiner.run();
        self.census_bins(census, "refine", &refiner.bins, &refiner.unbinned);

        if self.merge {
            info!("Merging bins.");
            let settings = MergeSettings {
                genome_floor: self.merge_singles.then_some(refiner.genome_floor).flatten(),
                max_bin_size: self.max_bin_size,
                seed: self.seeds.sample,
                ..self.merge_settings
            };
            let (merged, merges) = crate::refine::merger::merge_bins(
                &self.features(),
                std::mem::take(&mut refiner.bins),
                settings,
            );
            info!("Merged {merges} pairs of bins.");
            refiner.bins = merged;
            self.census_bins(census, "merge", &refiner.bins, &refiner.unbinned);
        }

        if self.eject {
            let ejected = crate::refine::eject::eject(
                &self.features(),
                &mut refiner.bins,
                self.levels,
                self.level_quantile,
                self.eject_factor,
                self.min_bin_size,
                self.seeds.sample,
            );
            info!(
                "Ejected {} contigs sitting outside their bin.",
                ejected.len()
            );
            refiner.unbinned.extend(ejected);
            self.census_bins(census, "eject", &refiner.bins, &refiner.unbinned);
        }

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
                min_bin_size: self.min_bin_size,
                genome_floor: refiner.genome_floor,
                duplication_bar: self.duplication.bar,
                min_contigs: MIN_RESCUE_CONTIGS,
                scope: self.dissolve_scope,
                rounds: self.dissolve_rounds,
                ladder: self.dissolve_ladder,
                n_neighbours: self.n_neighbours,
                completeness: self.min_completeness,
                contamination: self.max_completeness_contamination,
                improve: self.dissolve_improve,
                select: self.dissolve_select,
            };
            let ledger = crate::refine::dissolve::dissolve(
                &self.features(),
                self.quality.as_ref(),
                &mut refiner.bins,
                &mut refiner.unbinned,
                settings,
                |pool, round| self.evaluate_subset(pool, round),
            );
            info!("Dissolve pool: {ledger}");
            self.census_bins(census, "dissolve", &refiner.bins, &refiner.unbinned);

            // Off by default: it wins on multi sample and costs far more single sample,
            // because adopting the pool's refusals back into a bin is how a bin turns impure.
            if self.recruit_rescued {
                let (left_over, recruited) = crate::refine::recruit::recruit(
                    &self.features(),
                    &mut refiner.bins,
                    std::mem::take(&mut refiner.unbinned),
                    self.seeds.sample,
                );
                info!("Recruited {recruited} of the pool's leftovers into surviving bins.");
                refiner.unbinned = left_over;
                self.census_bins(census, "recruit_rescued", &refiner.bins, &refiner.unbinned);
            }
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

    /// Embed and cluster a subset of contigs. `contig_indices` are indices into the contig
    /// list as it stands after the initial length filter.
    fn partition_of(
        &self,
        graph: &crate::embedding::Graph,
        embeddings: Option<&ndarray::Array2<f64>>,
        contigs: &[usize],
    ) -> Result<HDBSCANResult> {
        if self.partition.reads_graph() {
            let lengths = self.features().contig_lengths(contigs);
            let run = |graph: &crate::embedding::Graph, nodes: &[usize], lengths: &[usize]| {
                find_best_partition(
                    graph,
                    embeddings,
                    nodes,
                    lengths,
                    self.node_size,
                    &self.scorer(),
                    self.seeds.sample,
                    self.seeds.partition,
                    self.partition,
                    self.partition_resolution,
                    self.partition_theta,
                )
            };
            match self
                .components
                .as_ref()
                .filter(|_| !self.scorer().needs_layout())
                .and_then(|component| Contraction::new(component, contigs))
            {
                None => run(graph, contigs, &lengths),
                Some(contraction) => {
                    let held = contraction.graph(graph);
                    let nodes = (0..contraction.len()).collect::<Vec<_>>();
                    run(&held, &nodes, &contraction.lengths(&lengths))
                        .map(|result| contraction.expand(result))
                }
            }
        } else {
            let embeddings = embeddings
                .ok_or_else(|| anyhow!("HDBSCAN clusters a layout but none was built"))?;
            find_best_clusters(
                embeddings,
                contigs,
                &self.scorer(),
                self.seeds.sample,
                self.largest_cluster,
            )
        }
    }

    /// The layout is only ever read by a score that needs one, so the graph arms under a
    /// graph score take the manifold alone and never pay for the SGD.
    fn embed(
        &self,
        contigs: &[usize],
    ) -> Result<(Option<ndarray::Array2<f64>>, crate::embedding::Graph)> {
        self.embed_with(contigs, self.n_neighbours)
    }

    fn embed_with(
        &self,
        contigs: &[usize],
        n_neighbours: usize,
    ) -> Result<(Option<ndarray::Array2<f64>>, crate::embedding::Graph)> {
        let features = self.features();
        if self.wants_layout() {
            let (embeddings, graph) =
                features.embed_with_graph(contigs, n_neighbours, self.seeds, &self.overrides)?;
            Ok((Some(embeddings), graph))
        } else {
            let graph = features.graph_of(contigs, n_neighbours, self.seeds, &self.overrides);
            Ok((None, graph))
        }
    }

    fn write_knn_report(&self, contigs: &[usize], path: &path::Path) -> Result<()> {
        let knn = self
            .features()
            .knn_of(contigs, self.n_neighbours, self.seeds, &self.overrides);
        let views = self.distance.views.selected();
        let labels = if views.is_empty() {
            vec![crate::embedding::metrics::VIEW_NAMES[0]]
        } else {
            views.iter().map(View::name).collect::<Vec<_>>()
        };
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

    fn wants_layout(&self) -> bool {
        !self.partition.reads_graph() || self.scorer().needs_layout()
    }

    fn evaluate_subset(
        &self,
        contig_indices: &HashSet<usize>,
        round: RoundParams,
    ) -> Result<HDBSCANResult> {
        let mut ordered_indices = contig_indices.iter().copied().collect::<Vec<_>>();
        ordered_indices.sort_unstable();
        let contig_id_map = ordered_indices
            .iter()
            .enumerate()
            .map(|(position, index)| (position, *index))
            .collect::<HashMap<_, _>>();

        let (subset_embeddings, subset_graph) =
            self.embed_with(&ordered_indices, round.n_neighbours)?;
        let mut hdbscan_result =
            self.partition_of(&subset_graph, subset_embeddings.as_ref(), &ordered_indices)?;
        debug!("HDBSCAN score {}", hdbscan_result.score);

        hdbscan_result.reindex_clusters(contig_id_map);

        conserved(
            hdbscan_result
                .cluster_map
                .values()
                .flatten()
                .copied()
                .chain(hdbscan_result.outliers.iter().copied()),
            contig_indices,
        )?;

        Ok(hdbscan_result)
    }

    fn scorer(&self) -> Objective<'_> {
        self.objective
            .build(&self.coverage_table.contig_lengths, self.min_bin_size)
    }

    fn features(&self) -> ContigFeatures<'_> {
        ContigFeatures::new(
            &self.coverage_table.table,
            &self.tnf_table.kmer_table,
            &self.coverage_table.contig_lengths,
        )
        .with_distance(self.distance)
        .with_homology(self.homology.as_ref())
        .with_sketches(self.sketches.as_ref())
        .with_markers(self.markers.as_ref())
        .with_bands(self.n_neighbours)
    }
}
