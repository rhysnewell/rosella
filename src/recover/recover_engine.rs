use std::{
    collections::{BTreeMap, HashMap, HashSet},
    path,
};

use anyhow::Result;
use log::{debug, info};
use ndarray::Array2;

use crate::{
    cli::RecoverArgs,
    clustering::{
        clusterer::{HDBSCANResult, find_best_clusters},
        objective::{Dbcv, ObjectiveChoice},
    },
    coverage::{
        coverage_calculator::{CoverageInputs, calculate_coverage},
        coverage_table::CoverageTable,
    },
    embedding::{
        features::ContigFeatures,
        metrics::{CoverageAggregation, DistanceSettings},
        umap::EmbedOverrides,
    },
    kmers::kmer_counting::{KmerFrequencyTable, count_kmers},
    refine::splitter::{RefineSettings, Refiner},
    seeds::Seeds,
};

pub const RECOVER_FASTA_EXTENSION: &str = ".fna";
pub const UNBINNED: &str = "unbinned";
pub(crate) const REFINING_BIN_SIZE: usize = 1000000;

pub fn embed_overrides(overrides: &crate::cli::EmbeddingOverrides) -> EmbedOverrides {
    EmbedOverrides {
        a: overrides.umap_a,
        b: overrides.umap_b,
        n_components: overrides.n_components,
        n_epochs: overrides.n_epochs,
        length_weight: overrides.length_weight,
    }
}

pub fn seeds(seed: u64, overrides: &crate::cli::SeedOverrides) -> Seeds {
    Seeds {
        knn: overrides.knn.unwrap_or(seed),
        init: overrides.init.unwrap_or(seed),
        layout: overrides.layout.unwrap_or(seed),
        sample: overrides.sample.unwrap_or(seed),
    }
}

pub fn distance_settings(distance: &crate::cli::DistanceParams) -> DistanceSettings {
    DistanceSettings {
        aggregation: CoverageAggregation::parse(&distance.coverage_aggregation)
            .expect("clap restricts the value"),
        length_scaled_variance: distance.length_scaled_variance,
    }
}

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
    pub(crate) filtered_contigs: HashSet<String>,
    pub(crate) max_bin_size: usize,
    pub(crate) max_retries: usize,
    pub(crate) overrides: EmbedOverrides,
    pub(crate) distance: DistanceSettings,
    pub(crate) largest_cluster: usize,
}

impl RecoverEngine {
    pub fn new(args: &RecoverArgs) -> Result<Self> {
        let output_directory = args.common.output_directory.clone();
        // check if output_directory contains .fna files, if so exit
        let output_directory_path = path::Path::new(&output_directory);
        if output_directory_path.exists() {
            let output_directory_files = output_directory_path.read_dir()?;
            for file in output_directory_files {
                let file = file?;
                let file_name = file.file_name();
                let file_name = file_name.to_str().unwrap();
                if file_name.ends_with(RECOVER_FASTA_EXTENSION) {
                    return Err(anyhow::anyhow!(
                        "Output directory contains .fna files. Please remove them before running rosella recover."
                    ));
                }
            }
        }

        let assembly = args.assembly.clone();
        std::fs::create_dir_all(&output_directory)?;
        info!("Calculating contig coverages.");
        let min_contig_size = args.binning.min_contig_size;
        let mut coverage_table = {
            let _timer = crate::timing::scope("coverage");
            calculate_coverage(&CoverageInputs {
                assembly: Some(&assembly),
                output_directory: &output_directory,
                threads: args.common.threads,
                coverage: &args.coverage,
                mapping: &args.mapping,
                filtering: &args.filtering,
                alignment: &args.alignment,
                trimming: &args.trimming,
            })?
        };
        let n_contigs = coverage_table.table.nrows();

        let filtered_contigs = {
            let _timer = crate::timing::scope("length_filter");
            coverage_table.filter_by_length(min_contig_size)?
        };

        assert_eq!(
            coverage_table.table.nrows(),
            n_contigs - filtered_contigs.len(),
            "Coverage table row count and total contigs minus filtered contigs do not match."
        );
        let mut tnf_table = {
            let _timer = crate::timing::scope("kmers");
            if let Some(kmer_table_path) = &args.common.kmer_frequency_file {
                info!("Reading TNF table.");
                KmerFrequencyTable::read(kmer_table_path)?
            } else {
                info!("Calculating TNF table.");
                count_kmers(&assembly, &output_directory, Some(n_contigs))?
            }
        };
        assert_eq!(
            n_contigs,
            tnf_table.kmer_table.nrows(),
            "Coverage table row count and TNF table row count do not match."
        );
        debug!("Filtering TNF table.");
        tnf_table.filter_by_name(&filtered_contigs)?;
        assert_eq!(
            coverage_table.table.nrows(),
            tnf_table.kmer_table.nrows(),
            "Coverage table and TNF table have different number of contigs."
        );

        info!(
            "{} valid contigs, {} filtered contigs.",
            coverage_table.table.nrows(),
            filtered_contigs.len()
        );
        let n_neighbours = args.binning.n_neighbours;
        let seeds = seeds(args.common.seed, &args.seeds);
        let min_bin_size = args.binning.min_bin_size;

        let n_contigs = coverage_table.table.nrows();
        let max_bin_size = args.binning.max_bin_size;
        let max_retries = if args.refine {
            args.binning.max_retries
        } else {
            0
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
            // filtered_contigs,
            filtered_contigs: HashSet::new(),
            max_bin_size,
            max_retries,
            overrides: embed_overrides(&args.overrides),
            distance: distance_settings(&args.distance),
            largest_cluster: args.binning.max_cluster_size,
            objective: ObjectiveChoice::parse(&args.binning.objective)
                .ok_or_else(|| anyhow::anyhow!("unknown objective {}", args.binning.objective))?,
        })
    }

    /// Runs through the rosella bin recovery pipeline
    pub fn run(self) -> Result<()> {
        info!("Embedding.");
        let all_contigs = (0..self.n_contigs).collect::<Vec<usize>>();
        let embeddings =
            self.features()
                .embed(&all_contigs, self.n_neighbours, self.seeds, &self.overrides)?;

        info!("Clustering.");
        let mut hdbscan_result = find_best_clusters(
            &embeddings,
            &all_contigs,
            &self.scorer(),
            self.seeds.sample,
            self.largest_cluster,
        )?;
        debug!("HDBSCAN score {}", hdbscan_result.score);
        debug!(
            "HDBSCAN outlier percentage: {}",
            hdbscan_result.outliers.len() as f64 / self.n_contigs as f64
        );

        info!("Rescuing unbinned.");
        self.evaluate_outliers(&mut hdbscan_result)?;
        info!(
            "HDBSCAN outlier percentage: {}",
            hdbscan_result.outliers.len() as f64 / self.n_contigs as f64
        );

        if self.max_retries > 0 {
            info!("Refining bins.");
        }
        let (cluster_map, outliers) = self.refine_clusters(hdbscan_result, &embeddings);

        let n_contigs = cluster_map.values().map(|v| v.len()).sum::<usize>() + outliers.len();
        let cluster_results = self.get_cluster_result(cluster_map, outliers, None);
        info!("Length of cluster results: {}", cluster_results.len());
        debug!(
            "cluster result len {} n_contigs {} contigs used {} coverage table and kmer table size {} {} n filtered contigs {}",
            cluster_results.len(),
            n_contigs,
            self.n_contigs,
            self.coverage_table.contig_lengths.len(),
            self.tnf_table.contig_names.len(),
            self.filtered_contigs.len()
        );
        if n_contigs != cluster_results.len() {
            bail!(
                "Number of contigs in cluster results ({}) does not match number of contigs in HDBSCAN result ({})",
                cluster_results.len(),
                n_contigs
            );
        }

        info!("Writing clusters.");
        {
            let _timer = crate::timing::scope("write");
            self.write_clusters(cluster_results, false)?;
        }

        crate::timing::report(
            path::Path::new(&self.output_directory).join(crate::timing::TIMINGS_FILE),
        )?;

        Ok(())
    }

    fn evaluate_outliers(&self, hdbscan_result: &mut HDBSCANResult) -> Result<()> {
        let outliers = std::mem::take(&mut hdbscan_result.outliers);
        let hdbscan_result_of_filtered_contigs = self.evaluate_subset(&outliers)?;

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

    /// Split the chimeric bins, then hand back the cluster map the writer expects.
    fn refine_clusters(
        &self,
        hdbscan_result: HDBSCANResult,
        embeddings: &Array2<f64>,
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
            overrides: self.overrides,
            largest_cluster: self.largest_cluster,
        };
        let scorer = self.scorer();
        let mut refiner = Refiner::new(
            self.features(),
            Some(embeddings),
            &scorer,
            settings,
            bins,
            unbinned,
        );
        refiner.run();

        let cluster_map = refiner
            .bins
            .iter()
            .map(|(bin_id, contigs)| (*bin_id, contigs.iter().copied().collect::<HashSet<_>>()))
            .collect::<HashMap<_, _>>();
        (cluster_map, refiner.unbinned.iter().copied().collect())
    }

    /// Embed and cluster a subset of contigs. `contig_indices` are indices into the contig
    /// list as it stands after the initial length filter.
    fn evaluate_subset(&self, contig_indices: &HashSet<usize>) -> Result<HDBSCANResult> {
        let mut ordered_indices = contig_indices.iter().copied().collect::<Vec<_>>();
        ordered_indices.sort_unstable();
        let contig_id_map = ordered_indices
            .iter()
            .enumerate()
            .map(|(position, index)| (position, *index))
            .collect::<HashMap<_, _>>();

        let subset_embeddings = self.features().embed(
            &ordered_indices,
            self.n_neighbours,
            self.seeds,
            &self.overrides,
        )?;
        let mut hdbscan_result = find_best_clusters(
            &subset_embeddings,
            &ordered_indices,
            &self.scorer(),
            self.seeds.sample,
            self.largest_cluster,
        )?;
        debug!("HDBSCAN score {}", hdbscan_result.score);

        hdbscan_result.reindex_clusters(contig_id_map);

        let n_clustered_contigs = hdbscan_result
            .cluster_map
            .values()
            .map(|contigs| contigs.len())
            .sum::<usize>()
            + hdbscan_result.outliers.len();
        if n_clustered_contigs != contig_indices.len() {
            return Err(anyhow!(
                "Number of clustered contigs does not match number of contigs in subset. {} != {}",
                n_clustered_contigs,
                contig_indices.len()
            ));
        }

        Ok(hdbscan_result)
    }

    fn scorer(&self) -> Dbcv<'_> {
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
    }
}
