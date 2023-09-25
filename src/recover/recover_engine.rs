use std::{collections::{HashMap, HashSet}, cmp::{Ordering, max, min}, path, io::{BufWriter, Write}, fs::OpenOptions};


use anyhow::Result;
use log::{info, debug};
use ndarray::{Dim, ArrayBase, OwnedRepr, Array2};
use needletail::{parse_fastx_file, parser::{write_fasta, LineEnding}};
use rayon::{prelude::*, slice::ParallelSliceMut, prelude::IntoParallelIterator};

#[cfg(not(feature = "no_flight"))]
use crate::{
    coverage::{coverage_calculator::calculate_coverage, coverage_table::CoverageTable},
    kmers::kmer_counting::{KmerFrequencyTable, count_kmers}, external_command_checker::check_for_flight
};
#[cfg(not(feature = "no_flight"))]
use std::{process::Command, io::{BufRead, Read}};

#[cfg(feature = "no_flight")]
use annembed::{prelude::*, fromhnsw::kgraph::KGraph};
#[cfg(feature = "no_flight")]
use crate::{
    coverage::{coverage_calculator::{calculate_coverage, MetabatDistance}, coverage_table::CoverageTable}, 
    sketch::{contig_sketcher::{ContigSketchResult, sketch_contigs}, sketch_distances::distance}, 
    embedding::embedder::{ContigInformation, EmbedderEngine}, 
    clustering::clusterer::{find_best_clusters, HDBSCANResult}, 
    kmers::kmer_counting::{KmerFrequencyTable, count_kmers, KmerCorrelation},
};

const RECOVER_FASTA_EXTENSION: &str = ".fna";
const DEBUG_BINS: bool = true;
const UNBINNED: &str = "unbinned";
const DEFAULT_B: f64 = 1.5;
const MAX_ITERATIONS: usize = 5;

pub fn run_recover(m: &clap::ArgMatches) -> Result<()> {
    let mut recover_engine = RecoverEngine::new(m)?;
    recover_engine.run()?;
    Ok(())
}

struct RecoverEngine {
    output_directory: String,
    assembly: String,
    coverage_table: CoverageTable,
    #[cfg(feature = "no_flight")]
    contig_sketches: ContigSketchResult,
    tnf_table: KmerFrequencyTable,
    n_neighbours: usize,
    filtering_rounds: usize,
    ef_construction: usize,
    max_layers: usize,
    n_contigs: usize,
    min_bin_size: usize,
    min_contig_size: usize,
    filtered_contigs: HashSet<String>,
}

impl RecoverEngine {
    pub fn new(m: &clap::ArgMatches) -> Result<Self> {
        // create output directory
        let output_directory = m.get_one::<String>("output-directory").unwrap().clone();
        // check if output_directory contains .fna files, if so exit
        let output_directory_path = path::Path::new(&output_directory);
        if output_directory_path.exists() {
            let output_directory_files = output_directory_path.read_dir()?;
            for file in output_directory_files {
                let file = file?;
                let file_name = file.file_name();
                let file_name = file_name.to_str().unwrap();
                if file_name.ends_with(RECOVER_FASTA_EXTENSION) {
                    return Err(anyhow::anyhow!("Output directory contains .fna files. Please remove them before running rosella recover."));
                }
            }
        }

        let assembly = m.get_one::<String>("assembly").unwrap().clone();
        // create the output directory but do not fail if it already exists
        std::fs::create_dir_all(&output_directory)?;
        info!("Calculating contig coverages.");
        let min_contig_size = m.get_one::<usize>("min-contig-size").unwrap().clone();
        let mut coverage_table = calculate_coverage(m)?;
        let n_contigs = coverage_table.table.nrows();

        #[cfg(feature = "no_flight")]
        let filtered_contigs = coverage_table.filter_by_length(min_contig_size)?;

        #[cfg(not(feature = "no_flight"))]
        let filtered_contigs = HashSet::new();

        assert_eq!(coverage_table.table.nrows(), n_contigs - filtered_contigs.len(), "Coverage table row count and total contigs minus filtered contigs do not match.");
        info!("Calculating TNF table.");
        let mut tnf_table = count_kmers(m, Some(n_contigs))?;
        assert_eq!(n_contigs, tnf_table.kmer_table.nrows(), "Coverage table row count and TNF table row count do not match.");
        debug!("Filtering TNF table.");
        tnf_table.filter_by_name(&filtered_contigs)?;
        assert_eq!(coverage_table.table.nrows(), tnf_table.kmer_table.nrows(), "Coverage table and TNF table have different number of contigs.");
        
        #[cfg(feature = "no_flight")]
        info!("Calculating contig sketches.");
        #[cfg(feature = "no_flight")]
        let contig_sketches = sketch_contigs(m)?;
        #[cfg(feature = "no_flight")]
        assert_eq!(coverage_table.table.nrows(), contig_sketches.contig_sketches.len(), "Coverage table and contig sketches have different number of contigs.");
        
        info!("{} valid contigs, {} filtered contigs.", coverage_table.table.nrows(), filtered_contigs.len());
        let n_neighbours = m.get_one::<usize>("n-neighbours").unwrap().clone();
        let filtering_rounds = m.get_one::<usize>("filtering-rounds").unwrap().clone();
        let ef_construction = m.get_one::<usize>("ef-construction").unwrap().clone();
        let max_layers = m.get_one::<usize>("max-layers").unwrap().clone();
        let min_bin_size = m.get_one::<usize>("min-bin-size").unwrap().clone();
        
        let n_contigs = coverage_table.table.nrows();
        Ok(
            Self {
                output_directory,
                assembly,
                coverage_table,
                #[cfg(feature = "no_flight")]
                contig_sketches,
                tnf_table,
                n_neighbours,
                filtering_rounds,
                ef_construction,
                max_layers,
                n_contigs,
                min_bin_size,
                min_contig_size,
                // filtered_contigs,
                filtered_contigs: HashSet::new(),
            }
        )
    }

    /// Runs through the rosella bin recovery pipeline
    pub fn run(&mut self) -> Result<()> {
        #[cfg(feature = "no_flight")]
        {
            // 1. Perform N rounds of filtering. Filtering involves removinf disconnnected contigs from the graph.
            //    The final round ensures each contig has at least one neighbour.
            info!("Filtering.");
            let kgraph = self.prefilter_nodes(self.filtering_rounds)?;
            
            // 2. Take the KGraph and embed using a UMAP-esque method.
            //    Embedding parameters need to be chosen more carefully.
            info!("Embedding.");
            let (embeddings, new_kgraph) = self.embed(&kgraph, DEFAULT_B, 0)?;
            let _kgraph = if let Some(new_kgraph) = new_kgraph {
                new_kgraph
            } else {
                kgraph
            };
            debug!("Embeddings: {:?}", embeddings);
            // 3. Cluster the embeddings using HDBSCAN.
            //    As with the embedding parameters, we need to better choose the clustering parameters.
            info!("Clustering.");
            let mut hdbscan_result = find_best_clusters(&embeddings, 2, 5)?;
            debug!("HDBSCAN score {}", hdbscan_result.score);
            debug!("HDBSCAN outlier percentage: {}", hdbscan_result.outliers.len() as f64 / self.n_contigs as f64);

            info!("Rescuing unbinned.");
            self.evaluate_outliers(&mut hdbscan_result)?;
            info!("HDBSCAN outlier percentage: {}", hdbscan_result.outliers.len() as f64 / self.n_contigs as f64);
            
            let n_contigs = hdbscan_result.cluster_map.values().map(|v| v.len()).sum::<usize>() + hdbscan_result.outliers.len();
            let cluster_results = self.get_cluster_result(hdbscan_result.cluster_map, hdbscan_result.outliers, Some(&embeddings), None);
            info!("Length of cluster results: {}", cluster_results.len());
            debug!("cluster result len {} n_contigs {} contigs used {} coverage table and kmer table size {} {} n filtered contigs {}", 
                cluster_results.len(), n_contigs, self.n_contigs, self.coverage_table.contig_lengths.len(), self.tnf_table.contig_names.len(), self.filtered_contigs.len());
            if n_contigs != cluster_results.len() {
                bail!("Number of contigs in cluster results ({}) does not match number of contigs in HDBSCAN result ({})", cluster_results.len(), n_contigs);
            }

            info!("Writing clusters.");
            self.write_clusters(cluster_results)?;
        }

        #[cfg(not(feature = "no_flight"))]
        {
            self.run_flight()?;
            let cluster_results = self.finalise_bins()?;

            info!("Writing clusters.");
            self.write_clusters(cluster_results)?;
        }

        Ok(())
    }

    #[cfg(feature = "no_flight")]
    fn evaluate_outliers(&self, hdbscan_result: &mut HDBSCANResult) -> Result<()> {
        let mut outliers = std::mem::take(&mut hdbscan_result.outliers);
        let hdbscan_result_of_filtered_contigs = self.evaluate_subset(
            &mut outliers,
            2,
            5,
            self.n_neighbours,
            DEFAULT_B
        )?;

        debug!("New HDBSCAN score {}", hdbscan_result_of_filtered_contigs.score);
        debug!("Number of clusters: {}", hdbscan_result_of_filtered_contigs.cluster_map.len());
        hdbscan_result.merge(hdbscan_result_of_filtered_contigs);

        Ok(())
    }

    #[cfg(not(feature = "no_flight"))]
    fn run_flight(&self) -> Result<()> {
        check_for_flight()?;

        // if rosella_bins.json exists just skip
        if std::path::Path::new(&format!("{}/rosella_bins.json", &self.output_directory)).exists() {
            return Ok(());
        }

        let mut flight_cmd = Command::new("flight");
        flight_cmd.arg("bin");
        // flight_cmd.arg("--assembly").arg(&self.assembly);
        flight_cmd.arg("--input").arg(format!("{}/coverage.tsv", &self.output_directory));
        flight_cmd.arg("--kmer_frequencies").arg(format!("{}/kmer_frequencies.tsv", &self.output_directory));
        flight_cmd.arg("--output").arg(&self.output_directory);
        flight_cmd.arg("--min_contig_size").arg(format!("{}", self.min_contig_size));
        flight_cmd.arg("--min_bin_size").arg(format!("{}", self.min_bin_size));
        flight_cmd.arg("--n_neighbors").arg(format!("{}", self.n_neighbours));
        flight_cmd.arg("--cores").arg(format!("{}", rayon::current_num_threads()));
        flight_cmd.stdout(std::process::Stdio::piped());
        flight_cmd.stderr(std::process::Stdio::piped());

        let mut child = match flight_cmd.spawn() {
            Ok(child) => child,
            Err(e) => {
                bail!("Error running flight: {}", e);
            }
        };
        
        if let Some(stderr) = child.stderr.take() {
            let stderr = std::io::BufReader::new(stderr);
            for line in stderr.lines() {
                let line = line?;
                let message = line.split("INFO: ").collect::<Vec<_>>();
                info!("{}", message[message.len() - 1]);
            }
        }

        Ok(())
    }

    #[cfg(not(feature = "no_flight"))]
    fn finalise_bins(&mut self) -> Result<Vec<ClusterResult>> {
        let mut bins =
            std::fs::File::open(format!("{}/rosella_bins.json", &self.output_directory))?;
        let mut data = String::new();
        bins.read_to_string(&mut data).unwrap();
        let mut cluster_map: HashMap<usize, HashSet<usize>> = serde_json::from_str(&data).unwrap();

        let mut removed_bins = Vec::new();
        for (bin, contigs) in cluster_map.iter() {
            let mut contigs = contigs.iter().cloned().collect::<Vec<_>>();
            contigs.sort_unstable();
            if *bin == 23 {
                debug!("Bin {} has {:?} contigs", bin, contigs);
                for contig in contigs.iter() {
                    debug!("Contig {} has length {}", self.coverage_table.contig_names[*contig], self.coverage_table.contig_lengths[*contig]);
                }
            }

            let bin_size = contigs.iter().map(|i| self.coverage_table.contig_lengths[*i]).sum::<usize>();
            if bin_size < self.min_bin_size || *bin == 0 {
                removed_bins.push(*bin);
            }
        }

        let mut removed_contigs = HashSet::new();
        for bin in removed_bins {
            match cluster_map.remove(&bin) {
                Some(contigs) => {
                    removed_contigs.par_extend(contigs);
                },
                None => {
                    bail!("Bin {} does not exist in cluster map", bin);
                }
            }
        }

        self.filtered_contigs.par_extend(self.coverage_table.get_contig_names(&removed_contigs));
        // self.contig_sketches.filter_by_index(&removed_contigs)?;
        // self.tnf_table.filter_by_index(&removed_contigs)?;
        // self.n_contigs = self.contig_sketches.contig_names.len();

        Ok(self.get_cluster_result(cluster_map, removed_contigs, None, None))
    }

    #[cfg(feature = "no_flight")]
    /// Embed and cluster a subset of contigs
    /// original_contig_indices are the indices of the contigs in the original contig list after intial filtering
    fn evaluate_subset(&self, original_contig_indices: &mut HashSet<usize>, min_cluster_size: usize, min_sample_size: usize, keep_n_edges: usize, b: f64) -> Result<HDBSCANResult> {
        
        // positional map, showing position index to original contig index
        let contig_id_map = original_contig_indices.iter().enumerate().map(|(i, &j)| (i, j)).collect::<HashMap<_, _>>();
        // let contig_hashset = original_contig_indices.iter().cloned().collect::<HashSet<_>>();
        // retrieve the data corresponding to these contigs via reference
        let (subset_kgraph, disconnected_nodes) = self.build_mutual_kgraph(keep_n_edges, &original_contig_indices)?;
        let (subset_embeddings, _) = self.embed(&subset_kgraph, b, MAX_ITERATIONS - 1)?;

        let mut hdbscan_result = find_best_clusters(&subset_embeddings, min_cluster_size, min_sample_size)?;
        debug!("HDBSCAN score {}", hdbscan_result.score);
        let disconnected_indices = disconnected_nodes.into_iter()
            .map(|i| subset_kgraph.get_data_id_from_idx(i).unwrap())
            .map(|graph_index| *contig_id_map.get(graph_index).unwrap()).collect::<Vec<_>>();

        hdbscan_result.reindex_clusters(contig_id_map);

        hdbscan_result.outliers.par_extend(disconnected_indices);
        let n_clustered_contigs = hdbscan_result.cluster_map.values().map(|v| v.len()).sum::<usize>() + hdbscan_result.outliers.len();

        if n_clustered_contigs != original_contig_indices.len() {
            return Err(anyhow!("Number of clustered contigs does not match number of contigs in subset. {} != {}", n_clustered_contigs, original_contig_indices.len()));
        }

        Ok(hdbscan_result)
    }

    #[cfg(feature = "no_flight")]
    fn embed(&self, kgraph: &KGraph<f64>, b: f64, iteration: usize) -> Result<(Array2<f64>, Option<KGraph<f64>>)> {
        let embedder_params = self.generate_embedder_params(b);
        let mut embedder = Embedder::new(&kgraph, embedder_params);
        let mut optional_new_kgraph = None;
        match embedder.embed() {
            Ok(_) => (),
            Err(e) => {
                debug!("Error embedding: {}", e);
                if iteration < MAX_ITERATIONS {
                    debug!("Embedding failed, generating new kgraph with more edges");
                    let indices_to_include = (0..self.n_contigs).collect::<HashSet<usize>>();
                    let (new_kgraph, _) = self.build_mutual_kgraph(iteration * 10 + 10, &indices_to_include)?;
                    optional_new_kgraph = Some(new_kgraph);
                    return self.embed(optional_new_kgraph.as_ref().unwrap(), b, iteration + 1);
                } else {
                    return Err(anyhow!("Error embedding: {}", e));
                }
            }
        }
        let embeddings = embedder.get_embedded_reindexed();
        Ok((embeddings, optional_new_kgraph))
    }

    #[cfg(feature = "no_flight")]
    /// use Hnsw and annembed to iteratively build mutual K-nn graphs
    /// until we reach a point where all nodes have at leat 1 neighbour
    fn prefilter_nodes(&mut self, rounds: usize) -> Result<KGraph<f64>> {
        let mut indices_to_include = (0..self.n_contigs).collect::<HashSet<usize>>();
        let (mut initial_kgraph, mut disconnected_nodes) = self.build_mutual_kgraph(0, &indices_to_include)?;
        let mut current_round = 1;
        let mut keep_n_edges = 0;
        while !disconnected_nodes.is_empty() {
            info!("Round {} - Filtering {} contigs.", current_round, disconnected_nodes.len());
            self.filter_contigs(&initial_kgraph, disconnected_nodes)?;
            indices_to_include = (0..self.n_contigs).collect::<HashSet<usize>>();
            current_round += 1;
            if current_round > rounds {
                // we've hit the max number of rounds, so ensure the next graph has
                // no disconnected nodes
                keep_n_edges += self.n_neighbours;
            }
            (initial_kgraph, disconnected_nodes) = self.build_mutual_kgraph(keep_n_edges, &indices_to_include)?;
        }
        Ok(initial_kgraph)
    }

    #[cfg(feature = "no_flight")]
    fn filter_contigs(&mut self, kgraph: &KGraph<f64>, disconnected_nodes: Vec<usize>) -> Result<()> {
        // use the kgraph to take the disconnected_contig indices and get the original indices
        // then filter the coverage table and contig sketches by index
        let original_indices = disconnected_nodes
            .iter()
            .filter_map(|idx| {
                kgraph.get_data_id_from_idx(*idx)
            })
            .map(|idx| {
                *idx
            })
            .collect::<HashSet<_>>();
        self.filtered_contigs.extend(self.coverage_table.filter_by_index(&original_indices)?);
        debug!("Filtered contigs: {}", self.filtered_contigs.len());
        self.contig_sketches.filter_by_index(&original_indices)?;
        self.tnf_table.filter_by_index(&original_indices)?;
        self.n_contigs = self.contig_sketches.contig_names.len();

        Ok(())
    }

    #[cfg(feature = "no_flight")]
    /// Build a mutual K-NN graph, but keep at least `keep_n_edges` edges per node
    /// Setting `keep_n_edges` to 0 can result in some nodes becoming disconnected, these nodes
    /// should either be removed or reconnected to the graph
    fn build_mutual_kgraph(&self, keep_n_edges: usize, indices_to_include: &HashSet<usize>) -> Result<(KGraph<f64>, Vec<usize>)> {
        let contig_data = self.retrieve_contig_data(indices_to_include);
        let embedder_engine = EmbedderEngine::new(contig_data);
        let n_neighbours = if indices_to_include.len() < self.n_neighbours * 10 {
            std::cmp::min(indices_to_include.len() / 2, self.n_neighbours)
        } else {
            self.n_neighbours
        };
        embedder_engine.build_mutual_kgraph(keep_n_edges, n_neighbours, self.n_contigs, self.max_layers, self.ef_construction)
    }

    #[cfg(feature = "no_flight")]
    fn retrieve_contig_data<'a>(&'a self, indices_to_include: &HashSet<usize>) -> Vec<Vec<ContigInformation<'a, Dim<[usize; 1]>>>> {
        let contig_data = self.coverage_table.table
            .rows()
            .into_iter()
            .zip(self.tnf_table.kmer_table.rows().into_iter())
            .zip(self.contig_sketches.contig_sketches.iter())
            .enumerate()
            .filter_map(|(row_index, ((coverage_row, tnf_row), sketch))| {
                if !indices_to_include.contains(&row_index) {
                    return None;
                }

                let contig_information = ContigInformation::new(
                    coverage_row,
                    tnf_row,
                    sketch,
                );
                Some(vec![contig_information])
            }).collect::<Vec<_>>();
        
        return contig_data;
    }

    #[cfg(feature = "no_flight")]
    fn generate_embedder_params(&self, b: f64) -> EmbedderParams {
        let mut params = EmbedderParams::default();
        let n_components = min(max(self.coverage_table.table.ncols() / 2, 2), 10);
        params.set_dim(n_components);
        // params.nb_grad_batch = self.nb_grad_batches;
        params.scale_rho = 1.0;
        params.beta = 0.5;
        // params.b = 1.0;
        params.b = b;

        params.grad_step = 2.;
        params.nb_sampling_by_edge = 10;
        params.nb_grad_batch = 15;
        params.grad_factor = 4;
        params.dmap_init = true;

        params
    }

    fn get_cluster_result(
        &self,
        cluster_map: HashMap<usize, HashSet<usize>>,
        outliers: HashSet<usize>,
        embeddings: Option<&ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>>,
        contig_index_map: Option<HashMap<usize, usize>>, // a map containing the key as the row index in the embeddings array and the value as the original contig index
                                                         // used when a cluster has been subset and re-embedded
    ) -> Vec<ClusterResult> {
        if DEBUG_BINS {
            // delete bins.txt
            let _ = std::fs::remove_file(format!("{}/bins.txt", self.output_directory));
        }

        let mut cluster_results = Vec::with_capacity(self.n_contigs);
        for (cluster_label, mut contig_indices) in cluster_map.into_iter() {
            contig_indices = match &contig_index_map {
                Some(index_map) => {
                    contig_indices
                        .into_iter()
                        .filter_map(|idx| {
                            index_map.get(&idx).map(|idx| {
                                *idx
                            })
                        })
                        .collect::<HashSet<_>>()
                },
                None => contig_indices,
            };
            // check the size of the cluster and if it is too small, set the cluster label to None
            let bin_size = contig_indices.iter().map(|i| self.coverage_table.contig_lengths[*i]).sum::<usize>();
            let cluster_label = if bin_size < self.min_bin_size {
                None
            } else {
                Some(cluster_label)
            };
            for contig_index in contig_indices.iter() {
                cluster_results.push(ClusterResult::new(*contig_index, cluster_label));
            }

            #[cfg(feature = "no_flight")]
            if DEBUG_BINS && embeddings.is_some() {
                let file = OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(format!("{}/bins.txt", self.output_directory))
                    .unwrap();
                let mut writer = BufWriter::new(file);
                let embeddings = embeddings.unwrap();
                
                
                debug!("Calculating metrics for bin {:?} with {} contigs", cluster_label, contig_indices.len());
                let bin_metrics = self.calculate_bin_metrics(&contig_indices, embeddings);
                debug!("Bin {:?} n_contigs {} metrics: {:?}", cluster_label, contig_indices.len(), bin_metrics);
                writeln!(writer, "{:?}\t{}\t{}\t{}\t{}\t{}\t{}", 
                    cluster_label, contig_indices.len(), bin_metrics.bin_size, bin_metrics.metabat_dist, bin_metrics.tnf_distance, bin_metrics.sketch_dist, bin_metrics.embedding_dist).unwrap();
                
            }
        }
        for outlier in outliers {
            cluster_results.push(ClusterResult::new(outlier, None));
        }
        cluster_results.par_sort_unstable();

        debug!("Cluster results: {:?}", &cluster_results[0..10]);
        cluster_results
    }

    #[cfg(feature = "no_flight")]
    /// calculate average distance between contigs in a bin
    fn calculate_bin_metrics(&self, contig_indices: &HashSet<usize>, embeddings: &ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>) -> BinMetrics {
        let bin_size = contig_indices.iter().map(|i| self.coverage_table.contig_lengths[*i]).sum::<usize>();
        
        // calculate all pairwise metabat distances
        let n_combinations = (contig_indices.len() * (contig_indices.len() - 1)) / 2;
        let indices = contig_indices.iter().cloned().collect::<Vec<_>>();
        let distances = (0..contig_indices.len())
            .into_par_iter()
            .flat_map(|i| {

                let contig_i = indices[i];
                let va_cov = self.coverage_table.table.row(contig_i);
                let va_sketch = &self.contig_sketches.contig_sketches[contig_i];
                let va_embedding = embeddings.row(contig_i);
                (i+1..contig_indices.len())
                    .into_par_iter()
                    .map(|j| {
                        let contig_j = &indices[j];
                        let vb_cov = self.coverage_table.table.row(*contig_j);
                        let vb_sketch = &self.contig_sketches.contig_sketches[*contig_j];
                        let vb_embedding = embeddings.row(*contig_j);
                        let mut metabat_dist = MetabatDistance::distance(va_cov, vb_cov);
                        if metabat_dist.is_nan() {
                            metabat_dist = 1.;
                        }

                        let mut tnf_distance = KmerCorrelation::distance(va_cov, vb_cov);
                        if tnf_distance.is_nan() {
                            tnf_distance = 1.;
                        }

                        let mut sketch_dist = 1. - distance(va_sketch, vb_sketch).unwrap().min_jaccard;
                        if sketch_dist.is_nan() {
                            sketch_dist = 1.;
                        }

                        // euclidean distance between embeddings
                        let mut embedding_dist = 0.;
                        for (a, b) in va_embedding.iter().zip(vb_embedding.iter()) {
                            embedding_dist += (a - b).powi(2);
                        }
                        embedding_dist = embedding_dist.sqrt();
                        if embedding_dist.is_nan() {
                            embedding_dist = 1.;
                        }

                        (metabat_dist as f32, tnf_distance as f32, sketch_dist as f32, embedding_dist as f32)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let metabat_dist = distances.iter().map(|(metabat_dist, _, _, _)| *metabat_dist).sum::<f32>() / n_combinations as f32;
        let tnf_dist = distances.iter().map(|(_, tnf_dist, _, _)| *tnf_dist).sum::<f32>() / n_combinations as f32;
        let sketch_dist = distances.iter().map(|(_, _, sketch_dist, _)| *sketch_dist).sum::<f32>() / n_combinations as f32;
        let embedding_dist = distances.iter().map(|(_, _, _, embedding_dist)| *embedding_dist).sum::<f32>() / n_combinations as f32;
        BinMetrics { bin_size, tnf_distance: tnf_dist, metabat_dist, sketch_dist, embedding_dist }
    }

    /// Take the cluster results and collect the contigs into bins
    fn write_clusters(&self, cluster_results: Vec<ClusterResult>) -> Result<()> {
        // open the assembly
        let mut reader = parse_fastx_file(path::Path::new(&self.assembly))?;

        let mut contig_idx = 0;
        let mut single_contig_bin_id = 0;
        let mut total_contig_idx = 0;
        debug!("Cluster result length: {}", cluster_results.len());
        while let Some(record) = reader.next() {
            let seqrec = record?;
            
            let contig_name = seqrec.id();
            let contig_name_str = std::str::from_utf8(contig_name)?;
            let contig_length = seqrec.seq().len();
            if contig_length < self.min_contig_size {
                let cluster_label = if seqrec.seq().len() < self.min_bin_size {
                    UNBINNED.to_string()
                } else {
                    single_contig_bin_id += 1;
                    #[cfg(not(feature = "no_flight"))]
                    {
                        contig_idx += 1;
                    }
                    format!("single_contig_{}", single_contig_bin_id)
                };

                let bin_path = path::Path::new(&self.output_directory).join(format!("rosella_bin_{}{}", cluster_label, RECOVER_FASTA_EXTENSION));
                let file = OpenOptions::new().append(true).create(true).open(bin_path)?;
                let mut writer = BufWriter::new(file);
                // write contig to bin
                write_fasta(contig_name, &seqrec.seq(), &mut writer, LineEnding::Unix)?;
                total_contig_idx += 1;
                continue;
            }
            // normalise sequence
            let cluster_result = &cluster_results[contig_idx];

            if total_contig_idx != cluster_result.contig_index {
                debug!("Contig index mismatch. {} != {}", total_contig_idx, cluster_result.contig_index);
            }
            
            // open writer for bin
            let cluster_label = match cluster_result.cluster_label {
                Some(cluster_label) => format!("{}", cluster_label),
                None => {
                    // contig is an outlier, so check it's length. If it is greater than
                    // the minimum bin size, then bin it, otherwise discard it.
                    if seqrec.seq().len() < self.min_bin_size {
                        UNBINNED.to_string()
                    } else {
                        single_contig_bin_id += 1;
                        format!("single_contig_{}", single_contig_bin_id)
                    }

                }
            };

            let bin_path = path::Path::new(&self.output_directory).join(format!("rosella_bin_{}.fna", cluster_label));
            let file = OpenOptions::new().append(true).create(true).open(bin_path)?;
            let mut writer = BufWriter::new(file);
            // write contig to bin
            write_fasta(contig_name, &seqrec.seq(), &mut writer, LineEnding::Unix)?;
            contig_idx += 1;
            total_contig_idx += 1;
        }

        debug!("Final contig idx: {}", contig_idx);

        Ok(())
    }
}

#[derive(Debug)]
struct BinMetrics {
    pub(crate) bin_size: usize,
    pub(crate) metabat_dist: f32,
    pub(crate) tnf_distance: f32,
    pub(crate) sketch_dist: f32,
    pub(crate) embedding_dist: f32,
}

#[derive(Debug, Eq, PartialEq, PartialOrd)]
struct ClusterResult {
    contig_index: usize,
    cluster_label: Option<usize>,
}

impl ClusterResult {
    pub fn new(contig_index: usize, cluster_label: Option<usize>) -> Self {
        Self {
            contig_index,
            cluster_label,
        }
    }
}

impl Ord for ClusterResult {
    fn cmp(&self, other: &Self) -> Ordering {
        self.contig_index.cmp(&other.contig_index)
    }
}