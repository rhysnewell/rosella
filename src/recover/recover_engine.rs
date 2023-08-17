use std::{collections::{HashMap, HashSet, VecDeque}, cmp::{Ordering, max, min}, path, io::{BufWriter, Write, BufRead, Read}, fs::OpenOptions, process::Command};

use annembed::{prelude::*, fromhnsw::kgraph::KGraph};
use anyhow::Result;
use bird_tool_utils::external_command_checker;
use log::{info, debug, error};
use ndarray::{Dim, ArrayBase, OwnedRepr, Array2};
use needletail::{parse_fastx_file, parser::{write_fasta, LineEnding}};
use petgraph::prelude::*;
use rayon::{prelude::*, slice::ParallelSliceMut, prelude::IntoParallelIterator};


use crate::{
    coverage::{coverage_calculator::{calculate_coverage, MetabatDistance}, coverage_table::CoverageTable}, 
    sketch::{contig_sketcher::{ContigSketchResult, sketch_contigs}, sketch_distances::distance}, 
    embedding::embedder::{ContigInformation, EmbedderEngine}, 
    clustering::{clusterer::{find_best_clusters, HDBSCANResult, build_kgraph_of_clusters}, PropagatedLabel}, 
    kmers::kmer_counting::{KmerFrequencyTable, count_kmers, KmerCorrelation}, graphs::nearest_neighbour_graph::intersect, external_command_checker::check_for_flight
};

const RECOVER_FASTA_EXTENSION: &str = ".fna";
const DEBUG_BINS: bool = true;
const UNBINNED: &str = "unbinned";
const INITIAL_SMALL_BIN_SIZE: usize = 100000;
const DEFAULT_B: f64 = 1.5;
const INITIAL_KEPT_BIN_SIZE: usize = 1000000;
const RECLUSTERING_MIN_BIN_SIZE: usize = 25000;
const MAXIMUM_SENSIBLE_BIN_SIZE: usize = 12000000;
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
    contig_sketches: ContigSketchResult,
    tnf_table: KmerFrequencyTable,
    n_neighbours: usize,
    filtering_rounds: usize,
    ef_construction: usize,
    max_layers: usize,
    n_contigs: usize,
    nb_grad_batches: usize,
    min_bin_size: usize,
    min_contig_size: usize,
    embeddings: Option<Array2<f64>>,
    filtered_contigs: HashSet<String>,
    // contig_nn: Hnsw<ContigInformation<'static, D>, ContigDistance>,
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
        let filtered_contigs = coverage_table.filter_by_length(min_contig_size)?;
        assert_eq!(coverage_table.table.nrows(), n_contigs - filtered_contigs.len(), "Coverage table row count and total contigs minus filtered contigs do not match.");
        info!("Calculating TNF table.");
        let mut tnf_table = count_kmers(m, Some(n_contigs))?;
        assert_eq!(n_contigs, tnf_table.kmer_table.nrows(), "Coverage table row count and TNF table row count do not match.");
        debug!("Filtering TNF table.");
        tnf_table.filter_by_name(&filtered_contigs)?;
        assert_eq!(coverage_table.table.nrows(), tnf_table.kmer_table.nrows(), "Coverage table and TNF table have different number of contigs.");
        info!("Calculating contig sketches.");
        let contig_sketches = sketch_contigs(m)?;
        assert_eq!(coverage_table.table.nrows(), contig_sketches.contig_sketches.len(), "Coverage table and contig sketches have different number of contigs.");
        info!("{} valid contigs, {} filtered contigs.", coverage_table.table.nrows(), filtered_contigs.len());
        let n_neighbours = m.get_one::<usize>("n-neighbours").unwrap().clone();
        let filtering_rounds = m.get_one::<usize>("filtering-rounds").unwrap().clone();
        let ef_construction = m.get_one::<usize>("ef-construction").unwrap().clone();
        let max_layers = m.get_one::<usize>("max-layers").unwrap().clone();
        let nb_grad_batches = m.get_one::<usize>("nb-grad-batches").unwrap().clone();
        let min_bin_size = m.get_one::<usize>("min-bin-size").unwrap().clone();
        
        let n_contigs = contig_sketches.contig_names.len();
        Ok(
            Self {
                output_directory,
                assembly,
                coverage_table,
                contig_sketches,
                tnf_table,
                n_neighbours,
                filtering_rounds,
                ef_construction,
                max_layers,
                n_contigs,
                nb_grad_batches,
                min_bin_size,
                min_contig_size,
                embeddings: None,
                // filtered_contigs,
                filtered_contigs: HashSet::new(),
            }
        )
    }

    /// Runs through the rosella bin recovery pipeline
    pub fn run(&mut self) -> Result<()> {
        // 1. Perform N rounds of filtering. Filtering involves removinf disconnnected contigs from the graph.
        //    The final round ensures each contig has at least one neighbour.
        info!("Filtering.");
        let kgraph = self.prefilter_nodes(self.filtering_rounds)?;
        // let kgraph = self.generate_multi_kgraph(3, 0, 0)?;
        
        // 2. Take the KGraph and embed using a UMAP-esque method.
        //    Embedding parameters need to be chosen more carefully.
        info!("Embedding.");
        let (embeddings, new_kgraph) = self.embed(&kgraph, DEFAULT_B, 0)?;
        let kgraph = if let Some(new_kgraph) = new_kgraph {
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

        // 4. Then we need to inspect each cluster individually and see if we can split it into smaller clusters.
        //    This involves embedding the contigs in each cluster and then clustering them again. If the resulting cluster looks
        //    decent we can keep it.
        
        // for _ in 0..self.filtering_rounds {
        // }
        
        // 5. After the initial embedding and clustering, this is where things get tricky. We need to filter out unclustered contigs
        //    or clusters that we think might be too small and then re-embed and re-cluster.
        // self.find_close_clusters(&mut hdbscan_result, &kgraph, false, false, true)?;
        // debug!("HDBSCAN outlier percentage: {}", hdbscan_result.outliers.len() as f64 / self.n_contigs as f64);
        // self.reembed_clusters(&mut hdbscan_result, &kgraph, &embeddings)?;
        // debug!("HDBSCAN outlier percentage: {}", hdbscan_result.outliers.len() as f64 / self.n_contigs as f64);

        info!("Rescuing unbinned.");
        self.evaluate_outliers(&mut hdbscan_result)?;
        // self.find_close_clusters(&mut hdbscan_result, &kgraph, false, false, true)?;

        // info!("Re-clustering.");
        // for _ in 0..self.filtering_rounds {
        //     self.find_close_clusters(&mut hdbscan_result, &kgraph, true, true, false)?;
        // }
        
        // for _ in 0..self.filtering_rounds {
        //     // self.find_close_clusters(&mut hdbscan_result, &kgraph, false, true, false)?;
        // }
        info!("HDBSCAN outlier percentage: {}", hdbscan_result.outliers.len() as f64 / self.n_contigs as f64);
        
        let n_contigs = hdbscan_result.cluster_map.values().map(|v| v.len()).sum::<usize>() + hdbscan_result.outliers.len();
        let cluster_results = self.get_cluster_result(hdbscan_result.cluster_map, hdbscan_result.outliers, Some(&embeddings), None);
        info!("Length of cluster results: {}", cluster_results.len());
        debug!("cluster result len {} n_contigs {} contigs used {} coverage table and kmer table size {} {} n filtered contigs {}", 
            cluster_results.len(), n_contigs, self.n_contigs, self.coverage_table.contig_lengths.len(), self.tnf_table.contig_names.len(), self.filtered_contigs.len());
        if n_contigs != cluster_results.len() {
            bail!("Number of contigs in cluster results ({}) does not match number of contigs in HDBSCAN result ({})", cluster_results.len(), n_contigs);
        }

        // self.run_flight()?;
        // let cluster_results = self.finalise_bins()?;
        
        info!("Writing clusters.");
        self.write_clusters(cluster_results)?;
        Ok(())
    }

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

    fn calculate_cluster_connectivity(&self, cluster: &HashSet<usize>, kgraph: &KGraph<f64>) -> Result<f64> {
        let mut total_connectivity = 0.0;
        let mut cluster_connectivity = 0.0;
        let mut n_contigs_connected_in_cluster = 0;
        let mut n_edges = 0;
        let cluster_hashset = cluster
            .iter()
            .map(|contig_index| kgraph.get_idx_from_dataid(contig_index).unwrap())
            .collect::<HashSet<usize>>();
        let kgraph_neighbours = kgraph.get_neighbours();

        for contig_i in cluster {
            let graph_index_i = kgraph.get_idx_from_dataid(&contig_i).unwrap();
            let neighbours_of_i = &kgraph_neighbours[graph_index_i];
            for neighbour in neighbours_of_i {
                if cluster_hashset.contains(&neighbour.node) {
                    cluster_connectivity += neighbour.weight;
                    n_contigs_connected_in_cluster += 1;
                    // break;
                }
                total_connectivity += neighbour.weight;
                n_edges += 1;
            }
        }
        Ok(n_contigs_connected_in_cluster as f64 / n_edges as f64)
    }

    fn reembed_clusters(
        &self,
        hdbscan_result: &mut HDBSCANResult,
        kgraph: &KGraph<f64>,
        embeddings: &ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>,
    ) -> Result<()> {
        let original_n_contigs = hdbscan_result.cluster_map.values().map(|v| v.len()).sum::<usize>() + hdbscan_result.outliers.len();
        hdbscan_result.renumber_clusters();

        // recluster every cluster
        let new_clusters = hdbscan_result.cluster_map.par_iter_mut()
            .filter_map(|(original_cluster_id, contigs)| {
                let n_contigs = contigs.len();

                let cluster_connectivity = self.calculate_cluster_connectivity(contigs, kgraph).unwrap();
                let cluster_metrics = self.calculate_bin_metrics(contigs, embeddings);
                debug!("Cluster {} connectivity {} n_contigs {}", original_cluster_id, cluster_connectivity, n_contigs);
                // we can filter here if we want to skip some clusters
                if n_contigs < 15 || cluster_metrics.metabat_dist < 0.05 {
                    // this will keep the cluster in the map but won't reembed it
                    return None;
                }

                let cluster_size = contigs.iter().map(|c| self.coverage_table.contig_lengths[*c]).sum::<usize>();
                match self.evaluate_subset(contigs, 2, 5, n_contigs, DEFAULT_B) {
                    Ok(new_cluster) => {
                        debug!("Cluster {} silhouette score {}", original_cluster_id, new_cluster.score);
                        // if new_cluster.score < Self::connectivity_score_to_silhouette_score_threshold(cluster_connectivity) 
                        if new_cluster.score < 0.95
                            && !(cluster_size > MAXIMUM_SENSIBLE_BIN_SIZE) {
                            return None;
                        }
                        Some((*original_cluster_id, new_cluster))
                    },
                    Err(e) => {
                        debug!("Error re-embedding cluster: {}", e);
                        None
                    }
                }

            })
            .collect::<HashMap<_, _>>();
        
        // hdbscan_result.cluster_map = HashMap::with_capacity(hdbscan_result.cluster_map.len() * 10);
        let mut max_cluster_id = hdbscan_result.cluster_map.keys().max().unwrap().clone() + 1;
        debug!("Max cluster id {}", max_cluster_id);
        new_clusters.into_iter()
            .for_each(|(original_cluster_id, new_clusters)| {
                // new_clusters.renumber_clusters();
                let original_cluster = hdbscan_result.cluster_map.remove(&original_cluster_id).unwrap();
                debug!("Original cluster {} silhouette score {} n_clusters {}", original_cluster_id, new_clusters.score, new_clusters.cluster_map.len());
                debug!("Original length {}", original_cluster.len());
                let mut sum_of_cluster_lengths = 0;
                for (_, cluster) in new_clusters.cluster_map {
                    debug!("New cluster {} length {}", max_cluster_id, cluster.len());
                    sum_of_cluster_lengths += cluster.len();
                    hdbscan_result.cluster_map.insert(max_cluster_id, cluster);
                    max_cluster_id += 1;
                }
                sum_of_cluster_lengths += new_clusters.outliers.len();
                assert_eq!(original_cluster.len(), sum_of_cluster_lengths, 
                    "Number of contigs changed after re-embedding clusters for cluster {}: {} -> {}", original_cluster_id, original_cluster.len(), sum_of_cluster_lengths);
                hdbscan_result.outliers.par_extend(new_clusters.outliers);
            });
        
        let new_n_contigs = hdbscan_result.cluster_map.values().map(|v| v.len()).sum::<usize>() + hdbscan_result.outliers.len();
        assert_eq!(original_n_contigs, new_n_contigs, "Number of contigs changed after re-embedding clusters. {} -> {}", original_n_contigs, new_n_contigs);
        Ok(())
    }

    fn connectivity_score_to_silhouette_score_threshold(connectivity_score: f64) -> f64 {
        // this is a linear function that maps the connectivity score to a silhouette score threshold
        0.35 * connectivity_score + 0.5
    }

    /// Now we are going to the cluster result and compare it to the original KGraph. We have two aims in this algorithm:
    ///     1. Try and assign unclustered contigs to their closest cluster if the answer is obvious
    ///     2. Try and combine small clusters that share a lot of edges in the KGraph
    /// The KGraph only contains the top edges between contigs, so it is not a complete graph. 
    fn find_close_clusters(&self, hdbscan_result: &mut HDBSCANResult, kgraph: &KGraph<f64>, rescue_outliers: bool, combine_clusters: bool, trim_bins: bool) -> Result<()> {
        let n_contigs_in_result = hdbscan_result.cluster_map.values().map(|v| v.len()).sum::<usize>() + hdbscan_result.outliers.len();

        hdbscan_result.renumber_clusters();
        let contig_to_cluster_map = hdbscan_result.get_contig_to_cluster_map();
        debug!("Number of clusters: {}", contig_to_cluster_map.len());
        // step 1. trim bins that have disconnected contigs
        if trim_bins {
            self.trim_bins(hdbscan_result, kgraph, &contig_to_cluster_map)?;
        }

        // step 2. find the closest cluster for each outlier
        if rescue_outliers {
            self.rescue_outliers(hdbscan_result, kgraph, &contig_to_cluster_map)?;
        }

        // step 3. find clusters that share a lot of edges and combine them
        //         we'll limit this to clusters that are smaller than the minimum bin size
        //         but thay can be combined into larger bins
        if combine_clusters {
            self.combine_clusters(hdbscan_result, kgraph, &contig_to_cluster_map)?;
        }

        

        let new_n_contigs_in_result = hdbscan_result.cluster_map.values().map(|v| v.len()).sum::<usize>() + hdbscan_result.outliers.len();
        if n_contigs_in_result != new_n_contigs_in_result {
            bail!("Number of contigs changed after adding outliers back to clusters. {} -> {}", n_contigs_in_result, new_n_contigs_in_result);
        }

        Ok(())
    }

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
                debug!("{}", line);
            }
        }

        Ok(())
    }

    fn finalise_bins(&mut self) -> Result<Vec<ClusterResult>> {
        let mut bins =
            std::fs::File::open(format!("{}/rosella_bins.json", &self.output_directory))?;
        let mut data = String::new();
        bins.read_to_string(&mut data).unwrap();
        let mut cluster_map: HashMap<usize, HashSet<usize>> = serde_json::from_str(&data).unwrap();

        let mut removed_bins = Vec::new();
        for (bin, contigs) in cluster_map.iter() {

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

    fn rescue_outliers(&self, hdbscan_result: &mut HDBSCANResult, kgraph: &KGraph<f64>, contig_to_cluster_map: &HashMap<usize, usize>) -> Result<()> {
        let n_contigs_in_result = hdbscan_result.cluster_map.values().map(|v| v.len()).sum::<usize>() + hdbscan_result.outliers.len();

        let outliers = std::mem::take(&mut hdbscan_result.outliers);
        // search the nearest neighbours of each outlier and record which cluster it belongs to
        let outlier_cluster_probabilities = outliers.par_iter()
            .map(|outlier| {
                let mut cluster_probabilities = HashMap::new();
                let neighbours = kgraph.get_out_edges_by_data_id(outlier).unwrap();
                let mut total_weight = 0.0;
                for edge in neighbours {
                    // take edge.node and find the contig index
                    let true_index = kgraph.get_data_id_from_idx(edge.node).unwrap();
                    match contig_to_cluster_map.get(&true_index) {
                        Some(cluster_id) => {
                            let cluster_weight = cluster_probabilities.entry(*cluster_id).or_insert(0.0);
                            *cluster_weight += edge.weight;
                            total_weight += edge.weight;
                        },
                        None => {
                            // this is an unclustered contig, we can't do anything with it
                        }
                    };
                }
                
                // now we have a map of cluster ids to weights, we can find the most likely cluster
                let mut max_cluster_id = None;
                let mut max_cluster_weight = 0.0;
                for (cluster_id, cluster_weight) in cluster_probabilities {
                    if cluster_weight > max_cluster_weight {
                        max_cluster_id = Some(cluster_id);
                        max_cluster_weight = cluster_weight;
                    }
                }

                max_cluster_weight /= total_weight;
                if max_cluster_weight < 0.5 {
                    max_cluster_id = None;
                }

                PropagatedLabel {
                    label: max_cluster_id,
                    weight: max_cluster_weight,
                    index: *outlier
                }

            })
            .collect::<Vec<_>>();

        // now we have a list of outlier cluster assignments, we can assign them to the clusters
        for outlier in outlier_cluster_probabilities {
            if let Some(cluster_id) = outlier.label {
                let cluster = hdbscan_result.cluster_map.get_mut(&cluster_id).unwrap();
                cluster.insert(outlier.index);
            } else {
                // this outlier doesn't belong to any cluster, so we add it back to the list of outliers
                hdbscan_result.outliers.insert(outlier.index);
            }
        }

        let new_n_contigs_in_result = hdbscan_result.cluster_map.values().map(|v| v.len()).sum::<usize>() + hdbscan_result.outliers.len();
        if n_contigs_in_result != new_n_contigs_in_result {
            bail!("Number of contigs changed after adding outliers back to clusters. {} -> {}", n_contigs_in_result, new_n_contigs_in_result);
        }

        Ok(())
    }

    fn get_cluster_size_bp(&self, cluster: &HashSet<usize>) -> usize {
        cluster.iter().map(|contig| self.coverage_table.contig_lengths[*contig]).sum()
    }

    fn combine_clusters(&self, hdbscan_result: &mut HDBSCANResult, kgraph: &KGraph<f64>, contig_to_cluster_map: &HashMap<usize, usize>) -> Result<()> {
        let n_contigs_in_result = hdbscan_result.cluster_map.values().map(|v: &HashSet<usize>| v.len()).sum::<usize>() + hdbscan_result.outliers.len();

        let mut cluster_map = std::mem::take(&mut hdbscan_result.cluster_map);
        let propagated_clusters = cluster_map.par_iter().map(|(cluster_index, contigs)| {
            let mut cluster_probabilities = HashMap::new();
            let mut total_connection_weight = 0.0;

            let cluster_size = self.get_cluster_size_bp(contigs);
            for contig in contigs {
                let neighbours = kgraph.get_out_edges_by_data_id(contig).unwrap();
                for edge in neighbours {
                    // take edge.node and find the contig index
                    let true_index = kgraph.get_data_id_from_idx(edge.node).unwrap();
                    match contig_to_cluster_map.get(&true_index) {
                        Some(cluster_id) => {
                            if cluster_index == cluster_id {
                                // this is an edge within the cluster, we don't care about it
                                continue;
                            }
                            let cluster_weight = cluster_probabilities.entry(*cluster_id).or_insert(0.0);
                            *cluster_weight += edge.weight;
                            total_connection_weight += edge.weight;
                        },
                        None => {
                            // this is an unclustered contig, we can't do anything with it
                        }
                    };
                }
            }
            
            // now we have a map of cluster ids to weights, we can find the most likely cluster
            let mut max_cluster_id = None;
            let mut max_cluster_weight = 0.0;
            for (cluster_id, cluster_weight) in cluster_probabilities {
                if cluster_weight > max_cluster_weight {
                    max_cluster_id = Some(cluster_id);
                    max_cluster_weight = cluster_weight;
                }
            }

            max_cluster_weight /= total_connection_weight;

            if (max_cluster_weight < 0.75 && cluster_size > INITIAL_SMALL_BIN_SIZE) || max_cluster_weight < 0.5 || max_cluster_weight.is_nan() {
                max_cluster_id = None;
            } else {
                debug!("Cluster {} has max cluster {} with weight {}", cluster_index, max_cluster_id.unwrap_or(0), max_cluster_weight);
            }

            PropagatedLabel {
                label: max_cluster_id,
                weight: max_cluster_weight,
                index: *cluster_index
            }
        }).collect::<Vec<_>>();

        // Now we have for each cluster a potential new cluster to combine it into
        // we need to combine these clusters, but we need to be aware of the fact that
        // we might be combining two clusters into a cluster that already exists
        // so we will treat each cluster as a graph node, and then find the connected components
        // of this graph. Each connected component will be a new cluster
        let mut graph = Graph::<usize, f64>::new();
        let mut cluster_id_to_node_index = HashMap::new();
        for cluster_link in propagated_clusters {
            if let Some(cluster_id) = cluster_link.label {
                // check if nodes exist in the graph already
                let n1 = cluster_id_to_node_index.entry(cluster_link.index).or_insert_with(|| graph.add_node(cluster_link.index));
                let n1 = *n1;
                let n2 = cluster_id_to_node_index.entry(cluster_id).or_insert_with(|| graph.add_node(cluster_id));
                let n2 = *n2;

                graph.add_edge(n1, n2, cluster_link.weight);
            }
        }

        let connected_components = petgraph::algo::connected_components(&graph);
        debug!("Found {} connected components", connected_components);
        // get the minimum spanning trees
        let mut new_cluster_map = HashMap::new();
        let mut max_cluster_id = cluster_map.keys().max().unwrap() + 1;
        let mut visited_nodes = HashSet::new();

        // visit all the nodes in the graph and find the clusters that they are linked to
        // update the visited nodes set to ensure we do not visit them again
        for current_node in graph.node_indices() {
            if visited_nodes.contains(&current_node) {
                continue;
            }

            // find all connected nodes
            let mut connected_nodes = HashSet::new();
            let mut queue = VecDeque::new();
            queue.push_back(current_node);
            while let Some(node) = queue.pop_front() {
                if visited_nodes.contains(&node) {
                    continue;
                }
                visited_nodes.insert(node);
                connected_nodes.insert(graph[node]);
                for neighbour in graph.neighbors(node) {
                    queue.push_back(neighbour);
                }
            }

            
            let mut cluster = HashSet::new();
            for cluster_id in connected_nodes {
                // remove from original cluster map
                if let Some(contigs) = cluster_map.remove(&cluster_id) {
                    cluster.extend(contigs);
                }
            }

            // now we have a cluster, we need to add it to the new cluster map
            if let Some(_) = new_cluster_map.insert(max_cluster_id, cluster.into_iter().collect::<HashSet<_>>()) {
                bail!("Cluster {} already exists in cluster map", max_cluster_id);
            }
            max_cluster_id += 1;
        }   

        // now combine the new and old cluster map
        cluster_map.extend(new_cluster_map);
        hdbscan_result.cluster_map = cluster_map;
        let new_n_contigs_in_result = hdbscan_result.cluster_map.values().map(|v| v.len()).sum::<usize>() + hdbscan_result.outliers.len();
        if n_contigs_in_result != new_n_contigs_in_result {
            bail!("Number of contigs changed after adding outliers back to clusters. {} -> {}", n_contigs_in_result, new_n_contigs_in_result);
        }

        Ok(())
    }

    /// Uses the KNN graph to find contigs within clusters that are partially disconnected from the rest of the cluster
    fn trim_bins(&self, hdbscan_result: &mut HDBSCANResult, kgraph: &KGraph<f64>, contig_to_cluster_map: &HashMap<usize, usize>) -> Result<()> {
        let n_contigs_in_result = hdbscan_result.cluster_map.values().map(|v: &HashSet<usize>| v.len()).sum::<usize>() + hdbscan_result.outliers.len();

        let mut cluster_map = std::mem::take(&mut hdbscan_result.cluster_map);
        let extra_outliers = cluster_map
            .par_iter_mut()
            .flat_map(|(_, contigs)| {
                
                // find which contigs are outliers
                let cluster_outliers = contigs.iter().filter_map(|contig| {
                    let contig_connectivity = self.calculate_contig_connectivity(*contig, contigs, kgraph);
                    if contig_connectivity < 0.1 {
                        Some(*contig)
                    } else {
                        None
                    }
                }).collect::<Vec<_>>();

                // remove the outliers from the cluster
                *contigs = contigs.iter().filter(|contig| !cluster_outliers.contains(contig)).map(|c| *c).collect::<HashSet<_>>();
                cluster_outliers
            })
            .collect::<Vec<_>>();

        // add the outliers to the outliers list
        hdbscan_result.outliers.par_extend(extra_outliers);
        hdbscan_result.cluster_map = cluster_map;
        let new_n_contigs_in_result = hdbscan_result.cluster_map.values().map(|v| v.len()).sum::<usize>() + hdbscan_result.outliers.len();
        if n_contigs_in_result != new_n_contigs_in_result {
            bail!("Number of contigs changed after adding outliers back to clusters. {} -> {}", n_contigs_in_result, new_n_contigs_in_result);
        }

        Ok(())
    }

    /// compares the connectivity of this contig to contigs within the cluster and the total connectivity of the contig
    /// returns the ratio of the connectivity to the cluster vs the connectivity to the rest of the graph
    fn calculate_contig_connectivity(&self, contig: usize, cluster: &HashSet<usize>, kgraph: &KGraph<f64>) -> f64 {
        let mut cluster_connectivity = 0.0;
        let mut total_connectivity = 0.0;

        let contig_idx = match kgraph.get_idx_from_dataid(&contig) {
            Some(idx) => idx,
            None => return 0.0
        };

        let neighbours = kgraph.get_neighbours();
        let neighbours_of_contig = &neighbours[contig_idx];
        for other_contig_idx in neighbours_of_contig {
            let other_contig = match kgraph.get_data_id_from_idx(other_contig_idx.node) {
                Some(contig) => contig,
                None => continue
            };

            if cluster.contains(&other_contig) {
                cluster_connectivity += other_contig_idx.weight;
            }
            total_connectivity += other_contig_idx.weight;
        }

        cluster_connectivity / total_connectivity
    }

    fn _find_close_clusters(&self, hdbscan_result: &mut HDBSCANResult, max_bin_size: usize, min_bin_size: usize) -> Result<()>  {
        // first we want to reinsert outliers as single contig clusters
        let kept_clusters = hdbscan_result.prepare_for_clustering_of_clusters(max_bin_size, min_bin_size, &self.coverage_table.contig_lengths, true);
        hdbscan_result.renumber_clusters();

        let contig_data_for_clusters = (0..hdbscan_result.cluster_map.len()).into_par_iter()
            .map(|cluster_id| {
                let cluster = hdbscan_result.cluster_map.get(&cluster_id).unwrap();
                let cluster_set = cluster.into_iter().cloned().collect::<HashSet<_>>();
                let cluster_points = self.retrieve_contig_data(&cluster_set);
                (cluster_id, cluster_points)
            })
            .collect::<HashMap<_, _>>();

        let (cluster_kgraph, _) = build_kgraph_of_clusters(
            &contig_data_for_clusters, 
            1, 
            self.n_neighbours, 
            self.max_layers, 
            self.ef_construction
        )?;

        debug!("Embedding clusters.");
        let (cluster_embeddings, _) = self.embed(&cluster_kgraph, DEFAULT_B, MAX_ITERATIONS - 1)?;
        debug!("Clustering clusters.");
        let cluster_hdbscan_result = find_best_clusters(&cluster_embeddings, 2, 5)?;
        debug!("Merging.");
        hdbscan_result.merge_clusters_from_result(cluster_hdbscan_result);
        hdbscan_result.merge_cluster_from_map(kept_clusters);

        Ok(())
    }

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

    /// get the bin metrics for each cluster i.e. bin size, and filter out clusters that are too small
    /// then re-embed and re-cluster
    fn filter_small_clusters(&self, hdbscan_result: &mut HDBSCANResult, filter_size: usize) {
        let filtered_clusters = hdbscan_result
            .cluster_map
            .iter()
            .filter_map(|(cluster, indices)| {
                let bin_size = indices.iter().map(|i| self.coverage_table.contig_lengths[*i]).sum::<usize>();
                if bin_size < filter_size || indices.len() < 2 {
                    Some(*cluster)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        for cluster in filtered_clusters {
            if let Some(filtered_indices) = hdbscan_result.cluster_map.remove(&cluster) {
                hdbscan_result.outliers.extend(filtered_indices);
            };
        }
    }

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

    fn generate_multi_kgraph(&self, n_graphs: usize, keep_n_edges: usize, indices_to_include: &HashSet<usize>) -> Result<(KGraph<f64>, Vec<usize>)> {
        
        let graphs = (0..n_graphs).into_par_iter()
            .map(|_| {
                self.build_mutual_kgraph(keep_n_edges, &indices_to_include)
            })
            .collect::<Result<Vec<(_, _)>>>()?;

        // intersect the graphs to create a single graph
        let mut return_graph = None;
        for (graph, _) in graphs.into_iter() {
            if let Some(main_graph) = return_graph {
                return_graph = Some(intersect(main_graph, graph, keep_n_edges)?);
            } else {
                return_graph = Some(graph);
            }
        }

        let return_graph = return_graph.unwrap();
        let disconnected_nodes = return_graph
            .get_neighbours()
            .par_iter()
            .enumerate()
            .filter_map(|(node_index, nbrs)| {
                if nbrs.len() == 0 {
                    debug!("Node {} has no neighbours.", node_index);
                    Some(node_index)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        
        // self.filter_contigs(&return_graph, disconnected_nodes)?;

        Ok((return_graph, disconnected_nodes))
        
    }

    /// use Hnsw and annembed to iteratively build mutual K-nn graphs
    /// until we reach a point where all nodes have at leat 1 neighbour
    fn prefilter_nodes(&mut self, rounds: usize) -> Result<KGraph<f64>> {
        let mut indices_to_include = (0..self.n_contigs).collect::<HashSet<usize>>();
        // let (mut initial_kgraph, mut disconnected_nodes) = self.generate_multi_kgraph(3, 0, &indices_to_include)?;
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
            // (initial_kgraph, disconnected_nodes) = self.generate_multi_kgraph(3, keep_n_edges, &indices_to_include)?;
            (initial_kgraph, disconnected_nodes) = self.build_mutual_kgraph(keep_n_edges, &indices_to_include)?;
        }
        Ok(initial_kgraph)
    }

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

    fn generate_embedder_params(&self, b: f64) -> EmbedderParams {
        let mut params = EmbedderParams::default();
        let n_components = min(max(self.coverage_table.table.ncols() / 2, 2), 10);
        params.set_dim(n_components);
        // params.nb_grad_batch = self.nb_grad_batches;
        params.scale_rho = 1.0;
        params.beta = 0.5;
        params.b = 1.0;
        // params.b = b;

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
        cluster_results
    }

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
        while let Some(record) = reader.next() {
            let seqrec = record?;
            
            let contig_name = seqrec.id();
            let contig_name_str = std::str::from_utf8(contig_name)?;
            let contig_length = seqrec.seq().len();
            if contig_length < self.min_contig_size || self.filtered_contigs.contains(contig_name_str) {
                let cluster_label = if seqrec.seq().len() < self.min_bin_size {
                    UNBINNED.to_string()
                } else {
                    single_contig_bin_id += 1;
                    format!("single_contig_{}", single_contig_bin_id)
                };

                let bin_path = path::Path::new(&self.output_directory).join(format!("rosella_bin_{}{}", cluster_label, RECOVER_FASTA_EXTENSION));
                let file = OpenOptions::new().append(true).create(true).open(bin_path)?;
                let mut writer = BufWriter::new(file);
                // write contig to bin
                write_fasta(contig_name, &seqrec.seq(), &mut writer, LineEnding::Unix)?;
                continue;
            }
            // normalise sequence
            // let seqrec = seqrec.normalize(false);

            let cluster_result = &cluster_results[contig_idx];
            
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
        }

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