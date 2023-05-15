use std::{collections::{HashMap, HashSet}, cmp::Ordering, path, io::{BufWriter, Write}, fs::{OpenOptions}};

use annembed::{prelude::*, fromhnsw::{kgraph_from_hnsw_all, kgraph::KGraph}};
use anyhow::Result;
use log::{info, debug, warn, error};
use ndarray::{Dim, ArrayBase, OwnedRepr, Array2};
use needletail::{parse_fastx_file, Sequence, parser::{write_fasta, LineEnding}};
use petal_clustering::{HDbscan, Fit};
use rayon::{prelude::*, slice::ParallelSliceMut, prelude::IntoParallelIterator};


use crate::{coverage::{coverage_calculator::{calculate_coverage, MetabatDistance}, coverage_table::CoverageTable, self}, kmers::kmer_counting::{count_kmers, KmerFrequencyTable}, sketch::{contig_sketcher::{ContigSketchResult, sketch_contigs}, sketch_distances::{SketchDistance, distance}}, embedding::embedder::{ContigDistance, ContigInformation, EmbedderEngine}, graphs::nearest_neighbour_graph::mutual_knn, clustering::cluster_utils::{condensed_pairwise_distance, silhouette_score}};

const RECOVER_FASTA_EXTENSION: &str = ".fna";
const DEBUG_BINS: bool = true;
const UNBINNED: &str = "unbinned";
const SMALL_DATASET_THRESHOLD: usize = 10000;

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
    n_neighbours: usize,
    max_nb_connections: usize,
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
        let filtered_contigs = coverage_table.filter_by_length(min_contig_size)?;
        info!("Calculating contig sketches.");
        let contig_sketches = sketch_contigs(m)?;
        assert_eq!(coverage_table.table.nrows(), contig_sketches.contig_sketches.len(), "Coverage table and contig sketches have different number of contigs.");

        let n_neighbours = m.get_one::<usize>("n-neighbours").unwrap().clone();
        let max_nb_connections = m.get_one::<usize>("max-nb-connections").unwrap().clone();
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
                n_neighbours,
                max_nb_connections,
                filtering_rounds,
                ef_construction,
                max_layers,
                n_contigs,
                nb_grad_batches,
                min_bin_size,
                min_contig_size,
                embeddings: None,
                filtered_contigs,
            }
        )
    }

    /// Runs through the rosella bin recovery pipeline
    pub fn run(&mut self) -> Result<()> {
        
        // 1. Perform N rounds of filtering. Filtering involves removinf disconnnected contigs from the graph.
        //    The final round ensures each contig has at least one neighbour.
        let kgraph = self.prefilter_nodes(self.filtering_rounds)?;
        
        // 2. Take the KGraph and embed using a UMAP-esque method.
        //    Embedding parameters need to be chosen more carefully.
        info!("Embedding.");
        let embeddings = self.embed(&kgraph)?;
        debug!("Embeddings: {:?}", embeddings);

        // 3. Cluster the embeddings using HDBSCAN.
        //    As with the embedding parameters, we need to better choose the clustering parameters.
        info!("Clustering.");
        let (mut cluster_map, mut outlier) = self.cluster(&embeddings, 2);


        // 4. After the initial embedding and clustering, this is where things get tricky. We need to filter out unclustered contigs
        //    or clusters that we think might be too small and then re-embed and re-cluster.


        // 5. Then we need to inspect each cluster individually and see if we can split it into smaller clusters.
        //    This involves embedding the contigs in each cluster and then clustering them again. If the resulting cluster looks
        //    decent we can keep it.

        let cluster_results = self.get_cluster_result(cluster_map, outlier, &embeddings, None);
        
        info!("Writing clusters.");
        self.write_clusters(cluster_results)?;

        Ok(())
    }

    fn embed(&self, kgraph: &KGraph<f64>) -> Result<Array2<f64>> {
        let embedder_params = self.generate_embedder_params();
        let mut embedder = Embedder::new(&kgraph, embedder_params);
        let _ = embedder.embed().unwrap();
        let embeddings = embedder.get_embedded_reindexed();
        Ok(embeddings)
    }

    /// use Hnsw and annembed to iteratively build mutual K-nn graphs
    /// until we reach a point where all nodes have at leat 1 neighbour
    fn prefilter_nodes(&mut self, rounds: usize) -> Result<KGraph<f64>> {
        let mut indices_to_include = (0..self.n_contigs).collect::<HashSet<usize>>();
        let (mut initial_kgraph, mut disconnected_contigs) = self.build_mutual_kgraph(0, &indices_to_include)?;
        let mut current_round = 1;
        let mut keep_n_edges = 0;
        while !disconnected_contigs.is_empty() {
            info!("Round {} - Filtering {} contigs.", current_round, disconnected_contigs.len());
            self.filter_contigs(&initial_kgraph, disconnected_contigs)?;
            indices_to_include = (0..self.n_contigs).collect::<HashSet<usize>>();
            current_round += 1;
            if current_round > rounds {
                // we've hit the max number of rounds, so esnure the next graph has
                // no disconnected nodes
                keep_n_edges += 1;
            }
            (initial_kgraph, disconnected_contigs) = self.build_mutual_kgraph(keep_n_edges, &indices_to_include)?;
        }
        Ok(initial_kgraph)
    }

    fn filter_contigs(&mut self, kgraph: &KGraph<f64>, disconnected_contigs: Vec<usize>) -> Result<()> {
        // use the kgraph to take the disconnected_contig indices and get the original indices
        // then filter the coverage table and contig sketches by index
        let original_indices = disconnected_contigs
            .iter()
            .filter_map(|idx| {
                kgraph.get_data_id_from_idx(*idx)
            })
            .map(|idx| {
                *idx
            })
            .collect::<HashSet<_>>();
        self.filtered_contigs.extend(self.coverage_table.filter_by_index(&original_indices)?);
        self.contig_sketches.filter_by_index(&original_indices)?;
        self.n_contigs = self.contig_sketches.contig_names.len();

        Ok(())
    }

    /// Build a mutual K-NN graph, but keep at least `keep_n_edges` edges per node
    /// Setting `keep_n_edges` to 0 can result in some nodes becoming disconnected, these nodes
    /// should either be removed or reconnected to the graph
    fn build_mutual_kgraph(&self, keep_n_edges: usize, indices_to_include: &HashSet<usize>) -> Result<(KGraph<f64>, Vec<usize>)> {
        let contig_data = self.retrieve_contig_data(indices_to_include);
        let embedder_engine = EmbedderEngine::new(contig_data);
        embedder_engine.build_mutual_kgraph(keep_n_edges, self.n_neighbours, self.n_contigs, self.max_layers, self.ef_construction)
    }

    fn retrieve_contig_data<'a>(&'a self, indices_to_include: &HashSet<usize>) -> Vec<Vec<ContigInformation<'a, Dim<[usize; 1]>>>> {
        let contig_data = self.coverage_table.table
            .rows()
            .into_iter()
            .zip(self.contig_sketches.contig_sketches.iter())
            .enumerate()
            .filter_map(|(row_index, (row, sketch))| {
                if !indices_to_include.contains(&row_index) {
                    return None;
                }

                let contig_information = ContigInformation::new(
                    row,
                    sketch,
                );
                Some(vec![contig_information])
            }).collect::<Vec<_>>();
        
        return contig_data;
    }

    fn generate_embedder_params(&self) -> EmbedderParams {
        let mut params = EmbedderParams::default();
        params.set_dim(2);
        params.nb_grad_batch = self.nb_grad_batches;
        params.scale_rho = 1.0;
        params.beta = 1.0;
        params.b = 0.4;
        params.grad_step = 1.;
        params.nb_sampling_by_edge = 10;
        params.dmap_init = true;

        params
    }

    /// Clusters the embeddings using HDBSCAN
    fn cluster(&self, embeddings: &ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>, starting_min_cluster_size: usize) -> (HashMap<usize, Vec<usize>>, Vec<usize>) {
        let n = embeddings.shape()[0]; // number of rows
        let end_size = starting_min_cluster_size + 10;

        let condensed_distances = condensed_pairwise_distance(embeddings);

        // for min_samples 5..15
        // calculate the best silhouette score, return the cluster map and outliers
        let mut cluster_results = (starting_min_cluster_size..end_size).into_par_iter()
            .flat_map(|min_cluster_size| {
                (2..15)
                    .into_par_iter()
                    .map(|min_samples| {
                        let mut clusterer = HDbscan::default();
                        clusterer.min_samples = min_samples;
                        clusterer.min_cluster_size = min_cluster_size;
                        clusterer.alpha = 1.0;
                        let (cluster_map, outliers) = clusterer.fit(&embeddings);
                        let mut s_score = silhouette_score(&condensed_distances, &cluster_map, n).expect("Failed to calculate silhouette score");
                        if s_score.is_nan() {
                            s_score = 0.0;
                        }
                        debug!("Min cluster size: {}, min samples: {}, silhouette score: {}", min_cluster_size, min_samples, s_score);
                        (cluster_map, outliers, s_score)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        
        // sort by silhouette score, closer to 1 is better
        cluster_results.sort_by(|a, b| {
            a.2.partial_cmp(&b.2).unwrap_or(Ordering::Equal).reverse()
        });

        let (cluster_map, outliers, score) = cluster_results.remove(0);
        info!("Best silhouette score: {}", score);
        (cluster_map, outliers)
    }

    fn get_cluster_result(
        &self,
        cluster_map: HashMap<usize, Vec<usize>>,
        outliers: Vec<usize>,
        embeddings: &ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>,
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
                        .collect::<Vec<_>>()
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

            if DEBUG_BINS {
                let file = OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(format!("{}/bins.txt", self.output_directory))
                    .unwrap();
                let mut writer = BufWriter::new(file);
                let (size, metabat, sketch, embedding) = self.calculate_bin_metrics(contig_indices.as_slice(), embeddings);
                writeln!(writer, "{:?}\t{}\t{}\t{}\t{}\t{}", cluster_label, contig_indices.len(), size, metabat, sketch, embedding).unwrap();
            }
        }
        for outlier in outliers {
            cluster_results.push(ClusterResult::new(outlier, None));
        }
        cluster_results.par_sort_unstable();
        cluster_results
    }

    /// calculate average distance between contigs in a bin
    fn calculate_bin_metrics(&self, contig_indices: &[usize], embeddings: &ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>) -> (usize, f32, f32, f32) {
        let bin_size = contig_indices.iter().map(|i| self.coverage_table.contig_lengths[*i]).sum::<usize>();
        
        // calculate all pairwise metabat distances
        let n_combinations = (contig_indices.len() * (contig_indices.len() - 1)) / 2;
        let distances = (0..contig_indices.len())
            .into_par_iter()
            .flat_map(|i| {
                let va_cov = self.coverage_table.table.row(i);
                let va_sketch = &self.contig_sketches.contig_sketches[i];
                let va_embedding = embeddings.row(i);
                (i+1..contig_indices.len())
                    .into_par_iter()
                    .map(move |j| {
                        let vb_cov = self.coverage_table.table.row(j);
                        let vb_sketch = &self.contig_sketches.contig_sketches[j];
                        let vb_embedding = embeddings.row(j);
                        let mut metabat_dist = MetabatDistance::distance(va_cov, vb_cov);
                        if metabat_dist.is_nan() {
                            metabat_dist = 1.;
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

                        (metabat_dist as f32, sketch_dist as f32, embedding_dist as f32)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let metabat_dist = distances.iter().map(|(metabat_dist, _, _)| *metabat_dist).sum::<f32>() / n_combinations as f32;
        let sketch_dist = distances.iter().map(|(_, sketch_dist, _)| *sketch_dist).sum::<f32>() / n_combinations as f32;
        let embedding_dist = distances.iter().map(|(_, _, embedding_dist)| *embedding_dist).sum::<f32>() / n_combinations as f32;
        (bin_size, metabat_dist, sketch_dist, embedding_dist)
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
            // // normalise sequence
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