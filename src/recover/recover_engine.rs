use std::{process::Command, collections::HashMap, cmp::Ordering, path, io::{BufWriter, Write}, fs::{File, OpenOptions}};

use annembed::{prelude::*, fromhnsw::{kgraph_from_hnsw_all, kgraph::KGraph}};
use anyhow::Result;
use hnsw_rs::prelude::{Hnsw, Distance};
use log::{info, debug, warn, error};
use ndarray::{Dimension, Dim, ArrayBase, OwnedRepr, Array2};
use needletail::{parse_fastx_file, Sequence, parser::{write_fasta, LineEnding}};
use petal_clustering::{HDbscan, Fit, Predict};
use rayon::{prelude::*, slice::ParallelSliceMut, prelude::IntoParallelIterator};


use crate::{coverage::{coverage_calculator::{calculate_coverage, MetabatDistance}, coverage_table::CoverageTable, self}, kmers::kmer_counting::{count_kmers, KmerFrequencyTable}, sketch::{contig_sketcher::{ContigSketchResult, sketch_contigs}, sketch_distances::{SketchDistance, distance}}, embedding::embedder::{ContigDistance, ContigInformation}};

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
    nb_layers: usize,
    ef_construction: usize,
    max_layers: usize,
    n_contigs: usize,
    nb_grad_batches: usize,
    min_bin_size: usize,
    min_contig_size: usize,
    embeddings: Option<Array2<f64>>,
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
                    error!("Output directory contains .fna files. Please remove them before running rosella recover.");
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
        coverage_table.filter(min_contig_size)?;
        info!("Calculating contig sketches.");
        let contig_sketches = sketch_contigs(m)?;
        assert_eq!(coverage_table.table.nrows(), contig_sketches.contig_sketches.len(), "Coverage table and contig sketches have different number of contigs.");

        let n_neighbours = m.get_one::<usize>("n-neighbours").unwrap().clone();
        let max_nb_connections = m.get_one::<usize>("max-nb-connections").unwrap().clone();
        let nb_layers = m.get_one::<usize>("nb-layers").unwrap().clone();
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
                nb_layers,
                ef_construction,
                max_layers,
                n_contigs,
                nb_grad_batches,
                min_bin_size,
                min_contig_size,
                embeddings: None,
            }
        )
    }

    /// Runs through the rosella bin recovery pipeline
    pub fn run(&mut self) -> Result<()> {
        
        let mut contig_nn: Hnsw<ContigInformation<'_, Dim<[usize; 1]>>, ContigDistance> = Hnsw::new(
            self.n_neighbours, 
            self.n_contigs, 
            self.max_layers, 
            self.ef_construction,
            ContigDistance
        );

        if self.n_contigs < SMALL_DATASET_THRESHOLD {
            warn!("Small dataset detected. Keeping pruned neighbourhoods.");
            contig_nn.set_keeping_pruned(true);
        }

        let contig_data = self.retrieve_contig_data();
        info!("Inserting contig data into HNSW.");
        let contig_data_for_insertion = contig_data
            .iter()
            .enumerate()
            .map(|(i, contig_information)| {
                (contig_information, i)
            })
            .collect::<Vec<_>>();
        contig_nn.parallel_insert(&contig_data_for_insertion);
        
        info!("Constructing kgraph.");
        let kgraph: KGraph<f64> = kgraph_from_hnsw_all(&contig_nn, self.n_neighbours).unwrap();
        let embedder_params = self.generate_embedder_params();
        let mut embedder = Embedder::new(&kgraph, embedder_params);
        info!("Embedding.");
        let _ = embedder.embed().unwrap();
        let embeddings = embedder.get_embedded_reindexed();
        // self.embeddings = Some(embeddings.clone());
        debug!("Embeddings: {:?}", embeddings);
        info!("Clustering.");
        let cluster_labels = self.cluster(embeddings);
        debug!("Cluster labels: {:?}", cluster_labels);
        
        info!("Writing clusters.");
        self.write_clusters(cluster_labels)?;

        Ok(())
    }

    fn generate_embedder_params(&self) -> EmbedderParams {
        let mut params = EmbedderParams::default();
        params.set_dim(2);
        params.nb_grad_batch = self.nb_grad_batches;
        params.scale_rho = 1.0;
        params.beta = 1.0;
        params.b = 0.5;
        params.grad_step = 1.;
        params.nb_sampling_by_edge = 10;
        params.dmap_init = true;

        params
    }

    /// Clusters the embeddings using HDBSCAN
    fn cluster(&self, embeddings: ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>) -> Vec<ClusterResult> {
        let mut clusterer = HDbscan::default();
        clusterer.min_samples = 10;
        clusterer.min_cluster_size = 10;
        clusterer.alpha = 1.0;
        let (cluster_map, outlier) = clusterer.fit(&embeddings);

        let mut cluster_results = self.get_cluster_result(cluster_map, outlier, &embeddings);
        cluster_results.par_sort_unstable();
        cluster_results
    }

    fn get_cluster_result(&self, cluster_map: HashMap<usize, Vec<usize>>, outliers: Vec<usize>, embeddings: &ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>) -> Vec<ClusterResult> {
        if DEBUG_BINS {
            // delete bins.txt
            let _ = std::fs::remove_file(format!("{}/bins.txt", self.output_directory));
        }

        let mut cluster_results = Vec::with_capacity(self.n_contigs);
        for (cluster_label, contig_indices) in cluster_map.iter() {
            // check the size of the cluster and if it is too small, set the cluster label to None
            let bin_size = contig_indices.iter().map(|i| self.coverage_table.contig_lengths[*i]).sum::<usize>();
            let cluster_label = if bin_size < self.min_bin_size {
                None
            } else {
                Some(*cluster_label)
            };
            for contig_index in contig_indices {
                cluster_results.push(ClusterResult::new(*contig_index, cluster_label));
            }

            if DEBUG_BINS {
                let file = OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(format!("{}/bins.txt", self.output_directory))
                    .unwrap();
                let mut writer = BufWriter::new(file);
                let (size, metabat, sketch, embedding) = self.calculate_bin_metrics(contig_indices, embeddings);
                writeln!(writer, "{:?}\t{}\t{}\t{}\t{}\t{}", cluster_label, contig_indices.len(), size, metabat, sketch, embedding).unwrap();
            }
        }
        for outlier in outliers {
            cluster_results.push(ClusterResult::new(outlier, None));
        }
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

    fn retrieve_contig_data<'a>(&'a self) -> Vec<Vec<ContigInformation<'a, Dim<[usize; 1]>>>> {
        let contig_data = self.coverage_table.table.rows().into_iter()
            .zip(self.contig_sketches.contig_sketches.iter())
            .map(|(row, sketch)| {
                let contig_information = ContigInformation::new(
                    row,
                    sketch,
                );
                vec![contig_information]
            }).collect::<Vec<_>>();
        
        return contig_data;
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
            let contig_length = seqrec.seq().len();
            if contig_length < self.min_contig_size {
                let bin_path = path::Path::new(&self.output_directory).join(format!("rosella_bin_{}{}", UNBINNED, RECOVER_FASTA_EXTENSION));
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