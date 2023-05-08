use std::{process::Command, collections::HashMap, cmp::Ordering, path, io::BufWriter, fs::{File, OpenOptions}};

use annembed::{prelude::*, fromhnsw::{kgraph_from_hnsw_all, kgraph::KGraph}};
use anyhow::Result;
use finch::serialization::Sketch;
use hnsw_rs::prelude::Hnsw;
use log::{info, debug};
use ndarray::{Dimension, Dim, ArrayBase, OwnedRepr};
use needletail::{parse_fastx_file, Sequence, parser::{write_fasta, LineEnding}};
use petal_clustering::{HDbscan, Fit, Predict};
use rayon::slice::ParallelSliceMut;

use crate::{coverage::{coverage_calculator::{calculate_coverage, MetabatDistance}, coverage_table::CoverageTable}, kmers::kmer_counting::{count_kmers, KmerFrequencyTable}, sketch::{contig_sketcher::{ContigSketchResult, sketch_contigs}, sketch_distances::SketchDistance}, embedding::embedder::{ContigDistance, ContigInformation}};


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
                if file_name.ends_with(".fna") {
                    panic!("Output directory contains .fna files. Please remove them before running rosella recover.");
                }
            }
        }

        let assembly = m.get_one::<String>("assembly").unwrap().clone();
        // create the output directory but do not fail if it already exists
        std::fs::create_dir_all(&output_directory)?;
        info!("Calculating contig coverages.");
        let coverage_table = calculate_coverage(m)?;
        info!("Calculating contig sketches.");
        let contig_sketches = sketch_contigs(m)?;

        let n_neighbours = m.get_one::<usize>("n-neighbours").unwrap().clone();
        let max_nb_connections = m.get_one::<usize>("max-nb-connections").unwrap().clone();
        let nb_layers = m.get_one::<usize>("nb-layers").unwrap().clone();
        let ef_construction = m.get_one::<usize>("ef-construction").unwrap().clone();
        let max_layers = m.get_one::<usize>("max-layers").unwrap().clone();
        let nb_grad_batches = m.get_one::<usize>("nb-grad-batches").unwrap().clone();
        let n_contigs = contig_sketches.contig_names.len();
        let min_bin_size = m.get_one::<usize>("min-bin-size").unwrap().clone();

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
            }
        )
    }

    /// Runs through the rosella bin recovery pipeline
    pub fn run(&mut self) -> Result<()> {
        
        let contig_nn: Hnsw<ContigInformation<'_, Dim<[usize; 1]>>, ContigDistance> = Hnsw::new(
            self.n_neighbours, 
            self.n_contigs, 
            self.max_layers, 
            self.ef_construction,
            ContigDistance
        );      

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
        debug!("Embeddings: {:?}", embeddings);
        info!("Clustering.");
        let cluster_labels = self.cluster(embeddings);
        debug!("Cluster labels: {:?}", cluster_labels);
        
        self.write_clusters(cluster_labels)?;

        Ok(())
    }

    fn generate_embedder_params(&self) -> EmbedderParams {
        let mut params = EmbedderParams::default();
        params.set_dim(2);
        params.nb_grad_batch = self.nb_grad_batches;
        params.scale_rho = 1.5;
        params.beta = 1.;
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

        let mut cluster_results = self.get_cluster_result(cluster_map, outlier);
        cluster_results.par_sort_unstable();
        cluster_results
    }

    fn get_cluster_result(&self, cluster_map: HashMap<usize, Vec<usize>>, outliers: Vec<usize>) -> Vec<ClusterResult> {
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
        }
        for outlier in outliers {
            cluster_results.push(ClusterResult::new(outlier, None));
        }
        cluster_results
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
                        "unbinned".to_string()
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