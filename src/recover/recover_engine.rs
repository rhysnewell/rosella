use std::process::Command;

use annembed::{prelude::*, fromhnsw::{kgraph_from_hnsw_all, kgraph::KGraph}};
use anyhow::Result;
use finch::serialization::Sketch;
use hnsw_rs::prelude::Hnsw;
use log::{info, debug};
use ndarray::{Dimension, Dim};

use crate::{coverage::{coverage_calculator::{calculate_coverage, MetabatDistance}, coverage_table::CoverageTable}, kmers::kmer_counting::{count_kmers, KmerFrequencyTable}, sketch::{contig_sketcher::{ContigSketchResult, sketch_contigs}, sketch_distances::SketchDistance}, embedding::embedder::{ContigDistance, ContigInformation}};


pub fn run_recover(m: &clap::ArgMatches) -> Result<()> {
    let mut recover_engine = RecoverEngine::new(m)?;
    recover_engine.run()?;
    Ok(())
}

struct RecoverEngine {
    output_directory: String,
    coverage_table: CoverageTable,
    contig_sketches: ContigSketchResult,
    n_neighbours: usize,
    max_nb_connections: usize,
    nb_layers: usize,
    ef_construction: usize,
    max_layers: usize,
    n_contigs: usize,
    nb_grad_batches: usize,
    // contig_nn: Hnsw<ContigInformation<'static, D>, ContigDistance>,
}

impl RecoverEngine {
    pub fn new(m: &clap::ArgMatches) -> Result<Self> {
        // create output directory
        let output_directory = m.get_one::<String>("output-directory").unwrap().clone();
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

        Ok(
            Self {
                output_directory,
                coverage_table,
                contig_sketches,
                n_neighbours,
                max_nb_connections,
                nb_layers,
                ef_construction,
                max_layers,
                n_contigs,
                nb_grad_batches,
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
        let contig_data_for_insertion = contig_data
            .iter()
            .enumerate()
            .map(|(i, contig_information)| {
                (contig_information, i)
            })
            .collect::<Vec<_>>();
        contig_nn.parallel_insert(&contig_data_for_insertion);
        
        info!("Constructing kgraph.");
        let kgraph: KGraph<f32> = kgraph_from_hnsw_all(&contig_nn, self.n_neighbours).unwrap();
        let embedder_params = self.generate_embedder_params();
        let mut embedder = Embedder::new(&kgraph, embedder_params);
        info!("Embedding.");
        let embed_result = embedder.embed().unwrap();
        let embeddings = embedder.get_embedded_reindexed();
        debug!("Embeddings: {:?}", embeddings);
        Ok(())
    }

    fn generate_embedder_params(&self) -> EmbedderParams {
        let mut params = EmbedderParams::default();
        params.set_dim(2);

        params
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

    fn run_flight(&self) -> Result<()> {
        let mut flight_command = Command::new("flight");
        // flight_command.arg("bin")
        Ok(())
    }
}

