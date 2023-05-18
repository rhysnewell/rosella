use annembed::fromhnsw::{kgraph::KGraph, kgraph_from_hnsw_all};
use anyhow::Result;
use finch::serialization::Sketch;
use hnsw_rs::prelude::{Distance, Hnsw};
use log::{debug, info};
use ndarray::{ArrayView, Dimension, Dim};
use rayon::prelude::*;

use crate::{coverage::coverage_calculator::MetabatDistance, sketch::sketch_distances::SketchDistance, graphs::nearest_neighbour_graph::{mutual_knn, intersect}, kmers::kmer_counting::KmerCorrelation};

pub struct EmbedderEngine<'a> {
    pub contig_information: Vec<Vec<ContigInformation<'a, Dim<[usize; 1]>>>>,
    pub distance: DepthDistance,
}

impl<'a> EmbedderEngine<'a> {
    pub fn new(contig_information: Vec<Vec<ContigInformation<'a, Dim<[usize; 1]>>>>) -> Self {
        Self {
            contig_information,
            distance: DepthDistance,
        }
    }

    /// Build a mutual K-NN graph, but keep at least `keep_n_edges` edges per node
    /// Setting `keep_n_edges` to 0 can result in some nodes becoming disconnected, these nodes
    /// should either be removed or reconnected to the graph
    pub fn build_mutual_kgraph(&self, keep_n_edges: usize, n_neighbours: usize, n_contigs: usize, max_layers: usize, ef_construction: usize) -> Result<(KGraph<f64>, Vec<usize>)> {
        let mut depth_nn: Hnsw<ContigInformation<'_, Dim<[usize; 1]>>, DepthDistance> = Hnsw::new(
            n_neighbours, 
            n_contigs, 
            max_layers, 
            ef_construction,
            DepthDistance
        );
        depth_nn.set_keeping_pruned(true);

        let mut tnf_nn: Hnsw<ContigInformation<'_, Dim<[usize; 1]>>, TNFDistance> = Hnsw::new(
            n_neighbours, 
            n_contigs, 
            max_layers, 
            ef_construction,
            TNFDistance
        );
        tnf_nn.set_keeping_pruned(true);

        info!("Inserting contig data into HNSW.");
        let contig_data_for_insertion = self.contig_information
            .iter()
            .enumerate()
            .map(|(i, contig_information)| {
                (contig_information, i)
            })
            .collect::<Vec<_>>();

        // Insert data into both NN graphs
        depth_nn.parallel_insert(&contig_data_for_insertion);
        tnf_nn.parallel_insert(&contig_data_for_insertion);
        
        info!("Constructing kgraph.");
        let depth_kgraph: KGraph<f64> = kgraph_from_hnsw_all(&depth_nn, n_neighbours).unwrap();
        let tnf_kgraph: KGraph<f64> = kgraph_from_hnsw_all(&tnf_nn, n_neighbours).unwrap();
        debug!("Intersection of depth and tnf kgraphs.");
        // Intersect both graphs, ensuring keep_n_edges are kept in the final graph
        let mut kgraph = intersect(depth_kgraph, tnf_kgraph, keep_n_edges + 1)?;
        kgraph = mutual_knn(kgraph, keep_n_edges)?;

        let disconnected_nodes = kgraph
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
        Ok((kgraph, disconnected_nodes))
    }
}

#[derive(Debug, Clone)]
pub struct ContigInformation<'a, D: Dimension> {
    pub coverage: ArrayView<'a, f64, D>,
    pub tnf: ArrayView<'a, f64, D>,
    pub sketch: &'a Sketch
}

impl<'a, D: Dimension> ContigInformation<'a, D> {
    pub fn new(coverage: ArrayView<'a, f64, D>, tnf: ArrayView<'a, f64, D>, sketch: &'a Sketch) -> Self {
        Self {
            coverage,
            tnf,
            sketch
        }
    }
}

#[derive(Debug, Clone)]
pub struct TNFDistance;
impl<'a, D: Dimension> Distance<ContigInformation<'a, D>> for TNFDistance {
    fn eval(&self, va: &[ContigInformation<D>], vb: &[ContigInformation<D>]) -> f32 {
        let tnf_correlation = va.iter()
            .map(|query_contig| {
                vb.iter()
                    .map(|ref_contig| {
                        KmerCorrelation::distance(query_contig.tnf.view(), ref_contig.tnf.view())
                    })
                    .sum::<f64>()
                    / vb.len() as f64
            })
            .sum::<f64>()
            / va.len() as f64;
        
        if tnf_correlation.is_nan() {
            return 0.0
        }
        tnf_correlation as f32
    }
}


#[derive(Debug, Clone)]
pub struct DepthDistance;

impl<'a, D: Dimension> Distance<ContigInformation<'a, D>> for DepthDistance {
    fn eval(&self, va: &[ContigInformation<D>], vb: &[ContigInformation<D>]) -> f32 {
        // mean metabat distance between all pairs of contigs
        let metabat_distance = va.iter()
            .map(|query_contig| {
                vb.iter()
                    .map(|ref_contig| {
                        MetabatDistance::distance(query_contig.coverage.view(), ref_contig.coverage.view())
                    })
                    .sum::<f64>()
                    / vb.len() as f64
            })
            .sum::<f64>()
            / va.len() as f64;
        
        // mean min_jaccard distance between all pairs of contigs
        let min_jaccard_distance = va.iter()
            .map(|query_contig| {
                vb.iter()
                    .map(|ref_contig| {
                        SketchDistance.eval(&[query_contig.sketch], &[ref_contig.sketch])
                    })
                    .sum::<f32>()
                    / vb.len() as f32
            })
            .sum::<f32>()
            / va.len() as f32;
        
        let result = metabat_distance * min_jaccard_distance as f64;
        if result.is_nan() {
            return 0.0
        }

        result as f32
    }
}