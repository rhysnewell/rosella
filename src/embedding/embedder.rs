use finch::serialization::Sketch;
use hnsw_rs::prelude::Distance;
use log::debug;
use ndarray::{ArrayView, Dimension, Dim};

use crate::{coverage::coverage_calculator::MetabatDistance, sketch::sketch_distances::SketchDistance};

#[derive(Debug, Clone)]
pub struct ContigInformation<'a, D: Dimension> {
    pub coverage: ArrayView<'a, f64, D>,
    pub sketch: &'a Sketch
}

impl<'a, D: Dimension> ContigInformation<'a, D> {
    pub fn new(coverage: ArrayView<'a, f64, D>, sketch: &'a Sketch) -> Self {
        Self {
            coverage,
            sketch
        }
    }
}

#[derive(Debug, Clone)]
pub struct ContigDistance;

impl<'a, D: Dimension> Distance<ContigInformation<'a, D>> for ContigDistance {
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
        
        result as f32
    }
}