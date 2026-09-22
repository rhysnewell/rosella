use std::collections::HashSet;

use log::debug;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Shape {
    pub coding_bases: u64,
    pub genes: u32,
}

impl Shape {
    pub fn add(&mut self, bases: usize) {
        self.coding_bases += bases as u64;
        self.genes += 1;
    }

    fn density(&self, length: usize) -> f64 {
        if length == 0 {
            return 0.0;
        }
        self.coding_bases as f64 / length as f64
    }

    fn mean_gene(&self) -> f64 {
        if self.genes == 0 {
            return 0.0;
        }
        self.coding_bases as f64 / f64::from(self.genes)
    }
}

/// Fewer genes than this and the mean is one gene's length, not the contig's gene shape.
const LEAST_GENES: u32 = 10;
const DENSITY_QUANTILE: f64 = 0.50;
const GENE_QUANTILE: f64 = 0.05;
const LEAST_ANCHORS: usize = 30;

fn quantile(sorted: &[f64], at: f64) -> f64 {
    let last = sorted.len().saturating_sub(1);
    sorted[(at * last as f64).round() as usize]
}

/// A small replicon is as gene dense as a chromosome, carries genes shorter than nearly any
/// chromosome's, and holds no single copy marker. Both bars are read off this assembly's own
/// marker carrying contigs, so nothing here is a tuned constant.
pub fn small_elements(shapes: &[Shape], carries: &[bool], lengths: &[usize]) -> HashSet<usize> {
    let anchors = (0..shapes.len())
        .filter(|at| carries[*at] && shapes[*at].genes >= LEAST_GENES)
        .collect::<Vec<_>>();
    if anchors.len() < LEAST_ANCHORS {
        debug!(
            "{} marker carrying contigs is too few to place the small element bars",
            anchors.len()
        );
        return HashSet::new();
    }

    let mut densities = anchors
        .iter()
        .map(|at| shapes[*at].density(lengths[*at]))
        .collect::<Vec<_>>();
    let mut genes = anchors
        .iter()
        .map(|at| shapes[*at].mean_gene())
        .collect::<Vec<_>>();
    densities.sort_unstable_by(f64::total_cmp);
    genes.sort_unstable_by(f64::total_cmp);
    let dense_bar = quantile(&densities, DENSITY_QUANTILE);
    let gene_bar = quantile(&genes, GENE_QUANTILE);

    let found = (0..shapes.len())
        .filter(|at| {
            !carries[*at]
                && shapes[*at].genes >= LEAST_GENES
                && shapes[*at].density(lengths[*at]) >= dense_bar
                && shapes[*at].mean_gene() < gene_bar
        })
        .collect::<HashSet<_>>();
    debug!(
        "{} small elements over {} anchors, coding density {dense_bar:.3}, mean gene {gene_bar:.0} bp",
        found.len(),
        anchors.len()
    );
    found
}
