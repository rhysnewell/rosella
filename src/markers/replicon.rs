use log::debug;
use ndarray::{Array1, ArrayView2};

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
const GENE_QUANTILE: f64 = 0.25;
const LEAST_ANCHORS: usize = 30;
const LEAST_CARRIERS: usize = 3;

/// A publishing policy, not a measured bar. Below it the gene shape is mostly random open
/// reading frames in eukaryotic sequence, so the bin would be noise rather than a replicon.
pub const LEAST_BASES: usize = 10_000;

fn quantile(sorted: &[f64], at: f64) -> f64 {
    let last = sorted.len().saturating_sub(1);
    sorted[(at * last as f64).round() as usize]
}

#[derive(Clone, Copy, Debug)]
pub struct Bars {
    density: f64,
    gene: f64,
}

/// Both bars come off this assembly's own marker carrying contigs, not a tuned number.
pub fn bars(shapes: &[Shape], carries: &[bool], lengths: &[usize]) -> Option<Bars> {
    let anchors = (0..shapes.len())
        .filter(|at| carries[*at] && shapes[*at].genes >= LEAST_GENES)
        .collect::<Vec<_>>();
    if anchors.len() < LEAST_ANCHORS {
        debug!(
            "{} marker carrying contigs is too few to place the replicon bars",
            anchors.len()
        );
        return None;
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
    let bars = Bars {
        density: quantile(&densities, DENSITY_QUANTILE),
        gene: quantile(&genes, GENE_QUANTILE),
    };
    debug!(
        "Replicon bars over {} anchors: coding density {:.3}, mean gene {:.0} bp",
        anchors.len(),
        bars.density,
        bars.gene
    );
    Some(bars)
}

fn shaped(shape: &Shape, length: usize, bars: Bars) -> bool {
    length >= LEAST_BASES
        && shape.genes >= LEAST_GENES
        && shape.density(length) >= bars.density
        && shape.mean_gene() < bars.gene
}

/// Carriers are the bin's own, so only a contig further out than all of them is foreign. With
/// too few carriers to measure a spread the bin is mostly elements, and gene shape decides.
fn foreign(
    members: &[usize],
    carries: &[bool],
    lengths: &[usize],
    composition: ArrayView2<f64>,
) -> Vec<usize> {
    let members = members
        .iter()
        .copied()
        .filter(|contig| lengths[*contig] >= LEAST_BASES)
        .collect::<Vec<_>>();
    if members.iter().filter(|contig| carries[**contig]).count() < LEAST_CARRIERS {
        return members
            .into_iter()
            .filter(|contig| !carries[*contig])
            .collect();
    }
    let mut total = Array1::<f64>::zeros(composition.ncols());
    let mut weight = 0.0;
    for contig in &members {
        let length = lengths[*contig] as f64;
        total.scaled_add(length, &composition.row(*contig));
        weight += length;
    }
    let distances = members
        .iter()
        .map(|contig| {
            let length = lengths[*contig] as f64;
            let row = composition.row(*contig);
            let centre = (&total - &(&row * length)) / (weight - length);
            (&row - &centre).mapv(|value| value * value).sum().sqrt()
        })
        .collect::<Vec<_>>();
    let spread = members
        .iter()
        .zip(&distances)
        .filter(|(contig, _)| carries[**contig])
        .map(|(_, distance)| *distance)
        .fold(f64::MIN, f64::max);
    members
        .iter()
        .zip(&distances)
        .filter(|(contig, distance)| !carries[**contig] && **distance > spread)
        .map(|(contig, _)| *contig)
        .collect()
}

pub fn small_replicons<'a>(
    shapes: &[Shape],
    carries: &[bool],
    lengths: &[usize],
    bins: impl IntoIterator<Item = &'a [usize]>,
    composition: ArrayView2<f64>,
) -> Vec<usize> {
    let Some(bars) = bars(shapes, carries, lengths) else {
        return Vec::new();
    };
    let mut found = bins
        .into_iter()
        .flat_map(|members| foreign(members, carries, lengths, composition))
        .filter(|contig| shaped(&shapes[*contig], lengths[*contig], bars))
        .collect::<Vec<_>>();
    found.sort_unstable();
    debug!("{} small replicons", found.len());
    found
}
