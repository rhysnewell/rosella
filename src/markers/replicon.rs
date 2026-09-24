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
/// Composition alone also ejects the bin's own islands, which is why depth must agree.
const PASSENGER_QUANTILE: f64 = 0.98;
/// Composition noise falls with length, so a contig is judged against carriers of its length.
const OCTAVES: usize = 5;
const DEPTH_FLOOR: f64 = 0.1;

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

#[derive(Debug, Default, PartialEq, Eq)]
pub struct Departures {
    pub replicons: Vec<usize>,
    pub passengers: Vec<usize>,
}

struct Measured {
    contig: usize,
    carrier: bool,
    composition: f64,
    depth: f64,
}

fn octave(length: usize) -> usize {
    ((length / LEAST_BASES).ilog2() as usize).min(OCTAVES - 1)
}

fn composition_distances(
    members: &[usize],
    lengths: &[usize],
    composition: ArrayView2<f64>,
) -> Vec<f64> {
    let mut total = Array1::<f64>::zeros(composition.ncols());
    let mut weight = 0.0;
    for contig in members {
        let length = lengths[*contig] as f64;
        total.scaled_add(length, &composition.row(*contig));
        weight += length;
    }
    members
        .iter()
        .map(|contig| {
            let length = lengths[*contig] as f64;
            let row = composition.row(*contig);
            let centre = (&total - &(&row * length)) / (weight - length);
            (&row - &centre).mapv(|value| value * value).sum().sqrt()
        })
        .collect()
}

/// Against the median of the others, so passengers cannot drag the reference toward themselves.
fn depth_offsets(members: &[usize], depths: ArrayView2<f64>) -> Vec<f64> {
    let rest = members.len() - 1;
    let mut squares = vec![0.0; members.len()];
    for sample in depths.columns() {
        let logs = members
            .iter()
            .map(|contig| (sample[*contig] + DEPTH_FLOOR).log2())
            .collect::<Vec<_>>();
        let mut order = (0..members.len()).collect::<Vec<_>>();
        order.sort_unstable_by(|left, right| logs[*left].total_cmp(&logs[*right]));
        for (rank, at) in order.iter().enumerate() {
            let kth = |k: usize| logs[order[if k < rank { k } else { k + 1 }]];
            let median = if rest % 2 == 1 {
                kth(rest / 2)
            } else {
                (kth(rest / 2 - 1) + kth(rest / 2)) / 2.0
            };
            squares[*at] += (logs[*at] - median).powi(2);
        }
    }
    squares
        .into_iter()
        .map(|sum| (sum / depths.ncols() as f64).sqrt())
        .collect()
}

/// An octave with too few carriers borrows the bar below it, which is the looser one.
fn limits(measured: &[Measured], lengths: &[usize], value: fn(&Measured) -> f64) -> [f64; OCTAVES] {
    let mut out = [f64::INFINITY; OCTAVES];
    let mut below = f64::INFINITY;
    for (band, limit) in out.iter_mut().enumerate() {
        let mut values = measured
            .iter()
            .filter(|measured| measured.carrier && octave(lengths[measured.contig]) == band)
            .map(value)
            .collect::<Vec<_>>();
        if values.len() >= LEAST_ANCHORS {
            values.sort_unstable_by(f64::total_cmp);
            below = quantile(&values, PASSENGER_QUANTILE);
        }
        *limit = below;
    }
    out
}

/// With too few carriers to measure a spread the bin is mostly elements, and gene shape decides.
pub fn departures<'a>(
    shapes: &[Shape],
    carries: &[bool],
    lengths: &[usize],
    bins: impl IntoIterator<Item = &'a [usize]>,
    composition: ArrayView2<f64>,
    depths: ArrayView2<f64>,
) -> Departures {
    let Some(bars) = bars(shapes, carries, lengths) else {
        return Departures::default();
    };
    let is_shaped = |contig: &usize| shaped(&shapes[*contig], lengths[*contig], bars);
    let mut replicons = Vec::new();
    let mut measured = Vec::new();
    for members in bins {
        let members = members
            .iter()
            .copied()
            .filter(|contig| lengths[*contig] >= LEAST_BASES)
            .collect::<Vec<_>>();
        if members.iter().filter(|contig| carries[**contig]).count() < LEAST_CARRIERS {
            replicons.extend(
                members
                    .into_iter()
                    .filter(|contig| !carries[*contig] && is_shaped(contig)),
            );
            continue;
        }
        let bin = members
            .iter()
            .zip(composition_distances(&members, lengths, composition))
            .zip(depth_offsets(&members, depths))
            .map(|((contig, composition), depth)| Measured {
                contig: *contig,
                carrier: carries[*contig],
                composition,
                depth,
            })
            .collect::<Vec<_>>();
        // Carriers are the bin's own, so only a contig further out than all of them is foreign.
        let spread = bin
            .iter()
            .filter(|measured| measured.carrier)
            .map(|measured| measured.composition)
            .fold(f64::MIN, f64::max);
        replicons.extend(
            bin.iter()
                .filter(|measured| !measured.carrier && measured.composition > spread)
                .map(|measured| measured.contig)
                .filter(is_shaped),
        );
        measured.extend(bin);
    }
    let composition_limits = limits(&measured, lengths, |measured| measured.composition);
    let depth_limits = limits(&measured, lengths, |measured| measured.depth);
    let (shaped_passengers, mut passengers): (Vec<_>, Vec<_>) = measured
        .iter()
        .filter(|measured| {
            let band = octave(lengths[measured.contig]);
            !measured.carrier
                && measured.composition > composition_limits[band]
                && measured.depth > depth_limits[band]
        })
        .map(|measured| measured.contig)
        .partition(is_shaped);
    replicons.extend(shaped_passengers);
    replicons.sort_unstable();
    replicons.dedup();
    passengers.sort_unstable();
    debug!(
        "{} small replicons and {} passengers",
        replicons.len(),
        passengers.len()
    );
    Departures {
        replicons,
        passengers,
    }
}
