use std::collections::{BTreeMap, HashSet};

use ndarray::ArrayView2;
use rayon::prelude::*;

use crate::coverage::read_stats::tie_share;
use crate::embedding::knn::KnnGraph;
use crate::markers::replicon::LEAST_BASES;
use crate::refine::owners::owners;
use crate::refine::report_context::{neighbour_weight, share};

// The partition labels every contig, so a short one whose neighbours mostly sit in other bins
// was placed on weak evidence.
const ADRIFT: f64 = 0.5;

// Two strains in one bin tie the reads of the sequence they share, and the copy away from the
// bin's depth belongs to the other strain. A repeat inside one genome ties too, at one depth.
const TWOFOLD: f64 = 1.0;
const TIED: f64 = 0.25;
const LEAST_ANCHORS: usize = 3;
const DEPTH_FLOOR: f64 = 0.1;

pub struct Tie {
    pub first: usize,
    pub second: usize,
    pub reads: u64,
}

fn logs(depths: ArrayView2<f64>, row: usize) -> Vec<f64> {
    depths
        .row(row)
        .iter()
        .map(|depth| (depth + DEPTH_FLOOR).log2())
        .collect()
}

fn apart(left: &[f64], right: &[f64]) -> f64 {
    let squares = left
        .iter()
        .zip(right)
        .map(|(left, right)| (left - right).powi(2))
        .sum::<f64>();
    (squares / left.len().max(1) as f64).sqrt()
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_unstable_by(f64::total_cmp);
    let middle = values.len() / 2;
    match values.len() % 2 {
        0 => (values[middle - 1] + values[middle]) / 2.0,
        _ => values[middle],
    }
}

// Only long contigs hold enough reads for their depth to stand for the genome.
fn offsets<'a>(
    bins: impl IntoIterator<Item = &'a [usize]>,
    lengths: &[usize],
    depths: ArrayView2<f64>,
) -> Vec<Option<(usize, f64)>> {
    let mut offsets = vec![None; lengths.len()];
    for (bin, members) in bins.into_iter().enumerate() {
        let anchors = members
            .iter()
            .filter(|row| lengths[**row] >= LEAST_BASES)
            .map(|row| logs(depths, *row))
            .collect::<Vec<_>>();
        if anchors.len() < LEAST_ANCHORS {
            continue;
        }
        let reference = (0..depths.ncols())
            .map(|sample| median(anchors.iter().map(|logs| logs[sample]).collect()))
            .collect::<Vec<_>>();
        for row in members {
            offsets[*row] = Some((bin, apart(&logs(depths, *row), &reference)));
        }
    }
    offsets
}

pub fn strays<'a>(
    bins: impl IntoIterator<Item = &'a [usize]>,
    lengths: &[usize],
    depths: ArrayView2<f64>,
    ties: &[Tie],
    reads: &[u64],
) -> Vec<usize> {
    let offsets = offsets(bins, lengths, depths);
    let mut strays = ties
        .iter()
        .filter_map(|tie| {
            let (bin, first) = offsets[tie.first]?;
            let (other, second) = offsets[tie.second]?;
            let gap = apart(&logs(depths, tie.first), &logs(depths, tie.second));
            if bin != other || gap < TWOFOLD {
                return None;
            }
            let out = if first > second {
                tie.first
            } else {
                tie.second
            };
            (lengths[out] < LEAST_BASES && tie_share(tie.reads, reads[out]) >= TIED).then_some(out)
        })
        .collect::<Vec<_>>();
    strays.sort_unstable();
    strays.dedup();
    strays
}

// A bin that already clears the bars gains nothing by losing a contig, and its own short contigs
// are the ones its completeness rests on.
pub fn adrift(
    bins: &BTreeMap<usize, Vec<usize>>,
    finished: &HashSet<usize>,
    lengths: &[usize],
    knn: &KnnGraph,
) -> Vec<usize> {
    let owner = owners(bins.iter().map(|(bin, members)| (*bin, members)));
    let mut adrift = owner
        .par_iter()
        .filter(|(contig, bin)| {
            lengths[**contig] < LEAST_BASES && !finished.contains(*bin) && {
                let mut weights = neighbour_weight(**contig, &owner, knn, lengths);
                share(&mut weights, **bin).is_some_and(|share| share < ADRIFT)
            }
        })
        .map(|(contig, _)| *contig)
        .collect::<Vec<_>>();
    adrift.sort_unstable();
    adrift
}
