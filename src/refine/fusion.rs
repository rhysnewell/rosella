use crate::embedding::features::ContigFeatures;
use crate::refine::bar::MIN_SPLIT_CONTIGS;
use crate::refine::gates::{SplitRejection, Trigger};

/// A cut on a fused bin has to leave the host genome standing at t1 completeness: the
/// largest piece keeps this share of the markers the whole bin held.
const KEEPS: f64 = 0.95;

/// And it has to be the cut the markers asked for: a strain pair cut anywhere else leaves
/// the second copies together.
const RESOLVES: f64 = 0.9;

pub fn fused(features: &ContigFeatures, indices: &[usize], bar: f64) -> Option<Trigger> {
    if indices.len() < MIN_SPLIT_CONTIGS {
        return None;
    }
    features
        .markers()?
        .fusion(indices)
        .filter(|fusion| fusion.fraction() >= bar)
        .map(|_| Trigger::Fused)
}

pub fn parts_on_markers(
    features: &ContigFeatures,
    indices: &[usize],
    pieces: &[Vec<usize>],
) -> Result<(), SplitRejection> {
    let Some(markers) = features.markers() else {
        return Err(SplitRejection::Scattered);
    };
    let whole = markers.distinct(indices);
    let largest = pieces
        .iter()
        .map(|piece| markers.distinct(piece))
        .max()
        .unwrap_or(0);
    if (largest as f64) < KEEPS * whole as f64 {
        return Err(SplitRejection::Scattered);
    }
    let before = markers.fusion(indices).map_or(0, |fusion| fusion.duplicated);
    let after = pieces
        .iter()
        .filter_map(|piece| markers.fusion(piece))
        .map(|fusion| fusion.duplicated)
        .sum::<usize>();
    match (after as f64) <= (1.0 - RESOLVES) * before as f64 {
        true => Ok(()),
        false => Err(SplitRejection::Unresolved),
    }
}
