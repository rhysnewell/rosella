use std::collections::HashSet;

use crate::embedding::features::ContigFeatures;
use crate::quality::ContigQuality;

pub const DISSOLVE_SELECT_NAMES: [&str; 2] = ["rounds", "ranked"];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Selection {
    Rounds,
    Ranked,
}

impl Selection {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "rounds" => Some(Self::Rounds),
            "ranked" => Some(Self::Ranked),
            _ => None,
        }
    }
}

/// No one k is right for every genome in a pool, so a round that runs first must not own the
/// contigs a later one would have grouped better. Rank every proposal and let the best claim first.
pub fn ranked(
    features: &ContigFeatures,
    quality: Option<&ContigQuality>,
    candidates: Vec<Vec<usize>>,
) -> Vec<(Vec<usize>, f64)> {
    let mut scored = candidates
        .into_iter()
        .map(|contigs| {
            let rank = worth(features, quality, &contigs);
            (contigs, rank)
        })
        .collect::<Vec<_>>();
    scored.sort_by(|one, other| {
        other
            .1
            .partial_cmp(&one.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| one.0.cmp(&other.0))
    });
    scored
}

fn worth(features: &ContigFeatures, quality: Option<&ContigQuality>, contigs: &[usize]) -> f64 {
    match quality {
        Some(quality) => quality.score(contigs).score(),
        None => features.bin_size(contigs) as f64,
    }
}

/// A proposal the winners already emptied is not the proposal that was scored, so what is left
/// of it is put back through the same bar rather than trusted on the rank it earned whole.
pub fn remaining(contigs: &[usize], claimed: &HashSet<usize>) -> Vec<usize> {
    contigs
        .iter()
        .copied()
        .filter(|contig| !claimed.contains(contig))
        .collect()
}
