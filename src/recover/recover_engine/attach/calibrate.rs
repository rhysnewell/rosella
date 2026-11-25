use std::ops::Range;
use std::sync::LazyLock;

use anyhow::Result;
use log::debug;
use ndarray::{Array2, Axis, concatenate, s};
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use super::Down;
use crate::embedding::knn::{KnnGraph, nearest_in};
use crate::embedding::metrics::prepared::PreparedAggregate;
use crate::learned::windows::Windows;
use crate::recover::offered_votes::Ballot;
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::owners::credit;

pub(super) static GAP: LazyLock<bool> =
    LazyLock::new(|| std::env::var_os("ROSELLA_GAP_PLACEMENT").is_some());

const PARENTS: usize = 2000;

type Cut = (usize, usize, Range<usize>);

impl RecoverEngine {
    // A piece cut from a binned long contig is a short contig whose home is known, so each band's
    // cutoff is read from pieces cut to that band's own lengths.
    pub(super) fn gap_cutoffs(
        &self,
        down: &Down,
        spans: &[Range<usize>],
    ) -> Result<Option<Vec<f32>>> {
        let Some(dir) = self.tables.windows.as_deref() else {
            return Ok(None);
        };
        let first = down.long_graph.indices.nrows();
        let lengths = &self.tables.coverage.contig_lengths;
        let names = &self.tables.coverage.contig_names;
        let bands = spans
            .iter()
            .map(|span| {
                down.order[span.clone()]
                    .iter()
                    .map(|contig| lengths[*contig])
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let Some(longest) = bands.iter().flatten().max().copied() else {
            return Ok(None);
        };
        let pool = (0..first)
            .filter(|contig| down.bin_of.contains_key(contig) && lengths[*contig] >= 2 * longest)
            .collect::<Vec<_>>();
        if pool.is_empty() {
            return Ok(None);
        }
        let seed = self.settings.seeds.seed;
        let mut rng = StdRng::seed_from_u64(seed);
        let parents = crate::seeds::sample_positions(pool.len(), PARENTS, seed)
            .into_iter()
            .map(|at| pool[at])
            .collect::<Vec<_>>();
        let mut cuts = Vec::<Cut>::with_capacity(parents.len() * bands.len());
        for parent in parents {
            for (band, sizes) in bands.iter().enumerate() {
                let length = sizes[rng.random_range(0..sizes.len())];
                let start = rng.random_range(0..=lengths[parent] - length);
                cuts.push((parent, band, start..start + length));
            }
        }
        let stretches = cuts
            .iter()
            .map(|(parent, _, span)| (names[*parent].as_str(), span.clone()))
            .collect::<Vec<_>>();
        let composition = crate::kmers::kmer_counting::stretches(
            &self.tables.assembly,
            &stretches,
            self.tables.tnf.kmer_size(),
        )?;
        let coverage =
            Windows::read(dir, &self.tables.coverage.sample_names)?.stretch_rows(&stretches)?;
        let composition = match &self.tables.learned {
            Some(learned) => {
                let named = stretches
                    .iter()
                    .map(|(name, _)| name.to_string())
                    .collect::<Vec<_>>();
                learned.embed(&named, composition.view(), coverage.view())?
            }
            None => composition,
        };
        let knn = self.pieces_nearest(down, &cuts, (composition, coverage))?;
        let scores = gaps(
            &knn,
            |neighbour| down.bin_of.get(&neighbour).copied(),
            |_| true,
            |_| None,
        );
        let mut by_band = vec![Vec::new(); bands.len()];
        for ((parent, band, _), score) in cuts.iter().zip(scores) {
            let home = score.is_some_and(|(bin, _)| down.bin_of.get(parent) == Some(&bin));
            by_band[*band].push((score.map_or(0.0, |(_, gap)| gap), home));
        }
        let worth = self.settings.bars.worth;
        Ok(Some(
            by_band
                .into_iter()
                .zip(&bands)
                .map(|(mut pieces, sizes)| {
                    let (cutoff, home) = cutoff(&mut pieces, worth);
                    debug!(
                        "Attach gap from {} bp: cutoff {cutoff:.3}, {home} of {} pieces home above it",
                        sizes.iter().min().unwrap_or(&0),
                        pieces.len()
                    );
                    cutoff
                })
                .collect(),
        ))
    }

    fn pieces_nearest(
        &self,
        down: &Down,
        cuts: &[Cut],
        (composition, coverage): (Array2<f64>, Array2<f64>),
    ) -> Result<KnnGraph> {
        let first = down.long_graph.indices.nrows();
        let tnf = concatenate(
            Axis(0),
            &[
                self.tables.tnf.kmer_table.slice(s![..first, ..]),
                composition.view(),
            ],
        )?;
        let depth = concatenate(
            Axis(0),
            &[
                self.tables.coverage.table.slice(s![..first, ..]),
                coverage.view(),
            ],
        )?;
        let rows = (0..first + cuts.len()).collect::<Vec<_>>();
        let prepared =
            PreparedAggregate::new(&depth, &tnf, &rows, &self.settings.features.distance);
        let k = down.long_graph.indices.ncols();
        let found = prepared.device_nearest(first, k + 1).unwrap_or_else(|| {
            nearest_in(
                down.long_graph,
                cuts.len(),
                k + 1,
                self.settings.knn.candidates,
                self.settings.seeds.knn,
                prepared.shifted(first),
            )
        });
        Ok(without_parents(&found, cuts, k))
    }
}

fn without_parents(found: &KnnGraph, cuts: &[Cut], k: usize) -> KnnGraph {
    let mut indices = Array2::from_elem((cuts.len(), k), u32::MAX);
    let mut dists = Array2::from_elem((cuts.len(), k), f32::INFINITY);
    for (row, (parent, _, _)) in cuts.iter().enumerate() {
        let (indices_row, dists_row) = (found.indices.row(row), found.dists.row(row));
        let kept = indices_row
            .iter()
            .zip(dists_row)
            .filter(|(at, _)| **at as usize != *parent)
            .take(k);
        for (slot, (at, distance)) in kept.enumerate() {
            indices[[row, slot]] = *at;
            dists[[row, slot]] = *distance;
        }
    }
    KnnGraph { indices, dists }
}

// A piece that goes home gains one and a stray costs the run's contamination worth, so the cutoff
// is the gap below which the pieces stop paying.
fn cutoff(pieces: &mut [(f32, bool)], worth: f64) -> (f32, usize) {
    pieces.sort_by(|a, b| b.0.total_cmp(&a.0));
    let (mut total, mut best, mut taken) = (0.0, 0.0, 0);
    for (at, (_, home)) in pieces.iter().enumerate() {
        total += if *home { 1.0 } else { -worth };
        if total > best {
            (best, taken) = (total, at + 1);
        }
    }
    let home = pieces[..taken].iter().filter(|(_, home)| *home).count();
    let cutoff = match taken {
        0 => f32::INFINITY,
        _ => pieces.get(taken).map_or(0.0, |(gap, _)| *gap),
    };
    (cutoff, home)
}

pub(super) fn gaps<'a>(
    knn: &KnnGraph,
    vote: impl Fn(usize) -> Option<usize> + Sync,
    counted: impl Fn(usize) -> bool + Sync,
    ballot: impl Fn(usize) -> Option<&'a Ballot> + Sync,
) -> Vec<Option<(usize, f32)>> {
    let (sigmas, rhos) = crate::embedding::fuzzy::scales(knn.dists.view());
    (0..knn.indices.nrows())
        .into_par_iter()
        .map(|row| {
            let mut mass = Vec::<(usize, f32)>::new();
            let mut total = 0.0;
            for (neighbour, distance) in knn.indices.row(row).iter().zip(knn.dists.row(row)) {
                if *neighbour == u32::MAX {
                    break;
                }
                let neighbour = *neighbour as usize;
                if !counted(neighbour) {
                    continue;
                }
                let weight = crate::embedding::fuzzy::membership(*distance, rhos[row], sigmas[row]);
                total += weight;
                if let Some(bin) = vote(neighbour) {
                    credit(&mut mass, bin, weight);
                }
            }
            let scale = f32::max(total, f32::MIN_POSITIVE);
            let voters = ballot(row).map_or(0, |ballot| {
                for bin in &ballot.bins {
                    credit(&mut mass, *bin, scale);
                }
                ballot.voters
            });
            mass.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
            let (bin, best) = *mass.first()?;
            let second = mass.get(1).map_or(0.0, |(_, held)| *held);
            Some((bin, (best - second) / (scale * (1 + voters) as f32)))
        })
        .collect()
}
