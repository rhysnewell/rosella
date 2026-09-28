use std::collections::HashMap;

use log::{info, warn};
use rand::{Rng, SeedableRng};

use crate::embedding::knn::{KnnGraph, nearest_in};
use crate::embedding::metrics::MIN_VAR;
use crate::recover::homing::{Sample, chances};
use crate::recover::recover_engine::RecoverEngine;
use crate::recover::recover_engine::attach::{best_bins, merge};

impl RecoverEngine {
    pub(super) fn home_chances(
        &self,
        bin_of: &HashMap<usize, usize>,
        proposals: &[(usize, Option<(usize, f32)>)],
        long_graph: &KnnGraph,
        among: &KnnGraph,
        band: &[usize],
    ) -> Option<Vec<f64>> {
        let first = long_graph.indices.nrows();
        let _timer = crate::timing::scope("calibrate");
        let lengths = &self.coverage_table.contig_lengths;
        let binned = (0..first)
            .filter(|contig| bin_of.contains_key(contig))
            .collect::<Vec<_>>();
        if binned.is_empty() || proposals.is_empty() {
            return None;
        }
        let drawn = crate::seeds::sample_positions(
            proposals.len(),
            super::weight::SAMPLE.min(proposals.len()),
            self.seeds.seed,
        );
        let parents = crate::seeds::sample_positions(
            binned.len(),
            drawn.len().min(binned.len()),
            self.seeds.seed.wrapping_add(crate::defaults::SEED_STRIDE),
        );
        let pieces = drawn
            .iter()
            .enumerate()
            .map(|(at, row)| {
                (
                    binned[parents[at % parents.len()]],
                    lengths[proposals[*row].0],
                )
            })
            .collect::<Vec<_>>();

        let names = &self.coverage_table.contig_names;
        let composition = match crate::kmers::kmer_counting::prefixes(
            &self.assembly,
            &pieces
                .iter()
                .map(|(parent, length)| (names[*parent].as_str(), *length))
                .collect::<Vec<_>>(),
            &self.tnf_table.kmer_sizes(),
        ) {
            Ok(composition) => composition,
            Err(error) => {
                warn!("Attaching uncalibrated, the pieces could not be read: {error}");
                return None;
            }
        };
        let coverage = self
            .coverage_table
            .table
            .select(ndarray::Axis(0), &self.depth_donors(&pieces, bin_of, first));
        let all = (0..self.n_contigs).collect::<Vec<_>>();
        let mut metric = self.features().prepared(&all);
        metric.extend(
            &coverage,
            &composition,
            MIN_VAR,
            &pieces.iter().map(|(_, length)| *length).collect::<Vec<_>>(),
        );
        let base = self.n_contigs;

        let long_without = |skip: &(dyn Fn(usize, usize) -> bool + Sync)| {
            nearest_in(
                long_graph,
                pieces.len(),
                long_graph.indices.ncols(),
                self.knn_candidates,
                self.seeds.knn,
                |piece, other| match skip(piece, other) {
                    true => f64::INFINITY,
                    false => metric.distance(base + piece, other),
                },
            )
        };
        let beside_parent = long_without(&|piece, other| other == pieces[piece].0);
        let without_bin =
            long_without(&|piece, other| bin_of.get(&other) == bin_of.get(&pieces[piece].0));
        let parked = nearest_in(
            among,
            pieces.len(),
            among.indices.ncols(),
            self.knn_candidates,
            self.seeds.knn,
            |piece, other| metric.distance(base + piece, band[other]),
        );
        let shares = |long: &KnnGraph| {
            best_bins(
                &merge(long, |row| row, &parked, pieces.len(), |at| band[at]),
                first,
                bin_of,
            )
        };
        let sample = |best: &Option<(usize, f32)>, length: usize| Sample {
            share: best.map_or(0.0, |(_, share)| f64::from(share)),
            length,
        };
        let homed = shares(&beside_parent)
            .iter()
            .zip(&pieces)
            .map(|(best, (parent, length))| {
                let right = best.is_some_and(|(bin, _)| bin_of.get(parent) == Some(&bin));
                (sample(best, *length), right)
            })
            .collect::<Vec<_>>();
        let homeless = shares(&without_bin)
            .iter()
            .zip(&pieces)
            .map(|(best, (_, length))| sample(best, *length))
            .collect::<Vec<_>>();
        let real = proposals
            .iter()
            .map(|(contig, best)| sample(best, lengths[*contig]))
            .collect::<Vec<_>>();
        let fit = chances(&homed, &homeless, &real)?;
        info!(
            "Attach calibrated on {} pieces: pieces return home {:.3} beside their bin, the share \
             of parked contigs with a home by length {:?}, {} over one half.",
            pieces.len(),
            homed.iter().filter(|(_, right)| *right).count() as f64 / homed.len() as f64,
            fit.prior
                .iter()
                .map(|(length, prior)| format!("{length} bp {prior:.3}"))
                .collect::<Vec<_>>(),
            fit.chances.iter().filter(|chance| **chance > 0.5).count(),
        );
        Some(fit.chances)
    }

    // A real short contig never carries exactly its neighbour's depth, so each piece takes the
    // depth of another long contig in its parent's bin and with it the bin's own depth spread.
    fn depth_donors(
        &self,
        pieces: &[(usize, usize)],
        bin_of: &HashMap<usize, usize>,
        first: usize,
    ) -> Vec<usize> {
        let mut members = HashMap::<usize, Vec<usize>>::new();
        for contig in 0..first {
            if let Some(bin) = bin_of.get(&contig) {
                members.entry(*bin).or_default().push(contig);
            }
        }
        let mut rng = rand::rngs::StdRng::seed_from_u64(self.seeds.seed);
        pieces
            .iter()
            .map(|(parent, _)| {
                let others = &members[&bin_of[parent]];
                match others.len() > 1 {
                    true => loop {
                        let donor = others[rng.random_range(0..others.len())];
                        if donor != *parent {
                            break donor;
                        }
                    },
                    false => *parent,
                }
            })
            .collect()
    }
}
