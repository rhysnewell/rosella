use std::collections::{BTreeMap, HashMap, HashSet};
use std::ops::Range;

use ndarray::Axis;

use super::votes::{long_shares, with_long};
use super::{Down, Proposal, Walked, best_bins, long_votes, merge};
use crate::embedding::knn::KnnGraph;
use crate::recover::recover_engine::RecoverEngine;

struct Band {
    first: Vec<Proposal>,
    open: Vec<usize>,
    knn: KnnGraph,
    bar: f32,
}

impl RecoverEngine {
    pub(super) fn walk_down(
        &self,
        down: &Down,
        bins: &BTreeMap<usize, Vec<usize>>,
        spans: &[Range<usize>],
        bars: &[f32],
        among: &KnnGraph,
    ) -> Walked {
        let lengths = &self.tables.coverage.contig_lengths;
        let mut walked = Walked {
            floor: self.settings.lengths.cutoff,
            ..Walked::default()
        };
        let mut held = Vec::new();
        for (span, bar) in spans.iter().zip(bars) {
            let band = &down.order[span.clone()];
            let open = down.open(span.clone());
            let (proposals, long, knn) = self.propose(&open, down, among);
            let joined = joining(&proposals, *bar);
            walked.evidence.extend(self.marker_evidence(bins, &joined));
            walked.joins.extend(joined);
            if self.stages.is_some() {
                walked.votes.extend(with_long(proposals.clone(), long));
            }
            held.push(Band {
                first: proposals,
                open,
                knn,
                bar: *bar,
            });
            walked.floor = band
                .iter()
                .map(|contig| lengths[*contig])
                .fold(walked.floor, usize::min);
        }
        self.vote_again(down, bins, &held, &mut walked);
        walked
    }

    // Marker repeats set the share bar, and a eukaryote barely carries markers. So a contig placed
    // in a eukaryotic bin goes unchecked and votes for no neighbour in a later round.
    fn vote_again(
        &self,
        down: &Down,
        bins: &BTreeMap<usize, Vec<usize>>,
        held: &[Band],
        walked: &mut Walked,
    ) {
        let eukaryotic = self
            .quality
            .eukaryotic(
                bins.values().map(Vec::as_slice),
                self.settings.bars.quality(),
            )
            .into_iter()
            .zip(bins.keys())
            .filter_map(|(called, bin)| called.then_some(*bin))
            .collect::<HashSet<_>>();
        let mut placed = walked
            .joins
            .iter()
            .map(|(contig, _)| *contig)
            .collect::<HashSet<_>>();
        let mut voters = down
            .bin_of
            .iter()
            .map(|(contig, bin)| (*contig, *bin))
            .chain(
                walked
                    .joins
                    .iter()
                    .filter(|(_, bin)| !eukaryotic.contains(bin))
                    .map(|(contig, bin)| (*contig, *bin)),
            )
            .collect::<HashMap<_, _>>();
        let mut rounds = 0;
        loop {
            let fresh = held
                .iter()
                .flat_map(|band| {
                    let best = best_bins(
                        &band.knn,
                        |neighbour| voters.get(&neighbour).copied(),
                        |row| down.ballots.get(band.open[row]),
                    );
                    band.first
                        .iter()
                        .zip(best)
                        .filter(|((contig, _), _)| !placed.contains(contig))
                        .filter_map(|((contig, _), best)| {
                            best.filter(|(_, share)| *share > band.bar)
                                .map(|(bin, _)| (*contig, bin))
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            if fresh.is_empty() {
                break;
            }
            rounds += 1;
            walked.evidence.extend(self.marker_evidence(bins, &fresh));
            placed.extend(fresh.iter().map(|(contig, _)| *contig));
            voters.extend(
                fresh
                    .iter()
                    .filter(|(_, bin)| !eukaryotic.contains(bin))
                    .map(|(contig, bin)| (*contig, *bin)),
            );
            walked.joins.extend(fresh);
        }
        log::debug!(
            "Attach rounds: {rounds} after the first, {} joins in all",
            walked.joins.len()
        );
    }

    fn propose(
        &self,
        open: &[usize],
        down: &Down,
        among: &KnnGraph,
    ) -> (Vec<Proposal>, Vec<f32>, KnnGraph) {
        let first = down.long_graph.indices.nrows();
        let band = open.iter().map(|at| down.order[*at]).collect::<Vec<_>>();
        let nearest = self.nearest_long(down.long_graph, &band);
        let _timer = crate::timing::scope("attach");
        let own = KnnGraph {
            indices: among.indices.select(Axis(0), open),
            dists: among.dists.select(Axis(0), open),
        };
        let knn = merge(&nearest, &own, band.len(), |at| down.order[at]);
        let long = match self.stages {
            Some(_) => long_shares(&knn, first, down.bin_of),
            None => Vec::new(),
        };
        let proposals = band
            .into_iter()
            .zip(best_bins(&knn, long_votes(first, down.bin_of), |row| {
                down.ballots.get(open[row])
            }))
            .collect();
        (proposals, long, knn)
    }
}

pub(super) fn joining(proposals: &[Proposal], bar: f32) -> Vec<(usize, usize)> {
    proposals
        .iter()
        .filter_map(|(contig, best)| {
            let (bin, share) = (*best)?;
            (share > bar).then_some((*contig, bin))
        })
        .collect()
}
