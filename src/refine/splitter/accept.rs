use std::collections::HashSet;

use log::debug;

use super::{Refiner, Round};
use crate::clustering::clusterer::{Partitioning, find_best_partition};
use crate::refine::bar::{MIN_SPLIT_CONTIGS, should_split};
use crate::refine::bin_stats::{AGGREGATE, BinStats, bin_stats};
use crate::refine::gates::{SplitRejection, Trigger};
use crate::refine::proposal::{
    Proposal, SplitOutcome, contigs, judge_split, leaves_two_standing, tighter,
};
use crate::refine::{bisect, peel};

impl Refiner<'_> {
    pub(super) fn propose(&self, bin_id: usize, round: &Round) -> Proposal {
        let indices = &self.refined.bins[&bin_id];
        let Some(stats) = self.cached.get(&bin_id) else {
            return Proposal::TooFewContigs;
        };
        let lengths = self.features.contig_lengths(indices);

        let partition = match self.settings.split {
            true => self.propose_partition(bin_id, indices, stats, &lengths, round),
            false => Proposal::NoTrigger,
        };
        if matches!(partition, Proposal::Accepted(..)) || !self.settings.peel {
            return partition;
        }
        let peeled = {
            let _timer = crate::timing::scope("peel");
            self.peel(indices, stats, &lengths)
        };
        match peeled {
            Some(outcome) => {
                debug!(
                    "Peel: {} contigs off bin {} of {}",
                    outcome.kept.len() - 1,
                    bin_id,
                    indices.len()
                );
                Proposal::Accepted(Trigger::Peeled, outcome)
            }
            None => partition,
        }
    }

    fn propose_partition(
        &self,
        bin_id: usize,
        indices: &[usize],
        stats: &BinStats,
        lengths: &[usize],
        round: &Round,
    ) -> Proposal {
        let bin_size = lengths.iter().sum::<usize>();
        let over_budget = match (
            self.judges.contamination.get(&bin_id),
            self.settings.max_contamination,
        ) {
            (Some(contamination), Some(budget)) => *contamination > budget,
            _ => false,
        };
        let trigger = should_split(
            stats,
            lengths,
            bin_size,
            over_budget,
            self.settings.max_bin_size,
            round.thresholds,
        );
        let Some(trigger) = trigger else {
            return match indices.len() < MIN_SPLIT_CONTIGS {
                true => Proposal::TooFewContigs,
                false => Proposal::NoTrigger,
            };
        };
        let Some(result) = self.cluster_bin(indices) else {
            return Proposal::NoClustering(trigger);
        };
        debug!(
            "Bin {} of {} contigs: {} clusters, validity {:.3}",
            bin_id,
            indices.len(),
            result.cluster_map.len(),
            result.score.unwrap_or_default()
        );
        let accepted = {
            let _timer = crate::timing::scope("accept");
            self.accept(indices, stats, result, round)
        };
        match accepted {
            Ok(outcome) => {
                debug!(
                    "Split: bin {} of {} contigs into {}",
                    bin_id,
                    indices.len(),
                    outcome.kept.len()
                );
                Proposal::Accepted(trigger, outcome)
            }
            Err(rejection) => Proposal::Rejected(trigger, rejection),
        }
    }

    fn stands(&self, outcome: &SplitOutcome) -> bool {
        let floor = self.settings.piece;
        debug!(
            "Cut: pieces {:?}, scattered {}, floor {floor}",
            outcome
                .kept
                .iter()
                .map(|piece| self.features.bin_size(piece))
                .collect::<Vec<_>>(),
            self.features.bin_size(&outcome.unbinned)
        );
        leaves_two_standing(&outcome.kept, floor, |piece| self.features.bin_size(piece))
    }

    // A lone contig counts as no spread, since a genome in one contig has none.
    fn peel(&self, indices: &[usize], stats: &BinStats, lengths: &[usize]) -> Option<SplitOutcome> {
        let peel = peel::candidate(indices, stats, lengths, self.settings.piece)?;
        let rest = self.stats_of(&peel.rest)?;
        let rest_bp = self.features.bin_size(&peel.rest) as f64;
        let lone_bp = self.features.bin_size(&peel.lone) as f64;
        if !tighter(rest.mean[AGGREGATE] * rest_bp / (rest_bp + lone_bp), stats) {
            return None;
        }
        let mut kept = peel
            .lone
            .into_iter()
            .map(|contig| vec![contig])
            .collect::<Vec<_>>();
        let keeps = self.keeps_a_genome(&peel.rest, &kept);
        let mut pieces = kept.iter().map(|_| None).collect::<Vec<_>>();
        kept.push(peel.rest);
        pieces.push(Some(rest));
        keeps.then_some(SplitOutcome {
            kept,
            stats: pieces,
            unbinned: Vec::new(),
        })
    }

    // Families shared with the peeled contigs are a single copy gene twice over, so the two
    // pieces are two organisms. Families that only complement them are one genome coming apart.
    fn keeps_a_genome(&self, rest: &[usize], lone: &[Vec<usize>]) -> bool {
        let Some(quality) = self.judges.quality else {
            return true;
        };
        // A remainder under the bin floor is dust. It carries no families to judge it by.
        if self.features.bin_size(rest) < self.settings.piece {
            return true;
        }
        let held = quality.features(rest);
        let peeled = lone
            .iter()
            .flat_map(|piece| quality.features(piece))
            .collect::<HashSet<_>>();
        peeled.is_empty() || !peeled.is_disjoint(&held)
    }

    fn cluster_bin(&self, indices: &[usize]) -> Option<Partitioning> {
        let seeds = self.settings.seeds;
        let graph = crate::refine::split_graph::bin_graph(
            &self.features,
            self.judges.assembly,
            indices,
            self.settings.knn.neighbours,
            seeds.knn,
            self.settings.knn.candidates,
        );
        find_best_partition(
            &graph,
            &self.features.contig_lengths(indices),
            seeds.partition,
            self.settings.partition.for_split(),
        )
        .ok()
    }

    fn accept(
        &self,
        indices: &[usize],
        stats: &BinStats,
        result: Partitioning,
        round: &Round,
    ) -> Result<SplitOutcome, SplitRejection> {
        let mut clusters = result
            .cluster_map
            .into_values()
            .map(|positions| contigs(indices, positions.into_iter()))
            .collect::<Vec<_>>();
        clusters.sort_unstable();
        let noise = contigs(indices, result.outliers.into_iter());

        let (kept, spare) =
            judge_split(clusters, noise, |cluster| self.features.bin_size(cluster))?;

        let pieces = kept
            .iter()
            .map(|piece| self.stats_of(piece))
            .collect::<Vec<_>>();
        if !self.pieces_are_tighter(&kept, &pieces, stats) {
            return Err(SplitRejection::NotTighter);
        }

        let outcome =
            self.place_leftovers(kept, pieces, spare, round.thresholds.quantile[AGGREGATE]);
        if !self.stands(&outcome) {
            return Err(SplitRejection::Shredded);
        }

        if !bisect::separates(
            &self.features,
            &outcome.kept,
            round.eligible,
            self.settings.seeds.seed,
        ) {
            return Err(SplitRejection::Unimodal);
        }
        Ok(outcome)
    }

    fn stats_of(&self, contigs: &[usize]) -> Option<BinStats> {
        bin_stats(&self.features, contigs, self.settings.seeds.seed)
    }

    // A chimeric bin falls apart into tighter pieces. A pure one does not.
    fn pieces_are_tighter(
        &self,
        kept: &[Vec<usize>],
        pieces: &[Option<BinStats>],
        whole: &BinStats,
    ) -> bool {
        let mut weighted = 0.0;
        let mut total = 0;
        for (piece, stats) in kept.iter().zip(pieces) {
            let Some(stats) = stats else {
                continue;
            };
            let size = self.features.bin_size(piece);
            weighted += stats.mean[AGGREGATE] * size as f64;
            total += size;
        }
        if total == 0 {
            return false;
        }

        tighter(weighted / total as f64, whole)
    }

    fn place_leftovers(
        &self,
        mut kept: Vec<Vec<usize>>,
        mut pieces: Vec<Option<BinStats>>,
        spare: Vec<usize>,
        level: f64,
    ) -> SplitOutcome {
        let spared = self.stats_of(&spare);
        let holds_together = spared
            .as_ref()
            .is_some_and(|stats| stats.mean[AGGREGATE] < level);

        if holds_together {
            kept.push(spare);
            pieces.push(spared);
            return SplitOutcome {
                kept,
                stats: pieces,
                unbinned: Vec::new(),
            };
        }

        SplitOutcome {
            kept,
            stats: pieces,
            unbinned: spare,
        }
    }
}
