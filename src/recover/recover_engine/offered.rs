use log::debug;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::Graph;
use crate::embedding::knn::KnnGraph;
use crate::progress::{self, Stage};
use crate::recover::census::Census;
use crate::recover::ladder::{Judge, combine_traced};
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::splitter::Refined;

pub const PARTITION: &str = "partition";

pub(super) struct SetRow {
    pub name: String,
    pub offered: (usize, usize),
    pub refined: Option<(usize, usize)>,
    pub whole: usize,
    pub trimmed: usize,
}

impl RecoverEngine {
    // A weaker binner's bin can beat a whole genome on one doubled marker in its remnant, so every
    // offered set goes through the refine cycle and the bins are arbitrated again.
    pub(super) fn refine_offered(
        &self,
        settled: Partitioning,
        graph: &Graph,
        knn: &KnnGraph,
        contigs: &[usize],
    ) -> (Partitioning, Vec<SetRow>) {
        if self.candidates.is_empty() {
            return (settled, Vec::new());
        }
        let lengths = self.lengths();
        let cutoff = self.settings.lengths.cutoff;
        let long = &contigs[..contigs.partition_point(|contig| lengths[*contig] >= cutoff)];
        let held = |bins: &Partitioning| {
            let bp = bins
                .cluster_map
                .values()
                .flatten()
                .map(|at| lengths[long[*at]]);
            (bins.cluster_map.len(), bp.sum::<usize>())
        };
        let mut rows = vec![SetRow {
            name: PARTITION.to_string(),
            offered: held(&settled),
            refined: None,
            whole: 0,
            trimmed: 0,
        }];
        let mut arms = vec![settled];
        let progress = progress::hold(progress::counted(
            Stage::RefiningBins,
            self.candidates.len() as u64,
        ));
        progress.set_message("offered sets");
        for (set, offered) in self.candidates.iter().zip(self.candidate_arms(long)) {
            let given = held(&offered);
            let refined = self.refine_clusters(offered, graph, knn, &mut Census::default());
            progress.inc(1);
            let refined = Partitioning {
                cluster_map: refined.bins,
                outliers: refined.unbinned,
                score: None,
                arm: self.settings.partition,
                seed: self.settings.seeds.partition,
            };
            rows.push(SetRow {
                name: set.name.clone(),
                offered: given,
                refined: Some(held(&refined)),
                whole: 0,
                trimmed: 0,
            });
            arms.push(refined);
        }
        progress::done(&progress);
        let judge = Judge {
            quality: &self.quality,
            contigs: long,
            lengths,
            bars: self.settings.bars,
        };
        // Refining every set the same way makes their cuts agree, so each set's raw bins vote
        // beside its refined ones.
        let sets = self.candidates.len();
        let voters = (0..=sets).chain(1..=sets).collect::<Vec<_>>();
        arms.extend(self.candidate_arms(long));
        let combined = combine_traced(&arms, &voters, &judge, None);
        for (offered_by, whole) in &combined.sources {
            let mut by = offered_by
                .iter()
                .map(|arm| voters[*arm])
                .collect::<Vec<_>>();
            by.sort_unstable();
            by.dedup();
            for voter in by {
                match whole {
                    true => rows[voter].whole += 1,
                    false => rows[voter].trimmed += 1,
                }
            }
        }
        for row in &rows {
            debug!(
                "Ensemble: {} gave {} bins whole and {} trimmed",
                row.name, row.whole, row.trimmed
            );
        }
        (combined.partitioning, rows)
    }

    // The refine cycle reads markers alone, which cannot see a genome's own marker-free contigs,
    // so its bins go back to the vote beside every offered set.
    pub(super) fn revote(&self, refined: &mut Refined, contigs: &[usize]) {
        let lengths = self.lengths();
        let cutoff = self.settings.lengths.cutoff;
        let long = &contigs[..contigs.partition_point(|contig| lengths[*contig] >= cutoff)];
        let mut arms = vec![Partitioning {
            cluster_map: std::mem::take(&mut refined.bins),
            outliers: std::mem::take(&mut refined.unbinned),
            score: None,
            arm: self.settings.partition,
            seed: self.settings.seeds.partition,
        }];
        arms.extend(self.candidate_arms(long));
        let judge = Judge {
            quality: &self.quality,
            contigs,
            lengths,
            bars: self.settings.bars,
        };
        let voters = (0..arms.len()).collect::<Vec<_>>();
        let combined = combine_traced(&arms, &voters, &judge, None).partitioning;
        refined.bins = combined.cluster_map;
        refined.unbinned = combined.outliers;
    }
}
