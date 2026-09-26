use anyhow::Result;

use crate::clustering::clusterer::Partitioning;
use crate::recover::ladder::{Judge, best_per_arm, combine};
use crate::recover::recover_engine::RecoverEngine;

impl RecoverEngine {
    // Short contigs pull their genome's long contigs out of the bins they contaminate, but on
    // some assemblies they drag long contigs in. A long-only partition lets each bin refuse them.
    pub(super) fn attract(&self, main: Partitioning, contigs: &[usize]) -> Result<Partitioning> {
        let lengths = &self.coverage_table.contig_lengths;
        let long = contigs
            .iter()
            .copied()
            .filter(|contig| lengths[*contig] >= self.cutoff)
            .collect::<Vec<_>>();
        if long.len() == contigs.len() {
            return Ok(main);
        }
        let _timer = crate::timing::scope("attract");
        let (graph, _) = self.embed(&long);
        let mut ladder = Vec::new();
        for step in 0..self.partition_seeds {
            for mut held in self.partition_of(
                &graph,
                &long,
                self.partition,
                true,
                self.seeds.partition + step as u64,
            )? {
                held.reindex_clusters(&long);
                ladder.push(held);
            }
        }
        let judge = Judge {
            quality: &self.quality,
            contigs,
            bars: self.bars(),
        };
        let mut arms = best_per_arm(ladder, &judge);
        arms.push(main);
        Ok(combine(&arms, &judge, None))
    }
}
