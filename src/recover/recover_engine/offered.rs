use crate::clustering::clusterer::Partitioning;
use crate::recover::ladder::{Judge, combine_traced};
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::splitter::Refined;

pub const PARTITION: &str = "partition";

impl RecoverEngine {
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
