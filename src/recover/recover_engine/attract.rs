use anyhow::Result;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::{Graph, knn::KnnGraph};
use crate::recover::ladder::{Judge, combine};
use crate::recover::recover_engine::RecoverEngine;

impl RecoverEngine {
    // Short contigs pull long contigs out of bins they contaminate, or on some assemblies drag
    // them in, so the long-only partition competes per bin. Short contigs make the weight wander.
    pub(super) fn partitioned(
        &mut self,
        contigs: &[usize],
    ) -> Result<(Graph, KnnGraph, Partitioning)> {
        let lengths = &self.coverage_table.contig_lengths;
        let long = contigs
            .iter()
            .copied()
            .filter(|contig| lengths[*contig] >= self.cutoff)
            .collect::<Vec<_>>();
        if long.len() == contigs.len() {
            return self.weighted_partition(contigs);
        }
        let (_, _, mut guard) = self.weighted_partition(&long)?;
        guard.reindex_clusters(&long);
        let _timer = crate::timing::scope("attract");
        let (graph, knn) = self.embed(contigs);
        let full = self.partition_all(&graph, contigs, None)?;
        let judge = Judge {
            quality: &self.quality,
            contigs,
            bars: self.bars(),
        };
        Ok((graph, knn, combine(&[guard, full], &judge, None)))
    }
}
