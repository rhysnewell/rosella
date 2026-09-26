use anyhow::Result;
use log::info;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::{Graph, knn::KnnGraph};
use crate::recover::recover_engine::RecoverEngine;

impl RecoverEngine {
    // Short contigs pull contaminants out of good bins on some assemblies and drag them in on
    // others, so each assembly keeps whichever partition its markers prefer.
    pub(super) fn partitioned(
        &mut self,
        contigs: &[usize],
    ) -> Result<(Graph, KnnGraph, Partitioning)> {
        let lengths = &self.coverage_table.contig_lengths;
        let long = contigs.partition_point(|contig| lengths[*contig] >= self.cutoff);
        if long == contigs.len() {
            return self.weighted_partition(contigs);
        }
        let settled = self.weighted_partition(&contigs[..long])?;
        let _timer = crate::timing::scope("attract");
        let (graph, knn) = self.embed(contigs);
        let full = self.partition_all(&graph, contigs, None)?;
        let kept = self.near_complete(&settled.2, contigs).len();
        let attracted = self.near_complete(&full, contigs).len();
        info!(
            "{kept} near complete bins from contigs of at least {} bp, {attracted} with the {} \
             shorter ones.",
            self.cutoff,
            contigs.len() - long
        );
        if attracted > kept {
            return Ok((graph, knn, full));
        }
        self.parked = contigs[long..].to_vec();
        Ok(settled)
    }
}
