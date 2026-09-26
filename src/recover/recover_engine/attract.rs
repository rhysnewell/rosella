use anyhow::Result;
use log::info;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::{Graph, knn::KnnGraph};
use crate::recover::recover_engine::RecoverEngine;

impl RecoverEngine {
    // Short contigs pull contaminants out of good bins on some assemblies and drag them in on
    // others, so each assembly keeps the partition worth more, as the settle judges its passes.
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
        let kept = self.pass_worth(&settled.2, contigs);
        let attracted = self.pass_worth(&full, contigs);
        info!(
            "Marker worth {kept:.0} from contigs of at least {} bp, {attracted:.0} with the {} \
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
