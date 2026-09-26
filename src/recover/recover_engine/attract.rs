use anyhow::Result;
use log::info;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::{Graph, knn::KnnGraph};
use crate::recover::recover_engine::RecoverEngine;

impl RecoverEngine {
    // Short contigs help some assemblies and hurt others. A gain no bigger than the worth the
    // weight alone moved across this run's passes is noise, so it keeps the long-only partition.
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
             shorter ones, against a spread of {:.0} across the passes.",
            self.cutoff,
            contigs.len() - long,
            self.worth_spread
        );
        if attracted - kept > self.worth_spread {
            return Ok((graph, knn, full));
        }
        self.parked = contigs[long..].to_vec();
        Ok(settled)
    }
}
