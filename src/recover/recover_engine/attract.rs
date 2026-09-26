use anyhow::Result;
use log::info;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::{Graph, knn::KnnGraph};
use crate::recover::recover_engine::RecoverEngine;

// One complete, clean bin in squared marker worth. A smaller gain cannot add a bin.
const ONE_GENOME: f64 = 100.0 * 100.0;

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
        let (graph, knn, settled) = self.weighted_partition(&contigs[..long])?;
        let kept = self.pass_worth(&settled, contigs);
        let bar = self.worth_spread.max(ONE_GENOME);
        // Past parity the graph is mostly short contigs rather than long ones they attract to.
        let admitted = long + (contigs.len() - long).min(long);
        let share = self.quality.hit_count(&contigs[long..admitted]) as f64
            / self.quality.hit_count(&contigs[..long]).max(1) as f64;
        let reach = kept * ((1.0 + share).powi(2) - 1.0);
        let floor = self.coverage_table.contig_lengths[contigs[admitted - 1]];
        info!(
            "{} shorter contigs down to {floor} bp carry {share:.3} of the long contigs' markers, \
             worth at most {reach:.0} against a bar of {bar:.0}.",
            admitted - long
        );
        self.parked = contigs[long..].to_vec();
        if self.long_only || reach <= bar {
            return Ok((graph, knn, settled));
        }
        // The long-only graph is rebuilt if it wins, so two graphs are never held at once.
        drop((graph, knn));
        let _timer = crate::timing::scope("attract");
        let (graph, knn) = self.embed(&contigs[..admitted]);
        let full = self.partition_all(&graph, &contigs[..admitted], None)?;
        let attracted = self.pass_worth(&full, contigs);
        info!("Marker worth {kept:.0} from long contigs alone, {attracted:.0} with the shorter.");
        if attracted - kept > bar {
            self.parked = contigs[admitted..].to_vec();
            return Ok((graph, knn, full));
        }
        drop((graph, knn));
        let (graph, knn) = self.embed(&contigs[..long]);
        Ok((graph, knn, settled))
    }
}
