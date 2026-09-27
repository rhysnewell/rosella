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
        let share = self.quality.hit_count(&contigs[long..]) as f64
            / self.quality.hit_count(&contigs[..long]).max(1) as f64;
        let reach = kept * ((1.0 + share).powi(2) - 1.0);
        info!(
            "{} shorter contigs carry {share:.3} of the long contigs' markers, worth at most \
             {reach:.0} against a bar of {bar:.0}.",
            contigs.len() - long
        );
        self.parked = contigs[long..].to_vec();
        if reach <= bar {
            return Ok((graph, knn, settled));
        }
        let glue = self.glue(&knn, &settled, contigs, long);
        // The long-only graph is rebuilt if it wins, so two graphs are never held at once.
        drop((graph, knn));
        let _timer = crate::timing::scope("attract");
        let admitted = contigs[..long]
            .iter()
            .chain(&glue.chosen)
            .copied()
            .collect::<Vec<_>>();
        let (graph, knn) = self.embed(&admitted);
        let mut full = self.partition_all(&graph, &admitted, None)?;
        let attracted = self.pass_worth(&full, &admitted);
        info!(
            "Marker worth {kept:.0} from long contigs alone, {attracted:.0} with {} shorter.",
            glue.chosen.len()
        );
        if attracted - kept > bar {
            full.reindex_clusters(&admitted);
            let chosen = glue
                .chosen
                .iter()
                .copied()
                .collect::<std::collections::HashSet<_>>();
            self.parked.retain(|contig| !chosen.contains(contig));
            if self.glue_attach {
                self.nearest = Some(glue.nearest);
            }
            let rows = self.n_contigs;
            return Ok((
                crate::embedding::lifted(&graph, &admitted, rows),
                knn.lifted(&admitted, rows),
                full,
            ));
        }
        drop((graph, knn));
        let (graph, knn) = self.embed(&contigs[..long]);
        Ok((graph, knn, settled))
    }
}
