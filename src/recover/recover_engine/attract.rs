use anyhow::Result;
use log::info;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::{
    Graph,
    knn::{KnnGraph, nearest_in},
};
use crate::recover::recover_engine::RecoverEngine;

// One complete, clean bin in squared marker worth. A smaller gain cannot add a bin.
const ONE_GENOME: f64 = 100.0 * 100.0;

pub(super) struct Nearest {
    pub(super) first: usize,
    pub(super) knn: KnnGraph,
}

impl RecoverEngine {
    // Short contigs pull long ones out of foreign bins, but most short contigs a partition places
    // are foreign themselves. So they shape the partition, leave it, and return only by attach.
    // A gain no bigger than the worth the weight alone moved across this run's passes is noise.
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
        // Rebuilt from its kNN afterwards, so two graphs are never held at once.
        drop(graph);
        let nearest = self.nearest_long(&knn, contigs, long);
        let full = {
            let _timer = crate::timing::scope("attract");
            let start = ndarray::concatenate(
                ndarray::Axis(0),
                &[knn.indices.view(), nearest.knn.indices.view()],
            )?;
            let (graph, _) = self.embed_from(contigs, &start);
            self.partition_all(&graph, contigs, None)?
        };
        let attracted = self.pass_worth(&full, contigs);
        info!(
            "Long contigs' marker worth {kept:.0} alone, {attracted:.0} with every shorter contig \
             in the graph."
        );
        if self.glue_attach {
            self.nearest = Some(nearest);
        }
        let graph = self.features().graph_from_knn(&contigs[..long], &knn);
        match attracted - kept > bar {
            true => Ok((graph, knn, long_only(full, long))),
            false => Ok((graph, knn, settled)),
        }
    }

    fn nearest_long(&self, knn: &KnnGraph, contigs: &[usize], long: usize) -> Nearest {
        let _timer = crate::timing::scope("nearest");
        let prepared = self.features().prepared(contigs);
        Nearest {
            first: contigs[long],
            knn: nearest_in(
                knn,
                contigs.len() - long,
                knn.indices.ncols(),
                self.knn_candidates,
                self.seeds.knn,
                |short, base| prepared.distance(long + short, base),
            ),
        }
    }

    fn embed_from(&self, contigs: &[usize], start: &ndarray::Array2<u32>) -> (Graph, KnnGraph) {
        let features = self.features();
        let built = features.knn_from(
            contigs,
            start,
            self.n_neighbours,
            self.seeds,
            self.knn_candidates,
        );
        let graph = features.graph_from_knn(contigs, &built);
        (graph, built)
    }
}

fn long_only(mut full: Partitioning, long: usize) -> Partitioning {
    full.cluster_map.retain(|_, members| {
        members.retain(|at| *at < long);
        !members.is_empty()
    });
    full.outliers.retain(|at| *at < long);
    full
}
