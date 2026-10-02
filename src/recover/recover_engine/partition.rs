use anyhow::Result;
use log::{debug, warn};

use crate::clustering::clusterer::{Partitioning, conserved, find_partitions};
use crate::clustering::graph_partition::Partition;
use crate::embedding::{Graph, KNN_ASSEMBLY, KNN_POOL, features::ContigFeatures, knn::KnnGraph};
use crate::recover::ladder::{Judge, best_per_arm, combine};
use crate::recover::partition_report::PartitionReport;
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::dissolve::{MIN_RESCUE_CONTIGS, PoolView, RoundParams};

impl RecoverEngine {
    pub(super) fn partition_all(
        &self,
        graph: &Graph,
        contigs: &[usize],
        report: Option<&mut PartitionReport>,
    ) -> Result<Partitioning> {
        let mut ladder = Vec::new();
        for step in 0..self.settings.cycle.partition_seeds {
            ladder.extend(self.partition_of(
                graph,
                contigs,
                self.settings.partition,
                true,
                self.settings.seeds.partition.wrapping_add(step as u64),
            )?);
        }
        Ok(self.pick_partition(ladder, contigs, report))
    }

    fn pick_partition(
        &self,
        ladder: Vec<Partitioning>,
        contigs: &[usize],
        mut partitions: Option<&mut PartitionReport>,
    ) -> Partitioning {
        let judge = Judge {
            quality: &self.quality,
            contigs,
            bars: self.settings.bars,
        };
        let report = self.reports.combine_report.as_ref().and_then(|path| {
            crate::recover::combine_report::CombineReport::create(
                path,
                self.names(),
                self.lengths(),
            )
            .map_err(|error| warn!("Could not write the combine report: {error}"))
            .ok()
        });
        if let Some(partitions) = partitions.as_deref_mut() {
            partitions.ladder(&ladder);
        }
        let arms = best_per_arm(ladder, &judge);
        let chosen = combine(&arms, &judge, report.as_ref());
        if let Some(report) = &report {
            report.flush();
        }
        if let Some(partitions) = partitions {
            partitions.chosen(&arms);
            partitions.add("combined", &chosen);
        }
        chosen
    }

    // `contigs` index the contig list as it stands after the initial length filter.
    fn partition_of(
        &self,
        graph: &Graph,
        contigs: &[usize],
        kind: Partition,
        rank_rungs: bool,
        partition_seed: u64,
    ) -> Result<Vec<Partitioning>> {
        find_partitions(
            graph,
            &self.features().contig_lengths(contigs),
            partition_seed,
            kind,
            rank_rungs,
        )
    }

    pub(super) fn embed(&self, features: &ContigFeatures, contigs: &[usize]) -> (Graph, KnnGraph) {
        let built = features.knn_of(
            contigs,
            self.settings.knn.neighbours,
            self.settings.knn.candidates,
            self.settings.seeds.knn,
            KNN_ASSEMBLY,
        );
        let graph = features.graph_from_knn(contigs, &built);
        (graph, built)
    }

    pub(super) fn evaluate_outliers(
        &self,
        partitioning: &mut Partitioning,
        induced: &KnnGraph,
    ) -> Result<()> {
        if partitioning.outliers.len() < MIN_RESCUE_CONTIGS {
            return Ok(());
        }
        let outliers = std::mem::take(&mut partitioning.outliers);
        let neighbours = self.settings.knn.neighbours;
        let (knn, order) = self.pool_neighbours(outliers, neighbours, PoolView::Combined, induced);
        let rescued = self
            .evaluate_subset(
                &knn,
                &order,
                RoundParams {
                    n_neighbours: neighbours,
                    ladder: false,
                },
            )?
            .swap_remove(0);

        debug!("New Partition score {:?}", rescued.score);
        debug!("Number of clusters: {}", rescued.cluster_map.len());
        partitioning.merge(rescued);
        Ok(())
    }

    pub(super) fn pool_neighbours(
        &self,
        order: Vec<usize>,
        n_neighbours: usize,
        view: PoolView,
        induced: &KnnGraph,
    ) -> (KnnGraph, Vec<usize>) {
        if view == PoolView::Combined
            && let Some(built) = induced.induced(&order)
        {
            return (built, order);
        }
        let features = match view {
            PoolView::Combined => self.features(),
            PoolView::Composition => self
                .features()
                .with_distance(self.settings.features.distance.composition_only()),
        };
        let knn = features.knn_of(
            &order,
            n_neighbours,
            self.settings.knn.candidates,
            self.settings.seeds.knn,
            KNN_POOL,
        );
        (knn, order)
    }

    pub(super) fn evaluate_subset(
        &self,
        knn: &KnnGraph,
        ordered_indices: &[usize],
        round: RoundParams,
    ) -> Result<Vec<Partitioning>> {
        let subset_graph = self.features().graph_from_knn(ordered_indices, knn);
        // Label propagation returns one labelling, and the rungs need a ladder to walk.
        let kind = match round.ladder && !self.settings.partition.runs_leiden() {
            true => {
                debug!("Pool rungs need a ladder, so this round runs Leiden");
                Partition::Leiden
            }
            false => self.settings.partition,
        };
        let mut results = self.partition_of(
            &subset_graph,
            ordered_indices,
            kind,
            !round.ladder,
            self.settings.seeds.partition,
        )?;
        if !round.ladder {
            results.truncate(1);
        }
        debug!("Partition score {:?}", results[0].score);

        let expected = ordered_indices.iter().copied().collect();
        for result in results.iter_mut() {
            result.reindex_clusters(ordered_indices);
            conserved(
                result
                    .cluster_map
                    .values()
                    .flatten()
                    .copied()
                    .chain(result.outliers.iter().copied()),
                &expected,
            )?;
        }

        Ok(results)
    }
}
