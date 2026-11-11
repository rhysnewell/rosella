use std::collections::HashSet;

use anyhow::Result;
use log::debug;

use crate::clustering::clusterer::{Partitioning, conserved, find_partitions};
use crate::clustering::graph_partition::Partition;
use crate::embedding::{Graph, KNN_ASSEMBLY, KNN_POOL, features::ContigFeatures, knn::KnnGraph};
use crate::progress::{self, Stage};
use crate::recover::ladder::{Judge, combine, top_rungs_per_arm};
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::dissolve::{MIN_RESCUE_CONTIGS, PoolView, RoundParams};

impl RecoverEngine {
    pub(super) fn partition_all(&self, graph: &Graph, contigs: &[usize]) -> Result<Partitioning> {
        let seeds = (0..self.settings.cycle.partition_seeds)
            .map(|step| self.settings.seeds.partition.wrapping_add(step as u64))
            .collect::<Vec<_>>();
        let ladder = self.partition_of(graph, contigs, self.settings.partition, true, &seeds)?;
        Ok(self.pick_partition(ladder, contigs))
    }

    fn pick_partition(&self, ladder: Vec<Partitioning>, contigs: &[usize]) -> Partitioning {
        let features = self.features();
        let counted = crate::markers::Fair(&self.quality);
        if self.quality.learns_checkm() {
            self.learn_checkm(&ladder, contigs);
        }
        let judge = Judge {
            quality: &counted,
            contigs,
            lengths: self.lengths(),
            bars: self.settings.bars,
            depth: &features,
            markers: &self.quality,
        };
        let traced = self
            .stages
            .as_ref()
            .map(|report| (report, report.next_pass()));
        if let Some((report, pass)) = traced {
            for (at, held) in ladder.iter().enumerate() {
                let name = format!("pass{pass}_rung{at}_{}_s{}", held.arm.name(), held.seed);
                report.partitioning(&name, held, contigs, &self.tables.coverage.contig_names);
            }
        }
        let mut arms = top_rungs_per_arm(ladder, &judge);
        let own = arms.len();
        let offered =
            !self.candidates.is_empty() && std::env::var_os("ROSELLA_LATE_OFFERS").is_none();
        // Another binner's bin can leave out a genome's own paralog, so the contest has to count
        // multi-copy families once, as the refine cycle does.
        if offered && let Some(best) = arms.first() {
            self.learn_multicopy(best, contigs);
        }
        if offered {
            arms.extend(self.candidate_arms(contigs));
        }
        let voters = (0..arms.len())
            .map(|arm| arm.saturating_sub(own - 1))
            .collect::<Vec<_>>();
        let combined = combine(&arms, &voters, &judge);
        if let Some((report, pass)) = traced {
            let names = &self.tables.coverage.contig_names;
            report.combined(&format!("pass{pass}"), &arms, &combined, contigs, names);
        }
        if offered {
            self.log_sources(&arms[own..], own, &combined.sources);
        }
        combined.partitioning
    }

    // Rungs repeat most bins, so each distinct bin reads once.
    fn learn_checkm(&self, ladder: &[Partitioning], contigs: &[usize]) {
        let distinct = ladder
            .iter()
            .flat_map(|partitioning| partitioning.cluster_map.values())
            .collect::<HashSet<_>>();
        let bins = distinct
            .into_iter()
            .map(|positions| positions.iter().map(|at| contigs[*at]).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        self.quality
            .learn_checkm(bins.iter().map(Vec::as_slice), self.settings.bars.quality());
    }

    fn log_sources(&self, offered: &[Partitioning], own: usize, sources: &[(Vec<usize>, bool)]) {
        let mut taken = vec![(0, 0); offered.len() + 1];
        for (arms, whole) in sources {
            let mut seen = arms
                .iter()
                .map(|arm| arm.saturating_sub(own - 1).min(offered.len()))
                .collect::<Vec<_>>();
            seen.dedup();
            for set in seen {
                match whole {
                    true => taken[set].0 += 1,
                    false => taken[set].1 += 1,
                }
            }
        }
        let names =
            std::iter::once("rosella").chain(self.candidates.iter().map(|set| set.name.as_str()));
        let offered =
            std::iter::once(None).chain(offered.iter().map(|held| Some(held.cluster_map.len())));
        for ((name, offered), (whole, trimmed)) in names.zip(offered).zip(taken) {
            let offered = offered.map_or_else(String::new, |count| format!(" of {count} offered"));
            debug!("Ensemble: {name} gave {whole} bins whole and {trimmed} trimmed{offered}");
        }
    }

    pub(super) fn candidate_arms(&self, contigs: &[usize]) -> Vec<Partitioning> {
        if self.candidates.is_empty() {
            return Vec::new();
        }
        let at = crate::rows::positions(contigs.iter().copied());
        self.candidates
            .iter()
            .map(|set| Partitioning {
                cluster_map: set
                    .bins
                    .iter()
                    .filter_map(|(_, members)| {
                        let mut held = members
                            .iter()
                            .filter_map(|contig| at.get(contig).copied())
                            .collect::<Vec<_>>();
                        held.sort_unstable();
                        (!held.is_empty()).then_some(held)
                    })
                    .enumerate()
                    .collect(),
                outliers: Vec::new(),
                score: None,
                arm: self.settings.partition,
                seed: self.settings.seeds.partition,
            })
            .collect()
    }

    fn partition_of(
        &self,
        graph: &Graph,
        contigs: &[usize],
        kind: Partition,
        rank_rungs: bool,
        seeds: &[u64],
    ) -> Result<Vec<Partitioning>> {
        find_partitions(
            graph,
            &self.features().contig_lengths(contigs),
            seeds,
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
        let progress = progress::hold(progress::spinning(Stage::RescuingUnbinned));
        progress.set_message(format!("{} contigs", outliers.len()));
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

        debug!("Pool partition score: {:?}", rescued.score);
        debug!("Pool clusters: {}", rescued.cluster_map.len());
        partitioning.merge(rescued);
        progress::done(&progress);
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
                debug!("Pool round: Leiden (rungs need a ladder)");
                Partition::Leiden
            }
            false => self.settings.partition,
        };
        let mut results = self.partition_of(
            &subset_graph,
            ordered_indices,
            kind,
            !round.ladder,
            &[self.settings.seeds.partition],
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
