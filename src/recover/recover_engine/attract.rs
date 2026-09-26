use std::collections::{HashMap, HashSet};

use anyhow::Result;
use ndarray::Array2;

use crate::clustering::clusterer::{Partitioning, find_partitions};
use crate::coverage::coverage_table::CoverageTable;
use crate::embedding::{KNN_ASSEMBLY, features::ContigFeatures};
use crate::kmers::kmer_counting::KmerFrequencyTable;
use crate::recover::ladder::{Judge, best_per_arm, combine};
use crate::recover::recover_engine::RecoverEngine;

// Contigs under the cutoff give their genome's long contigs neighbours of their own, which pulls
// them out of the bins they contaminate. They sit in the graph only and never enter a bin.
pub struct Attractors {
    coverage: Array2<f64>,
    composition: Array2<f64>,
    lengths: Vec<usize>,
    binnable: Vec<Option<usize>>,
    pub(crate) links: Option<Vec<crate::assembly_graph::Link>>,
    replace: bool,
}

impl Attractors {
    pub fn split(
        coverage: &mut CoverageTable,
        composition: &mut KmerFrequencyTable,
        cutoff: usize,
        replace: bool,
    ) -> Result<Option<Self>> {
        if coverage.contig_lengths.iter().all(|length| *length >= cutoff) {
            return Ok(None);
        }
        let mut next = 0;
        let binnable = coverage
            .contig_lengths
            .iter()
            .map(|length| {
                (*length >= cutoff).then(|| {
                    next += 1;
                    next - 1
                })
            })
            .collect();
        let view = Self {
            coverage: coverage.table.clone(),
            composition: composition.kmer_table.clone(),
            lengths: coverage.contig_lengths.clone(),
            binnable,
            links: None,
            replace,
        };
        let removed = coverage.filter_by_length(cutoff)?;
        composition.filter_by_name(&removed)?;
        Ok(Some(view))
    }

    pub fn binnable(&self, held: Partitioning) -> Partitioning {
        let keep = |members: HashSet<usize>| {
            members
                .into_iter()
                .filter_map(|row| self.binnable[row])
                .collect::<HashSet<_>>()
        };
        let mut clusters = held.cluster_map.into_iter().collect::<Vec<_>>();
        clusters.sort_unstable_by_key(|(label, _)| *label);
        let mut cluster_map = HashMap::new();
        for (_, members) in clusters {
            let kept = keep(members);
            if !kept.is_empty() {
                cluster_map.insert(cluster_map.len(), kept);
            }
        }
        Partitioning {
            cluster_map,
            outliers: keep(held.outliers),
            ..held
        }
    }
}

impl RecoverEngine {
    pub(super) fn attract(&self, main: Partitioning, contigs: &[usize]) -> Result<Partitioning> {
        let Some(view) = &self.attractors else {
            return Ok(main);
        };
        let _timer = crate::timing::scope("attract");
        let features = ContigFeatures::new(&view.coverage, &view.composition, &view.lengths)
            .with_distance(self.distance)
            .with_links(view.links.as_deref(), self.link_weight);
        let every = (0..view.lengths.len()).collect::<Vec<_>>();
        let knn = features.knn_of(
            &every,
            self.n_neighbours,
            self.seeds,
            self.knn_candidates,
            KNN_ASSEMBLY,
        );
        let graph = features.graph_from_knn(&every, &knn);
        let mut ladder = Vec::new();
        for step in 0..self.partition_seeds {
            ladder.extend(
                find_partitions(
                    &graph,
                    &view.lengths,
                    self.ladder_band(),
                    self.seeds.partition + step as u64,
                    self.partition,
                    true,
                    self.leiden,
                )?
                .into_iter()
                .map(|held| view.binnable(held)),
            );
        }
        let judge = Judge {
            quality: &self.quality,
            contigs,
            bars: self.bars(),
        };
        let mut arms = best_per_arm(ladder, &judge);
        if !view.replace {
            arms.push(main);
        }
        Ok(combine(&arms, &judge, None))
    }
}
