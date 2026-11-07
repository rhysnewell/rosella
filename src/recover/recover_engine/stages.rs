use std::collections::HashSet;

use anyhow::Result;
use log::debug;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::Graph;
use crate::embedding::knn::KnnGraph;
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::dissolve::{PoolSearcher, PoolView, RoundParams};
use crate::refine::finished::Finished;
use crate::refine::splitter::{RefineSettings, Refined, Refiner};

impl RecoverEngine {
    pub(super) fn refine_clusters(
        &self,
        partitioning: Partitioning,
        assembly: &Graph,
        induced: &KnnGraph,
    ) -> Refined {
        let settings = &self.settings;
        let refine = RefineSettings {
            piece: settings.bars.piece(),
            max_bin_size: settings.lengths.max_bin,
            knn: settings.knn,
            max_retries: settings.cycle.refine_rounds,
            seeds: settings.seeds,
            max_contamination: settings.split_contamination,
            partition: settings.partition,
        };
        let mut refined = Refiner::new(
            self.features(),
            refine,
            partitioning.cluster_map,
            partitioning.outliers,
            assembly,
            &self.quality,
        )
        .with_contamination(
            self.given
                .as_ref()
                .map(|given| given.contamination.clone())
                .unwrap_or_default(),
        )
        .run();

        let finished = self.dissolve_stage(&mut refined, induced);
        self.shed_stage(&mut refined, finished);

        for contigs in refined.bins.values_mut() {
            contigs.sort_unstable();
        }
        refined
    }

    fn shed_stage(&self, refined: &mut Refined, finished: Finished) {
        if finished.mostly() {
            debug!(
                "Shed: skipped ({:.4} of pool bins already clear the bars)",
                finished.share()
            );
            return;
        }
        let dropped = crate::refine::shed::shed(
            &mut refined.bins,
            &mut refined.unbinned,
            &self.quality,
            self.settings.bars.quality(),
        );
        debug!("Shed: {dropped} contigs");
    }

    fn dissolve_stage(&self, refined: &mut Refined, induced: &KnnGraph) -> Finished {
        let genomes = self.quality.genomes();
        debug!("Finished share read against {genomes} genomes");
        let settings = crate::refine::dissolve::DissolveSettings {
            bars: self.settings.bars,
            rounds: self.settings.cycle.dissolve_rounds,
            passes: self.settings.cycle.dissolve_passes,
            n_neighbours: self.settings.knn.neighbours,
            max_bin_size: self.settings.lengths.max_bin,
            genomes,
        };
        let fair = crate::markers::Fair(&self.quality);
        let ledger = crate::refine::dissolve::dissolve(
            crate::refine::dissolve::PoolInputs {
                features: &self.features(),
                quality: &fair,
                settings,
            },
            &mut refined.bins,
            &mut refined.unbinned,
            &EngineSearch {
                engine: self,
                induced,
            },
        );
        debug!("Dissolve pool: {ledger}");
        ledger.finished()
    }
}

struct EngineSearch<'a> {
    engine: &'a RecoverEngine,
    induced: &'a KnnGraph,
}

impl PoolSearcher for EngineSearch<'_> {
    fn neighbours(
        &self,
        pool: &HashSet<usize>,
        n_neighbours: usize,
        view: PoolView,
    ) -> Result<(KnnGraph, Vec<usize>)> {
        Ok(self.engine.pool_neighbours(
            crate::refine::ranking::sorted(pool.iter().copied()),
            n_neighbours,
            view,
            self.induced,
        ))
    }

    fn partition(
        &self,
        knn: &KnnGraph,
        order: &[usize],
        round: RoundParams,
    ) -> Result<Vec<Partitioning>> {
        self.engine.evaluate_subset(knn, order, round)
    }
}
