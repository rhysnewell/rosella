use std::collections::HashSet;

use anyhow::Result;
use log::debug;
use ndarray::s;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::Graph;
use crate::embedding::knn::KnnGraph;
use crate::quality::Scorer;
use crate::recover::coherence::Coherence;
use crate::recover::fold::Role;
use crate::recover::recover_engine::RecoverEngine;
use crate::recover::strains::{Depths, split_strains};
use crate::refine::dissolve::{PoolSearcher, PoolView, RoundParams};
use crate::refine::finished::Finished;
use crate::refine::splitter::{RefineSettings, Refined, Refiner};

impl RecoverEngine {
    pub(super) fn refine_clusters(
        &self,
        partitioning: Partitioning,
        assembly: &Graph,
        induced: &KnnGraph,
    ) -> (Refined, Finished) {
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
        self.trace_bins("stage3_split", &refined.bins);
        self.split_strains(&mut refined);

        let finished = self.dissolve_stage(&mut refined, induced);
        self.trace_bins("stage4_dissolved", &refined.bins);
        self.shed_stage(&mut refined, finished);
        self.trace_bins("stage5_shed", &refined.bins);

        for contigs in refined.bins.values_mut() {
            contigs.sort_unstable();
        }
        (refined, finished)
    }

    // The vote puts offered bins back after the refine cycle, so its winners have met neither rule.
    pub(super) fn finish_offered(&self, refined: &mut Refined, finished: Finished) {
        self.split_strains(refined);
        self.shed_stage(refined, finished);
        for contigs in refined.bins.values_mut() {
            contigs.sort_unstable();
        }
    }

    fn split_strains(&self, refined: &mut Refined) {
        let lengths = &self.tables.coverage.contig_lengths;
        let depths = self.tables.coverage.table.slice(s![.., ..;2]);
        let hosts = refined
            .bins
            .values()
            .filter(|members| Role::of(self.quality.score(members), false) == Role::Host)
            .map(Vec::as_slice)
            .collect::<Vec<_>>();
        let Some(coherence) = Coherence::of(&hosts, depths, lengths) else {
            return;
        };
        let depths = Depths {
            depths,
            lengths,
            read_lengths: self.tables.read_lengths.as_deref(),
        };
        let splits = split_strains(
            &mut refined.bins,
            &self.quality,
            |members| self.quality.carriers(members),
            &depths,
            &coherence,
            &self.features(),
        );
        debug!("Strain splits: {splits}");
    }

    fn shed_stage(&self, refined: &mut Refined, finished: Finished) {
        let mut bars = self.settings.bars.quality();
        if finished.mostly() {
            bars.completeness = 0.0;
        }
        let dropped = crate::refine::shed::shed(
            &mut refined.bins,
            &mut refined.unbinned,
            &self.quality,
            bars,
        );
        debug!("Shed: {dropped} contigs (share {:.4})", finished.share());
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
