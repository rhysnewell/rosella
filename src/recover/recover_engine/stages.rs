use std::collections::{HashMap, HashSet};

use anyhow::{Result, anyhow};
use log::{debug, warn};

use crate::clustering::clusterer::Partitioning;
use crate::embedding::Graph;
use crate::embedding::features::ContigFeatures;
use crate::embedding::knn::KnnGraph;
use crate::recover::census::Census;
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::dissolve::{PoolSearcher, PoolView, RoundParams};
use crate::refine::finished::Finished;
use crate::refine::owners::owners;
use crate::refine::report_context::Inputs as ReportInputs;
use crate::refine::splitter::{RefineSettings, Refined, Refiner};

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Stage {
    Dissolve,
    Join,
    Shed,
}

pub const SHIPPED_ORDER: &str = "dissolve,join,shed";

const STAGES: [Stage; 3] = [Stage::Dissolve, Stage::Join, Stage::Shed];

impl Stage {
    fn name(self) -> &'static str {
        match self {
            Self::Dissolve => "dissolve",
            Self::Join => "join",
            Self::Shed => "shed",
        }
    }
}

pub fn parse_order(order: &str) -> Result<Vec<Stage>> {
    let held = order
        .split(',')
        .map(str::trim)
        .filter(|name| !name.is_empty())
        .map(|name| {
            STAGES
                .into_iter()
                .find(|stage| stage.name() == name)
                .ok_or_else(|| anyhow!("{name} is not a stage of the refine cycle"))
        })
        .collect::<Result<Vec<_>>>()?;
    for stage in STAGES {
        if !held.contains(&stage) {
            warn!(
                "The stage order leaves out {}, so it never runs.",
                stage.name()
            );
        }
    }
    Ok(held)
}

struct Cycle<'a> {
    refined: &'a mut Refined,
    induced: &'a KnnGraph,
    finished: Finished,
}

impl RecoverEngine {
    pub(super) fn refine_clusters(
        &self,
        partitioning: Partitioning,
        assembly: &Graph,
        induced: &KnnGraph,
        census: &mut Census,
    ) -> Refined {
        let settings = &self.settings;
        let refine = RefineSettings {
            min_bin_size: settings.bars.min_bin_size,
            max_bin_size: settings.lengths.max_bin,
            knn: settings.knn,
            max_retries: settings.cycle.refine_rounds,
            seeds: settings.seeds,
            max_contamination: None,
            partition: settings.partition,
            split: settings.cycle.split,
            peel: settings.cycle.peel,
        };
        let mut refined = Refiner::new(
            self.features(),
            refine,
            partitioning.cluster_map,
            partitioning.outliers,
        )
        .with_assembly(assembly)
        .with_quality(&self.quality)
        .with_cuts(self.reports.cut_report.is_some())
        .run();
        self.census_bins(census, "refine", &refined.bins, &refined.unbinned);

        let mut seen: HashMap<Stage, usize> = HashMap::new();
        let mut cycle = Cycle {
            refined: &mut refined,
            induced,
            finished: Finished::default(),
        };
        for stage in &settings.cycle.order {
            let pass = seen.entry(*stage).or_insert(0);
            self.run_stage(*stage, &mut cycle, census, *pass);
            *pass += 1;
        }

        for contigs in refined.bins.values_mut() {
            contigs.sort_unstable();
        }
        refined
    }

    fn run_stage(&self, stage: Stage, cycle: &mut Cycle<'_>, census: &mut Census, pass: usize) {
        let before = cycle
            .refined
            .cuts
            .is_some()
            .then(|| cycle.refined.bins.clone());
        let ran = match stage {
            Stage::Dissolve if self.settings.cycle.dissolve => {
                cycle.finished = self.dissolve_stage(cycle.refined, cycle.induced);
                true
            }
            Stage::Join if self.settings.cycle.join => {
                self.join_stage(cycle.refined);
                true
            }
            Stage::Shed => {
                self.shed_stage(cycle.refined, cycle.induced, cycle.finished);
                true
            }
            _ => false,
        };
        if !ran {
            return;
        }
        let label = stage_label(stage.name(), pass);
        let refined = &mut *cycle.refined;
        self.census_bins(census, &label, &refined.bins, &refined.unbinned);
        if let (Some(before), Some(log)) = (before, refined.cuts.as_mut()) {
            let after = owners(
                refined
                    .bins
                    .iter()
                    .map(|(label, members)| (*label, members)),
            );
            let lengths = self.lengths();
            log.diff(&label, before.values(), &after, |contig| lengths[contig]);
        }
    }

    fn report_inputs<'a>(
        &'a self,
        features: &'a ContigFeatures<'a>,
        knn: &'a KnnGraph,
    ) -> ReportInputs<'a> {
        ReportInputs {
            features,
            quality: &self.quality,
            knn,
            lengths: self.lengths(),
            names: self.names(),
        }
    }

    fn shed_stage(&self, refined: &mut Refined, knn: &KnnGraph, finished: Finished) {
        if finished.mostly() {
            debug!(
                "Shed skipped: {:.4} of the bins the pool was handed already cleared the bars.",
                finished.share()
            );
            return;
        }
        if let Some(path) = &self.reports.shed_report
            && let Err(error) = crate::refine::shed_report::write(
                path,
                &refined.bins,
                &self.report_inputs(&self.features(), knn),
                &self.quality,
            )
        {
            warn!(
                "Could not write the shed report to {}: {error}",
                path.display()
            );
        }
        let dropped = crate::refine::shed::shed(
            &mut refined.bins,
            &mut refined.unbinned,
            &self.quality,
            self.settings.bars.quality(),
        );
        debug!("Shed {dropped} contigs the bin already held a marker copy for.");
    }

    fn dissolve_stage(&self, refined: &mut Refined, induced: &KnnGraph) -> Finished {
        // Stale by a round, since merge and both eject arms move the bins it was
        // measured on. Recomputing it here was measured and lost bins.
        let settings = crate::refine::dissolve::DissolveSettings {
            bars: self.settings.bars,
            genome_floor: refined.genome_floor,
            rounds: self.settings.cycle.dissolve_rounds,
            passes: self.settings.cycle.dissolve_passes,
            n_neighbours: self.settings.knn.neighbours,
            max_bin_size: self.settings.lengths.max_bin,
        };
        let report = self.reports.pool_report.as_ref().and_then(|path| {
            crate::refine::pool_report::PoolReport::create(path, self.names())
                .map_err(|error| {
                    warn!(
                        "Could not write the pool report to {}: {error}",
                        path.display()
                    )
                })
                .ok()
        });
        let ledger = crate::refine::dissolve::dissolve(
            crate::refine::dissolve::PoolInputs {
                features: &self.features(),
                quality: &self.quality,
                settings,
                oracle: &self.oracle,
                report: report.as_ref(),
            },
            &mut refined.bins,
            &mut refined.unbinned,
            &EngineSearch {
                engine: self,
                induced,
            },
        );
        if let Some(report) = report.as_ref() {
            report.flush();
        }
        debug!("Dissolve pool: {ledger}");
        ledger.finished()
    }

    fn join_stage(&self, refined: &mut Refined) {
        let _timer = crate::timing::scope("join");
        let bars = self.settings.bars;
        let ledger = crate::refine::join::join(
            &self.features(),
            &self.quality,
            &mut refined.bins,
            crate::refine::join::JoinSettings {
                completeness: bars.completeness,
                contamination: bars.contamination,
                max_bin_size: self.settings.lengths.max_bin,
            },
        );
        debug!("Join: {ledger}");
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

// A stage that runs twice needs a row of its own rather than one that overwrites the first.
pub fn stage_label(stage: &str, pass: usize) -> String {
    match pass {
        0 => stage.to_string(),
        _ => format!("{stage}_{}", pass + 1),
    }
}
