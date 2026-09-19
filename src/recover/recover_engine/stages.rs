use anyhow::{Result, anyhow};
use log::{debug, warn};

use crate::embedding::knn::KnnGraph;
use crate::recover::census::Census;
use crate::recover::recover_engine::{MIN_RESCUE_CONTIGS, RecoverEngine};
use crate::refine::rung::Bars;
use crate::refine::splitter::Refiner;

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Stage {
    Dissolve,
    Join,
    Recruit,
}

pub const SHIPPED_ORDER: &str = "dissolve,join,recruit";

pub fn parse_order(order: &str) -> Result<Vec<Stage>> {
    order
        .split(',')
        .map(str::trim)
        .filter(|name| !name.is_empty())
        .map(|name| match name {
            "dissolve" => Ok(Stage::Dissolve),
            "join" => Ok(Stage::Join),
            "recruit" => Ok(Stage::Recruit),
            other => Err(anyhow!("{other} is not a stage of the refine cycle")),
        })
        .collect()
}

impl RecoverEngine {
    pub(super) fn run_stage(
        &self,
        stage: Stage,
        refiner: &mut Refiner,
        induced: &KnnGraph,
        bars: Bars,
        census: &mut Census,
        pass: usize,
    ) {
        match stage {
            Stage::Dissolve if self.dissolve => {
                self.dissolve_stage(refiner, induced, bars, census, pass)
            }
            Stage::Join if self.join => self.join_stage(refiner, bars, census, pass),
            Stage::Recruit if self.recruit => {
                self.recruit_stage(refiner, induced, bars, census, pass)
            }
            _ => {}
        }
    }

    fn dissolve_stage(
        &self,
        refiner: &mut Refiner,
        induced: &KnnGraph,
        bars: Bars,
        census: &mut Census,
        pass: usize,
    ) {
        // Stale by a round, since merge and both eject arms move the bins it was
        // measured on. Recomputing it here was measured and lost bins.
        let settings = crate::refine::dissolve::DissolveSettings {
            bars,
            hold: self.dissolve_hold,
            genome_floor: refiner.genome_floor,
            min_contigs: MIN_RESCUE_CONTIGS,
            rounds: self.dissolve_rounds,
            passes: self.dissolve_passes,
            n_neighbours: self.n_neighbours,
            max_bin_size: self.max_bin_size,
            reembed: self.dissolve_reembed,
        };
        let report = self.pool_report.as_ref().and_then(|path| {
            crate::refine::pool_report::PoolReport::create(path, &self.coverage_table.contig_names)
                .map_err(|error| warn!("No pool report at {}: {error}", path.display()))
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
            &mut refiner.bins,
            &mut refiner.unbinned,
            crate::refine::dissolve::PoolSearch::new(
                |pool, n_neighbours, view| self.pool_neighbours(pool, n_neighbours, view, induced),
                |knn, order, round| self.evaluate_subset(knn, order, round),
            ),
        );
        if let Some(report) = report.as_ref() {
            report.flush();
        }
        debug!("Dissolve pool: {ledger}");
        self.census_bins(
            census,
            stage_label("dissolve", pass),
            &refiner.bins,
            &refiner.unbinned,
        );
    }

    fn join_stage(&self, refiner: &mut Refiner, bars: Bars, census: &mut Census, pass: usize) {
        let _timer = crate::timing::scope("join");
        let ledger = crate::refine::join::join(
            &self.features(),
            &self.quality,
            &mut refiner.bins,
            crate::refine::join::JoinSettings {
                completeness: bars.completeness,
                contamination: self.contamination_bar,
                max_bin_size: self.max_bin_size,
            },
        );
        debug!("Join: {ledger}");
        self.census_bins(
            census,
            stage_label("join", pass),
            &refiner.bins,
            &refiner.unbinned,
        );
    }

    fn recruit_stage(
        &self,
        refiner: &mut Refiner,
        induced: &KnnGraph,
        bars: Bars,
        census: &mut Census,
        pass: usize,
    ) {
        let _timer = crate::timing::scope("recruit");
        let ledger = crate::refine::recruit::recruit(
            &self.features(),
            &self.quality,
            induced,
            &mut refiner.bins,
            crate::refine::recruit::RecruitSettings {
                floor: bars.completeness * self.recruit_floor,
                confidence: self.recruit_confidence,
                completeness: bars.completeness,
                contamination: self.contamination_bar,
                max_bin_size: self.max_bin_size,
                passes: crate::tuning::JOIN_PASSES,
            },
        );
        debug!("Recruit: {ledger}");
        self.census_bins(
            census,
            stage_label("recruit", pass),
            &refiner.bins,
            &refiner.unbinned,
        );
    }
}

/// The census keys on a static name, so a stage that runs twice needs a second one rather
/// than a row that overwrites the first.
fn stage_label(stage: &'static str, pass: usize) -> &'static str {
    match (stage, pass) {
        ("dissolve", 0) => "dissolve",
        ("dissolve", _) => "dissolve_2",
        ("join", 0) => "join",
        ("join", _) => "join_2",
        ("recruit", 0) => "recruit",
        _ => "recruit_2",
    }
}
