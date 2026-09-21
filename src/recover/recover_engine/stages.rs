use anyhow::{Result, anyhow};
use log::{debug, warn};

use crate::embedding::knn::KnnGraph;
use crate::recover::census::Census;
use crate::recover::recover_engine::{MIN_RESCUE_CONTIGS, RecoverEngine};
use crate::refine::finished::Finished;
use crate::refine::rung::Bars;
use crate::refine::splitter::Refiner;

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Stage {
    Dissolve,
    Join,
    Recruit,
    Audit,
    Shed,
}

pub const SHIPPED_ORDER: &str = "dissolve,join,recruit,audit,shed";

const NAMES: [(&str, Stage); 5] = [
    ("dissolve", Stage::Dissolve),
    ("join", Stage::Join),
    ("recruit", Stage::Recruit),
    ("audit", Stage::Audit),
    ("shed", Stage::Shed),
];

pub fn parse_order(order: &str) -> Result<Vec<Stage>> {
    let held = order
        .split(',')
        .map(str::trim)
        .filter(|name| !name.is_empty())
        .map(|name| {
            NAMES
                .iter()
                .find(|(known, _)| *known == name)
                .map(|(_, stage)| *stage)
                .ok_or_else(|| anyhow!("{name} is not a stage of the refine cycle"))
        })
        .collect::<Result<Vec<_>>>()?;
    for (name, stage) in NAMES {
        if !held.contains(&stage) {
            warn!("The stage order leaves out {name}, so it never runs.");
        }
    }
    Ok(held)
}

impl RecoverEngine {
    pub(super) fn run_stage(
        &self,
        stage: Stage,
        refiner: &mut Refiner,
        induced: &KnnGraph,
        bars: Bars,
        census: &mut Census,
        finished: &mut Finished,
        pass: usize,
    ) {
        match stage {
            Stage::Dissolve if self.dissolve => {
                self.dissolve_stage(refiner, induced, bars, census, finished, pass)
            }
            Stage::Join if self.join => self.join_stage(refiner, bars, census, pass),
            Stage::Recruit if self.recruit => {
                self.recruit_stage(refiner, induced, bars, census, pass)
            }
            Stage::Audit => self.audit_stage(refiner, induced, census, pass),
            Stage::Shed => self.shed_stage(refiner, induced, bars, census, *finished, pass),
            _ => {}
        }
    }

    fn audit_stage(&self, refiner: &mut Refiner, knn: &KnnGraph, census: &mut Census, pass: usize) {
        if let Some(path) = &self.audit_report
            && let Err(error) = crate::refine::audit_report::write(
                path,
                &refiner.bins,
                &self.features(),
                &self.quality,
                knn,
                &self.coverage_table.contig_lengths,
                &self.coverage_table.contig_names,
            )
        {
            warn!("No audit report at {}: {error}", path.display());
        }
        let evicted = crate::refine::audit::audit(
            &mut refiner.bins,
            &mut refiner.unbinned,
            knn,
            &self.coverage_table.contig_lengths,
        );
        debug!("Audit unbinned {evicted} contigs.");
        self.census_bins(
            census,
            &stage_label("audit", pass),
            &refiner.bins,
            &refiner.unbinned,
        );
    }

    fn shed_stage(
        &self,
        refiner: &mut Refiner,
        knn: &KnnGraph,
        bars: Bars,
        census: &mut Census,
        finished: Finished,
        pass: usize,
    ) {
        if self.finished_gate && finished.mostly() {
            debug!(
                "Shed skipped: {:.4} of the bins the pool was handed already cleared the bars.",
                finished.share()
            );
            self.census_bins(
                census,
                &stage_label("shed", pass),
                &refiner.bins,
                &refiner.unbinned,
            );
            return;
        }
        if let Some(path) = &self.shed_report
            && let Err(error) = crate::refine::shed_report::write(
                path,
                &refiner.bins,
                crate::refine::shed_report::Inputs {
                    features: &self.features(),
                    markers: &self.quality,
                    knn,
                    lengths: &self.coverage_table.contig_lengths,
                    names: &self.coverage_table.contig_names,
                    spacings: self.shed_length_multiple,
                },
            )
        {
            warn!("No shed report at {}: {error}", path.display());
        }
        let features = self.features();
        let top = refiner
            .genome_floor
            .unwrap_or(bars.min_bin_size)
            .max(bars.min_bin_size);
        let rung = bars.at(top, 0);
        let settled = |members: &[usize]| {
            self.dissolve_hold
                .holds(&features, &self.quality, members, bars, rung)
        };
        let dropped = crate::refine::shed::shed(
            &mut refiner.bins,
            &mut refiner.unbinned,
            &self.quality,
            crate::quality::Bars {
                completeness: bars.completeness,
                contamination: self.contamination_bar,
            },
            self.shed_length_multiple,
            &settled,
            self.shed_split.then_some(crate::refine::shed::Split {
                features: &features,
                min_bin_size: self.min_bin_size,
                seed: self.seeds.partition,
            }),
        );
        debug!("Shed {dropped} contigs the bin already held a marker copy for.");
        self.census_bins(
            census,
            &stage_label("shed", pass),
            &refiner.bins,
            &refiner.unbinned,
        );
    }

    fn dissolve_stage(
        &self,
        refiner: &mut Refiner,
        induced: &KnnGraph,
        bars: Bars,
        census: &mut Census,
        finished: &mut Finished,
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
            rung_walk: self.dissolve_rung_walk,
            finished_gate: self.finished_gate,
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
        *finished = ledger.finished();
        debug!("Dissolve pool: {ledger}");
        self.census_bins(
            census,
            &stage_label("dissolve", pass),
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
            &stage_label("join", pass),
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
            &stage_label("recruit", pass),
            &refiner.bins,
            &refiner.unbinned,
        );
    }
}

/// A stage that runs twice needs a row of its own rather than one that overwrites the first.
pub fn stage_label(stage: &str, pass: usize) -> String {
    match pass {
        0 => stage.to_string(),
        _ => format!("{stage}_{}", pass + 1),
    }
}
