use std::collections::BTreeSet;

use anyhow::Result;
use log::{info, warn};
use rayon::prelude::*;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::metrics::DistanceSettings;
use crate::embedding::weight::{
    Contigs, NEIGHBOURS, STEPS, UNITS, centre, recall_each, unit_at, unit_centre,
};
use crate::embedding::{Graph, knn::KnnGraph};
use crate::quality::{Bars, Scorer, edge_spread};
use crate::recover::partition_report::{PULL, PartitionReport, SAMPLE_PASS, keep_pass, pass_file};
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::ranking::sorted;

const SAMPLE: usize = 2_000;
// Each pass is a neighbour search and a partition, so a weight that wanders is cut off here.
const SETTLE_PASSES: usize = 8;
// Two plateaus with the same mean can differ in the last bits, which ran a pass twice.
const SAME_WEIGHT: f64 = 1e-9;

pub(super) struct Worth {
    pub kept: f64,
    pub spread: f64,
}

pub(super) struct Settled {
    pub graph: Graph,
    pub knn: KnnGraph,
    pub partitioning: Partitioning,
    pub weight: Option<f64>,
    pub unit: Option<f64>,
    pub parked: Vec<usize>,
    pub worth: Worth,
}

type Pass = (Graph, KnnGraph, Partitioning);

struct Walk {
    weight: f64,
    spread: f64,
    kept: usize,
    passes: usize,
}

fn pass_name(at: usize) -> String {
    format!("pass{}", at + 1)
}

impl RecoverEngine {
    // Gold says short contigs in the partition move foreign long contigs home but add no bins,
    // and the short contigs a partition places are mostly foreign. So they wait for attach.
    pub(super) fn partitioned(&self, contigs: &[usize]) -> Result<Settled> {
        let lengths = self.lengths();
        let cutoff = self.settings.lengths.cutoff;
        let long = contigs.partition_point(|contig| lengths[*contig] >= cutoff);
        let mut settled = self.weighted_partition(&contigs[..long])?;
        if long < contigs.len() {
            settled.worth.kept = self.pass_worth(&settled.partitioning, contigs);
            settled.parked = contigs[long..].to_vec();
        }
        Ok(settled)
    }

    fn weighted_partition(&self, contigs: &[usize]) -> Result<Settled> {
        let mut report = self
            .reports
            .partition_report
            .as_ref()
            .map(|_| PartitionReport::new(contigs));
        let (graph, knn) = self.embed(&self.features(), contigs);
        let first = self.partition_all(&graph, contigs, report.as_mut(), SAMPLE_PASS)?;
        let near_complete = self.near_complete(&first, contigs);
        let mut distance = self.settings.features.distance.clone();
        if self.settings.derive_unit {
            distance.variance_unit = self.derived_unit(&near_complete)?;
        }
        let start = self.derived_weight(&near_complete, &distance)?;
        let mut partitioned = match start {
            None => (graph, knn, first),
            Some(weight) => self.pass(contigs, &distance, weight, report.as_mut(), 0)?,
        };
        let mut worth = Worth {
            kept: 0.0,
            spread: 0.0,
        };
        let mut passes = vec![SAMPLE_PASS.to_string()];
        let (weight, kept) = match start {
            None => (None, SAMPLE_PASS.to_string()),
            Some(start) => {
                let settled =
                    self.settle(contigs, &distance, start, &mut partitioned, &mut report)?;
                worth.spread = settled.spread;
                passes.extend((0..settled.passes).map(pass_name));
                (Some(settled.weight), pass_name(settled.kept))
            }
        };
        if let Some(path) = &self.reports.combine_report
            && let Err(error) = keep_pass(path, &kept, &passes)
        {
            warn!("Could not keep the combine report of the kept pass: {error}");
        }
        if let (Some(report), Some(path)) = (&mut report, &self.reports.partition_report) {
            report.keep(&kept);
            report.write(path, self.names(), self.lengths())?;
            keep_pass(&pass_file(path, PULL), &kept, &passes)?;
        }
        let (graph, knn, partitioning) = partitioned;
        Ok(Settled {
            graph,
            knn,
            partitioning,
            weight,
            unit: distance.variance_unit,
            parked: Vec::new(),
            worth,
        })
    }

    fn settle(
        &self,
        contigs: &[usize],
        distance: &DistanceSettings,
        start: f64,
        partitioned: &mut Pass,
        report: &mut Option<PartitionReport>,
    ) -> Result<Walk> {
        let mut near_complete = self.near_complete(&partitioned.2, contigs);
        let mut trace = vec![(start, self.pass_worth(&partitioned.2, contigs))];
        let mut best = trace[0];
        let mut kept = 0;
        while trace.len() < SETTLE_PASSES
            && let Some(weight) = self.derived_weight(&near_complete, distance)?
            && trace
                .iter()
                .all(|(used, _)| (used - weight).abs() > SAME_WEIGHT)
        {
            let next = self.pass(contigs, distance, weight, report.as_mut(), trace.len())?;
            let worth = self.pass_worth(&next.2, contigs);
            let spread = self.spread(&partitioned.2, &next.2, contigs);
            // The one standard error rule (Breiman et al. 1984): an edge inside its own noise
            // keeps the pass already held.
            let clears = worth - best.1 > spread;
            info!(
                "Coverage weight {weight:.3} moves marker worth by {:+.0} against a spread of \
                 {spread:.0}, {}.",
                worth - best.1,
                if clears { "kept" } else { "refused" }
            );
            trace.push((weight, worth));
            // A refused pass says nothing about where the weight belongs, so the walk ends there.
            if !clears {
                break;
            }
            best = (weight, worth);
            kept = trace.len() - 1;
            near_complete = self.near_complete(&next.2, contigs);
            *partitioned = next;
        }
        info!(
            "Coverage weight and marker worth per pass {trace:.3?}, kept {:.3}.",
            best.0
        );
        Ok(Walk {
            weight: best.0,
            spread: best.1 - trace.iter().map(|(_, worth)| *worth).fold(best.1, f64::min),
            kept,
            passes: trace.len(),
        })
    }

    fn pass(
        &self,
        contigs: &[usize],
        distance: &DistanceSettings,
        weight: f64,
        report: Option<&mut PartitionReport>,
        at: usize,
    ) -> Result<Pass> {
        let pass = pass_name(at);
        let mut report = report;
        if let Some(report) = report.as_deref_mut() {
            report.pass(&pass);
        }
        let mut settings = self.settings.features.clone();
        settings.distance = DistanceSettings {
            aggregate_weight: Some(weight),
            ..distance.clone()
        };
        let (graph, knn) = self.embed(&self.tables.features(&settings), contigs);
        let partitioning = self.partition_all(&graph, contigs, report, &pass)?;
        Ok((graph, knn, partitioning))
    }

    // Squared so one whole genome outweighs its markers split across two bins. Long contigs only,
    // so a pass holding short contigs is judged on the markers the long-only pass had.
    pub(super) fn pass_worth(&self, partitioning: &Partitioning, contigs: &[usize]) -> f64 {
        partitioning
            .cluster_map
            .par_iter()
            .map(|(_, members)| {
                let long = self.long(members.iter().map(|at| contigs[*at]));
                let worth = self.quality.score(&long).score(self.settings.bars.worth);
                worth.max(0.0).powi(2)
            })
            .collect::<Vec<_>>()
            .into_iter()
            .sum()
    }

    // A bin both passes hold reads the same under any resample of the markers, so it cancels.
    fn spread(&self, kept: &Partitioning, next: &Partitioning, contigs: &[usize]) -> f64 {
        let bins = |partitioning: &Partitioning| {
            partitioning
                .cluster_map
                .values()
                .map(|members| self.long(members.iter().map(|at| contigs[*at])))
                .filter(|members| !members.is_empty())
                .collect::<BTreeSet<_>>()
        };
        let (before, after) = (bins(kept), bins(next));
        let points = |bin: &Vec<usize>| self.quality.points(bin, self.settings.bars.worth);
        edge_spread(
            after
                .difference(&before)
                .map(|bin| (1.0, points(bin)))
                .chain(before.difference(&after).map(|bin| (-1.0, points(bin)))),
        )
    }

    fn long(&self, members: impl Iterator<Item = usize>) -> Vec<usize> {
        let lengths = &self.tables.coverage.contig_lengths;
        sorted(members.filter(|contig| lengths[*contig] >= self.settings.lengths.cutoff))
    }

    // Near complete bins stand in for labels.
    fn near_complete(&self, partitioning: &Partitioning, contigs: &[usize]) -> Vec<Vec<usize>> {
        let bars = Bars {
            completeness: self.settings.bars.completeness,
            contamination: self.settings.bars.contamination,
        };
        partitioning
            .cluster_map
            .par_iter()
            .map(|(_, members)| sorted(members.iter().map(|at| contigs[*at])))
            .filter(|members| self.quality.score(members).clears(bars))
            .collect()
    }

    // Every contig split in two must leave halves long enough to be binned in their own right.
    fn judged(&self, near_complete: &[Vec<usize>]) -> Option<Vec<usize>> {
        let lengths = &self.tables.coverage.contig_lengths;
        let pool = sorted(
            near_complete
                .iter()
                .flatten()
                .copied()
                .filter(|contig| lengths[*contig] >= 2 * self.settings.lengths.cutoff),
        );
        let order = crate::seeds::sample_positions(
            pool.len(),
            SAMPLE.min(pool.len()),
            self.settings.seeds.seed,
        );
        if order.len() <= NEIGHBOURS {
            info!(
                "{} near complete bins hold {} contigs to judge the coverage distance on, too few \
                 to move it off the sample count.",
                near_complete.len(),
                order.len()
            );
            return None;
        }
        Some(order.iter().map(|at| pool[*at]).collect())
    }

    fn derived_weight(
        &self,
        near_complete: &[Vec<usize>],
        distance: &DistanceSettings,
    ) -> Result<Option<f64>> {
        let _timer = crate::timing::scope("weight");
        if distance.evidence.is_some() {
            return Ok(None);
        }
        let Some(chosen) = self.judged(near_complete) else {
            return Ok(None);
        };
        let weight = self
            .recall_curves(&chosen, std::slice::from_ref(distance))?
            .map(|curves| centre(&curves[0]));
        info!(
            "Coverage weight {} from {} contigs in {} near complete bins.",
            weight.map_or("unchanged".to_string(), |weight| format!("{weight:.3}")),
            chosen.len(),
            near_complete.len()
        );
        Ok(weight)
    }

    fn derived_unit(&self, near_complete: &[Vec<usize>]) -> Result<Option<f64>> {
        let _timer = crate::timing::scope("weight");
        let Some(chosen) = self.judged(near_complete) else {
            return Ok(None);
        };
        let units = (0..UNITS)
            .map(|at| DistanceSettings {
                variance_unit: Some(unit_at(at)),
                ..self.settings.features.distance.clone()
            })
            .collect::<Vec<_>>();
        let Some(curves) = self.recall_curves(&chosen, &units)? else {
            return Ok(None);
        };
        let unit = unit_centre(&curves);
        let best = curves
            .iter()
            .enumerate()
            .map(|(at, curve)| {
                let recall = curve.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                format!("{:.3}:{recall:.4}", unit_at(at))
            })
            .collect::<Vec<_>>();
        info!(
            "Variance unit {unit:.3} from {} contigs, best recall per unit [{}].",
            chosen.len(),
            best.join(", ")
        );
        Ok(Some(unit))
    }

    fn recall_curves(
        &self,
        chosen: &[usize],
        settings: &[DistanceSettings],
    ) -> Result<Option<Vec<[f64; STEPS]>>> {
        let names = chosen
            .iter()
            .map(|contig| self.tables.coverage.contig_names[*contig].as_str())
            .collect::<Vec<_>>();
        let [first, second] = self.tables.halves_of(&names)?;
        let contigs = Contigs {
            lengths: chosen
                .iter()
                .map(|contig| self.lengths()[*contig])
                .collect(),
            coverage: chosen
                .iter()
                .map(|contig| {
                    crate::embedding::features::row_slice(&self.tables.coverage.table, *contig)
                })
                .collect(),
            whole: chosen
                .iter()
                .map(|contig| {
                    crate::embedding::features::row_slice(&self.tables.tnf.kmer_table, *contig)
                })
                .collect(),
            first: (0..chosen.len())
                .map(|at| crate::embedding::features::row_slice(&first, at))
                .collect(),
            second: (0..chosen.len())
                .map(|at| crate::embedding::features::row_slice(&second, at))
                .collect(),
        };
        Ok(recall_each(&contigs, settings, self.settings.seeds.seed))
    }
}
