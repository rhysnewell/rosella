use std::collections::BTreeSet;

use anyhow::Result;
use log::info;
use rayon::prelude::*;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::weight::{Contigs, NEIGHBOURS, STEPS, centre, recall};
use crate::embedding::{Graph, knn::KnnGraph};
use crate::quality::{Bars, Scorer, edge_spread};
use crate::recover::partition_report::PartitionReport;
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
    pub parked: Vec<usize>,
    pub worth: Worth,
}

type Pass = (Graph, KnnGraph, Partitioning);

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
        let first = self.partition_all(&graph, contigs, report.as_mut())?;
        let start = self.derived_weight(&self.near_complete(&first, contigs))?;
        let mut partitioned = match start {
            None => (graph, knn, first),
            Some(weight) => self.pass(contigs, weight, report.as_mut())?,
        };
        let mut worth = Worth {
            kept: 0.0,
            spread: 0.0,
        };
        let weight = match start {
            None => None,
            Some(start) => {
                let (weight, spread) =
                    self.settle(contigs, start, &mut partitioned, &mut report)?;
                worth.spread = spread;
                Some(weight)
            }
        };
        if let (Some(report), Some(path)) = (&report, &self.reports.partition_report) {
            report.write(path, self.names(), self.lengths())?;
        }
        let (graph, knn, partitioning) = partitioned;
        Ok(Settled {
            graph,
            knn,
            partitioning,
            weight,
            parked: Vec::new(),
            worth,
        })
    }

    fn settle(
        &self,
        contigs: &[usize],
        start: f64,
        partitioned: &mut Pass,
        report: &mut Option<PartitionReport>,
    ) -> Result<(f64, f64)> {
        let mut near_complete = self.near_complete(&partitioned.2, contigs);
        let mut trace = vec![(start, self.pass_worth(&partitioned.2, contigs))];
        let mut best = trace[0];
        while trace.len() < SETTLE_PASSES
            && let Some(weight) = self.derived_weight(&near_complete)?
            && trace
                .iter()
                .all(|(used, _)| (used - weight).abs() > SAME_WEIGHT)
        {
            let next = self.pass(contigs, weight, report.as_mut())?;
            near_complete = self.near_complete(&next.2, contigs);
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
            if clears {
                best = (weight, worth);
                *partitioned = next;
            }
        }
        info!(
            "Coverage weight and marker worth per pass {trace:.3?}, kept {:.3}.",
            best.0
        );
        let spread = best.1 - trace.iter().map(|(_, worth)| *worth).fold(best.1, f64::min);
        Ok((best.0, spread))
    }

    fn pass(
        &self,
        contigs: &[usize],
        weight: f64,
        report: Option<&mut PartitionReport>,
    ) -> Result<Pass> {
        let mut report = report;
        if let Some(report) = report.as_deref_mut() {
            report.pass("derived");
        }
        let mut settings = self.settings.features;
        settings.distance.aggregate_weight = Some(weight);
        let (graph, knn) = self.embed(&self.tables.features(settings), contigs);
        let partitioning = self.partition_all(&graph, contigs, report)?;
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
    fn derived_weight(&self, near_complete: &[Vec<usize>]) -> Result<Option<f64>> {
        let _timer = crate::timing::scope("weight");
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
                "{} near complete bins hold {} contigs to judge the coverage weight on, too few to \
                 move it off the sample count.",
                near_complete.len(),
                order.len()
            );
            return Ok(None);
        }
        let chosen = order.iter().map(|at| pool[*at]).collect::<Vec<_>>();
        let weight = self.recall_curve(&chosen)?.map(|curve| centre(&curve));
        info!(
            "Coverage weight {} from {} contigs in {} near complete bins.",
            weight.map_or("unchanged".to_string(), |weight| format!("{weight:.3}")),
            order.len(),
            near_complete.len()
        );
        Ok(weight)
    }

    fn recall_curve(&self, chosen: &[usize]) -> Result<Option<[f64; STEPS]>> {
        let names = chosen
            .iter()
            .map(|contig| self.tables.coverage.contig_names[*contig].as_str())
            .collect::<Vec<_>>();
        let [first, second] = self
            .tables
            .halves
            .composition(&names, self.tables.tnf.kmer_size())?;
        let contigs = Contigs {
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
        Ok(recall(
            &contigs,
            self.settings.features.distance.presence_fraction,
            self.settings.seeds.seed,
        ))
    }
}
