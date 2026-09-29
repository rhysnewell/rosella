use anyhow::Result;
use log::info;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::weight::{Contigs, NEIGHBOURS, STEPS, centre, recall};
use crate::embedding::{Graph, knn::KnnGraph};
use crate::quality::{Bars, Scorer};
use crate::recover::partition_report::PartitionReport;
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::select::sorted;

const SAMPLE: usize = 2_000;
// Each pass is a neighbour search and a partition, so a weight that wanders is cut off here.
const SETTLE_PASSES: usize = 8;
// Two plateaus with the same mean can differ in the last bits, which ran a pass twice.
const SAME_WEIGHT: f64 = 1e-9;

impl RecoverEngine {
    pub(super) fn weighted_partition(
        &mut self,
        contigs: &[usize],
    ) -> Result<(Graph, KnnGraph, Partitioning)> {
        let mut report = self
            .partition_report
            .as_ref()
            .map(|_| PartitionReport::new(contigs));
        let (graph, knn) = self.embed(contigs);
        let first = self.partition_all(&graph, contigs, report.as_mut())?;
        let mut partitioned = match self.derived_weight(&self.near_complete(&first, contigs))? {
            None => (graph, knn, first),
            Some(weight) => {
                self.distance.aggregate_weight = Some(weight);
                self.pass(contigs, report.as_mut())?
            }
        };
        if let Some(start) = self.distance.aggregate_weight {
            let mut near_complete = self.near_complete(&partitioned.2, contigs);
            let mut trace = vec![(start, self.pass_worth(&partitioned.2, contigs))];
            let mut best = trace[0];
            while trace.len() < SETTLE_PASSES
                && let Some(weight) = self.derived_weight(&near_complete)?
                && trace
                    .iter()
                    .all(|(used, _)| (used - weight).abs() > SAME_WEIGHT)
            {
                self.distance.aggregate_weight = Some(weight);
                let next = self.pass(contigs, report.as_mut())?;
                near_complete = self.near_complete(&next.2, contigs);
                trace.push((weight, self.pass_worth(&next.2, contigs)));
                if trace[trace.len() - 1].1 > best.1 {
                    best = trace[trace.len() - 1];
                    partitioned = next;
                }
            }
            self.distance.aggregate_weight = Some(best.0);
            info!(
                "Coverage weight and marker worth per pass {trace:.3?}, kept {:.3}.",
                best.0
            );
            self.worth_spread =
                best.1 - trace.iter().map(|(_, worth)| *worth).fold(best.1, f64::min);
        }
        if let (Some(report), Some(path)) = (&report, &self.partition_report) {
            report.write(
                path,
                &self.coverage_table.contig_names,
                &self.coverage_table.contig_lengths,
            )?;
        }
        Ok(partitioned)
    }

    fn pass(
        &self,
        contigs: &[usize],
        report: Option<&mut PartitionReport>,
    ) -> Result<(Graph, KnnGraph, Partitioning)> {
        let mut report = report;
        if let Some(report) = report.as_deref_mut() {
            report.pass("derived");
        }
        let (graph, knn) = self.embed(contigs);
        let partitioning = self.partition_all(&graph, contigs, report)?;
        Ok((graph, knn, partitioning))
    }

    // Squared so one whole genome outweighs its markers split across two bins. Long contigs only,
    // so a pass holding short contigs is judged on the markers the long-only pass had.
    pub(super) fn pass_worth(&self, partitioning: &Partitioning, contigs: &[usize]) -> f64 {
        let lengths = &self.coverage_table.contig_lengths;
        partitioning
            .cluster_map
            .values()
            .map(|members| {
                let long = members
                    .iter()
                    .map(|at| contigs[*at])
                    .filter(|contig| lengths[*contig] >= self.cutoff);
                let worth = self.quality.score(&sorted(long)).score(self.worth);
                worth.max(0.0).powi(2)
            })
            .sum()
    }

    // Near complete bins stand in for labels.
    fn near_complete(&self, partitioning: &Partitioning, contigs: &[usize]) -> Vec<Vec<usize>> {
        let bars = Bars {
            completeness: self.min_completeness,
            contamination: self.contamination_bar,
        };
        partitioning
            .cluster_map
            .values()
            .map(|members| sorted(members.iter().map(|at| contigs[*at])))
            .filter(|members| self.quality.score(members).clears(bars))
            .collect()
    }

    // Every contig split in two must leave halves long enough to be binned in their own right.
    fn derived_weight(&self, near_complete: &[Vec<usize>]) -> Result<Option<f64>> {
        let _timer = crate::timing::scope("weight");
        let lengths = &self.coverage_table.contig_lengths;
        let pool = sorted(
            near_complete
                .iter()
                .flatten()
                .copied()
                .filter(|contig| lengths[*contig] >= 2 * self.cutoff),
        );
        let order =
            crate::seeds::sample_positions(pool.len(), SAMPLE.min(pool.len()), self.seeds.seed);
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
            .map(|contig| self.coverage_table.contig_names[*contig].as_str())
            .collect::<Vec<_>>();
        let [first, second] = crate::kmers::kmer_counting::halves(
            &self.assembly,
            &names,
            self.tnf_table.kmer_size(),
        )?;
        let contigs = Contigs {
            coverage: chosen
                .iter()
                .map(|contig| {
                    crate::embedding::features::row_slice(&self.coverage_table.table, *contig)
                })
                .collect(),
            whole: chosen
                .iter()
                .map(|contig| {
                    crate::embedding::features::row_slice(&self.tnf_table.kmer_table, *contig)
                })
                .collect(),
            first: first.rows().into_iter().map(|row| row.to_vec()).collect(),
            second: second.rows().into_iter().map(|row| row.to_vec()).collect(),
        };
        Ok(recall(
            &contigs,
            self.distance.presence_fraction,
            self.seeds.seed,
        ))
    }
}
