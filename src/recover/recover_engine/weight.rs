use anyhow::Result;
use log::info;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::weight::{Contigs, NEIGHBOURS, derive};
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
        let mut partitioned = match self.derived_weight(&self.near_complete(&first))? {
            None => (graph, knn, first),
            Some(weight) => {
                self.distance.aggregate_weight = Some(weight);
                self.pass(contigs, report.as_mut())?
            }
        };
        if let Some(start) = self.distance.aggregate_weight {
            let mut used = vec![start];
            let mut latest = None;
            let mut fixed = false;
            while let Some(weight) = self
                .derived_weight(&self.near_complete(&latest.as_ref().unwrap_or(&partitioned).2))?
            {
                fixed = (weight - used[used.len() - 1]).abs() <= SAME_WEIGHT;
                if fixed
                    || used.len() == SETTLE_PASSES
                    || used.iter().any(|w| (w - weight).abs() <= SAME_WEIGHT)
                {
                    break;
                }
                self.distance.aggregate_weight = Some(weight);
                latest = Some(self.pass(contigs, report.as_mut())?);
                used.push(weight);
            }
            let kept = if fixed { used[used.len() - 1] } else { start };
            if fixed && let Some(latest) = latest {
                partitioned = latest;
            }
            self.distance.aggregate_weight = Some(kept);
            info!(
                "Coverage weight per pass {used:.3?}, kept {kept:.3}{}.",
                if fixed {
                    " at a fixed point"
                } else {
                    " with no fixed point"
                }
            );
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

    // Near complete bins stand in for labels.
    fn near_complete(&self, partitioning: &Partitioning) -> Vec<Vec<usize>> {
        let bars = Bars {
            completeness: self.min_completeness,
            contamination: self.contamination_bar,
        };
        partitioning
            .cluster_map
            .values()
            .map(|members| sorted(members.iter().copied()))
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
                .filter(|contig| lengths[*contig] >= 2 * self.min_contig_size),
        );
        let kept = crate::seeds::consistent_sample(&pool, SAMPLE, self.seeds.seed);
        let chosen = crate::seeds::sample_positions(kept.len(), kept.len(), self.seeds.seed)
            .into_iter()
            .map(|at| kept[at])
            .collect::<Vec<_>>();
        if chosen.len() <= NEIGHBOURS {
            info!(
                "{} near complete bins hold {} contigs to judge the coverage weight on, too few to \
                 move it off the sample count.",
                near_complete.len(),
                chosen.len()
            );
            return Ok(None);
        }
        let names = chosen
            .iter()
            .map(|contig| self.coverage_table.contig_names[*contig].as_str())
            .collect::<Vec<_>>();
        let [first, second] = crate::kmers::kmer_counting::halves(
            &self.assembly,
            &names,
            &self.tnf_table.kmer_sizes(),
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
        let weight = derive(&contigs, self.distance.presence_fraction, self.seeds.seed);
        info!(
            "Coverage weight {} from {} contigs in {} near complete bins.",
            weight.map_or("unchanged".to_string(), |weight| format!("{weight:.3}")),
            chosen.len(),
            near_complete.len()
        );
        Ok(weight)
    }
}
