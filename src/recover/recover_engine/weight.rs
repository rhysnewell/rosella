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
            while let Some(weight) = self.derived_weight(&self.near_complete(
                &latest.as_ref().unwrap_or(&partitioned).2,
            ))? {
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
                if fixed { " at a fixed point" } else { " with no fixed point" }
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
                .filter(|contig| lengths[*contig] >= 2 * self.cutoff),
        );
        let drawn = if self.weight_blocks {
            pool.len()
        } else {
            SAMPLE.min(pool.len())
        };
        let order = crate::seeds::sample_positions(pool.len(), drawn, self.seeds.seed);
        if order.len() <= NEIGHBOURS {
            info!(
                "{} near complete bins hold {} contigs to judge the coverage weight on, too few to \
                 move it off the sample count.",
                near_complete.len(),
                order.len()
            );
            return Ok(None);
        }
        let blocks = order.len().div_ceil(SAMPLE);
        let size = order.len().div_ceil(blocks);
        let mut curves = Vec::with_capacity(blocks);
        for block in order.chunks(size) {
            let chosen = block.iter().map(|at| pool[*at]).collect::<Vec<_>>();
            if let Some(curve) = self.recall_curve(&chosen)? {
                curves.push((chosen.len(), curve));
            }
        }
        let weight = match curves.as_slice() {
            [] => None,
            [(_, curve)] => Some(centre(curve)),
            _ => {
                let counted = curves.iter().map(|(len, _)| *len).sum::<usize>();
                let mut total = [0.0; STEPS];
                for (len, curve) in &curves {
                    for (sum, value) in total.iter_mut().zip(curve) {
                        *sum += value * *len as f64;
                    }
                }
                Some(centre(&total.map(|sum| sum / counted as f64)))
            }
        };
        info!(
            "Coverage weight {} from {} contigs in {} blocks in {} near complete bins.",
            weight.map_or("unchanged".to_string(), |weight| format!("{weight:.3}")),
            order.len(),
            blocks,
            near_complete.len()
        );
        Ok(weight)
    }

    fn recall_curve(&self, chosen: &[usize]) -> Result<Option<[f64; STEPS]>> {
        if chosen.len() <= NEIGHBOURS {
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
        Ok(recall(&contigs, self.distance.presence_fraction, self.seeds.seed))
    }
}
