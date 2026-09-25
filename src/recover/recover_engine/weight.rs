use anyhow::Result;
use log::info;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::weight::{Contigs, NEIGHBOURS, Plateau, derive};
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
            Some(plateau) => {
                self.distance.aggregate_weight = Some(plateau.centre());
                self.pass(contigs, report.as_mut())?
            }
        };
        if let Some(start) = self.distance.aggregate_weight {
            let mut used = vec![start];
            let mut latest = None;
            let mut settled = false;
            while let Some(plateau) = self
                .derived_weight(&self.near_complete(&latest.as_ref().unwrap_or(&partitioned).2))?
            {
                settled = plateau.holds(used[used.len() - 1]);
                let weight = plateau.centre();
                if settled
                    || used.len() == SETTLE_PASSES
                    || used.iter().any(|w| (w - weight).abs() <= SAME_WEIGHT)
                {
                    break;
                }
                self.distance.aggregate_weight = Some(weight);
                latest = Some(self.pass(contigs, report.as_mut())?);
                used.push(weight);
            }
            let kept = if settled { used[used.len() - 1] } else { start };
            if settled && let Some(latest) = latest {
                partitioned = latest;
            }
            self.distance.aggregate_weight = Some(kept);
            info!(
                "Coverage weight per pass {used:.3?}, kept {kept:.3}{}.",
                if settled {
                    " on its own plateau"
                } else {
                    " with no pass on its own plateau"
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
    fn derived_weight(&self, near_complete: &[Vec<usize>]) -> Result<Option<Plateau>> {
        let _timer = crate::timing::scope("weight");
        let lengths = &self.coverage_table.contig_lengths;
        let pool = sorted(
            near_complete
                .iter()
                .flatten()
                .copied()
                .filter(|contig| lengths[*contig] >= 2 * self.min_contig_size),
        );
        let chosen =
            crate::seeds::sample_positions(pool.len(), SAMPLE.min(pool.len()), self.seeds.seed)
                .into_iter()
                .map(|at| pool[at])
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
        let plateau = derive(&contigs, self.distance.presence_fraction, self.seeds.seed);
        info!(
            "Coverage weight {} from {} contigs in {} near complete bins.",
            plateau
                .as_ref()
                .map_or("unchanged".to_string(), |plateau| format!(
                    "{:.3}",
                    plateau.centre()
                )),
            chosen.len(),
            near_complete.len()
        );
        Ok(plateau)
    }
}
