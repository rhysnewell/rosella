use anyhow::Result;
use log::info;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::weight::{Contigs, derive};
use crate::embedding::{Graph, knn::KnnGraph};
use crate::quality::{Bars, Scorer};
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::select::sorted;

const SAMPLE: usize = 2_000;

impl RecoverEngine {
    pub(super) fn weighted_partition(
        &mut self,
        contigs: &[usize],
    ) -> Result<(Graph, KnnGraph, Partitioning)> {
        let (graph, knn) = self.embed(contigs);
        let first = self.partition_all(&graph, contigs)?;
        let Some(weight) = self.derived_weight(&first)? else {
            return Ok((graph, knn, first));
        };
        self.distance.aggregate_weight = Some(weight);
        let (graph, knn) = self.embed(contigs);
        let partitioning = self.partition_all(&graph, contigs)?;
        Ok((graph, knn, partitioning))
    }

    // Near complete bins stand in for labels. Every contig split in two must leave halves long
    // enough to be binned in their own right.
    fn derived_weight(&self, partitioning: &Partitioning) -> Result<Option<f64>> {
        let _timer = crate::timing::scope("weight");
        let bars = Bars {
            completeness: self.min_completeness,
            contamination: self.contamination_bar,
        };
        let lengths = &self.coverage_table.contig_lengths;
        let near_complete = partitioning
            .cluster_map
            .values()
            .map(|members| sorted(members.iter().copied()))
            .filter(|members| self.quality.score(members).clears(bars))
            .collect::<Vec<_>>();
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
        if chosen.len() <= crate::embedding::weight::NEIGHBOURS {
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
