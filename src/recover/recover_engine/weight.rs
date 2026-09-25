use anyhow::Result;
use log::info;

use crate::clustering::clusterer::Partitioning;
use crate::coverage::scatter::{self, Source};
use crate::embedding::weight::{Contigs, NEIGHBOURS, derive};
use crate::embedding::{Graph, knn::KnnGraph};
use crate::quality::{Bars, Scorer};
use crate::recover::partition_report::PartitionReport;
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::select::sorted;
use rand::{SeedableRng, rngs::StdRng};

const SAMPLE: usize = 2_000;

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
        let near_complete = self.near_complete(&first);
        let rescaled = self.fit_scatter(&near_complete);
        let weight = self.derived_weight(&near_complete)?;
        if let Some(weight) = weight {
            self.distance.aggregate_weight = Some(weight);
        }
        let partitioned = match weight.is_some() || rescaled {
            false => (graph, knn, first),
            true => {
                if let Some(report) = report.as_mut() {
                    report.pass("derived");
                }
                let (graph, knn) = self.embed(contigs);
                let partitioning = self.partition_all(&graph, contigs, report.as_mut())?;
                (graph, knn, partitioning)
            }
        };
        if let (Some(report), Some(path)) = (&report, &self.partition_report) {
            report.write(
                path,
                &self.coverage_table.contig_names,
                &self.coverage_table.contig_lengths,
            )?;
        }
        Ok(partitioned)
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

    fn fit_scatter(&mut self, near_complete: &[Vec<usize>]) -> bool {
        let Some(source) = self.coverage_variance.clone() else {
            return false;
        };
        let _timer = crate::timing::scope("scatter");
        let lengths = &self.coverage_table.contig_lengths;
        let models = match source {
            Source::Bins => scatter::fit(&self.coverage_table.table, lengths, near_complete),
            Source::Neighbours => self.neighbour_scatter(),
            Source::Given(models) => models.into_iter().map(Some).collect(),
        };
        for (sample, model) in models.iter().enumerate() {
            match model {
                Some(model) => info!(
                    "Sample {sample} depth scatter: sampling {:.3e}, bias {:.3e}.",
                    model.sampling, model.bias
                ),
                None => info!("Sample {sample} has too few contigs to fit its depth scatter."),
            }
        }
        if models.iter().all(Option::is_none) {
            return false;
        }
        scatter::apply(
            &mut self.coverage_table.table,
            lengths,
            &models,
            self.distance.variance_floor,
        );
        self.distance.variance_floor = 0.0;
        true
    }

    fn neighbour_scatter(&self) -> Vec<Option<scatter::Scatter>> {
        let lengths = &self.coverage_table.contig_lengths;
        let long = (0..self.n_contigs)
            .filter(|contig| lengths[*contig] >= 2 * self.min_contig_size)
            .collect::<Vec<_>>();
        let pool =
            crate::seeds::sample_positions(long.len(), SAMPLE.min(long.len()), self.seeds.seed)
                .into_iter()
                .map(|at| long[at])
                .collect::<Vec<_>>();
        if pool.len() <= NEIGHBOURS {
            return vec![None; self.coverage_table.table.ncols() / 2];
        }
        let whole = pool
            .iter()
            .map(|contig| {
                crate::embedding::features::row_slice(&self.tnf_table.kmer_table, *contig)
            })
            .collect::<Vec<_>>();
        let neighbours = crate::embedding::weight::noise::nearest(&whole, NEIGHBOURS);
        scatter::fit_neighbours(
            &self.coverage_table.table,
            lengths,
            &pool,
            &neighbours,
            NEIGHBOURS,
            &mut StdRng::seed_from_u64(self.seeds.seed),
        )
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
        let weight = derive(&contigs, self.distance, self.seeds.seed);
        info!(
            "Coverage weight {} from {} contigs in {} near complete bins.",
            weight.map_or("unchanged".to_string(), |weight| format!("{weight:.3}")),
            chosen.len(),
            near_complete.len()
        );
        Ok(weight)
    }
}
