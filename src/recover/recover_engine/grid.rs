use anyhow::Result;
use log::info;

use crate::embedding::fuzzy::{Kernel, reshaped_graph};
use crate::embedding::knn::KnnGraph;
use crate::recover::ladder::{Judge, combine};
use crate::recover::partition_report::PartitionReport;
use crate::recover::recover_engine::RecoverEngine;

// Multiples of the median neighbour distance, so every assembly spans the same contrast.
const BANDWIDTHS: [f32; 5] = [0.125, 0.25, 0.5, 1.0, 2.0];
const CUTOFFS: [usize; 3] = [50, 80, 100];

fn percentiles(knn: &KnnGraph) -> Option<impl Fn(usize) -> f32> {
    let mut distances = knn
        .dists
        .iter()
        .copied()
        .filter(|distance| distance.is_finite())
        .collect::<Vec<_>>();
    if distances.is_empty() {
        return None;
    }
    distances.sort_unstable_by(f32::total_cmp);
    Some(move |percent: usize| distances[(distances.len() - 1) * percent / 100])
}

impl RecoverEngine {
    pub(super) fn grid(
        &self,
        knn: &KnnGraph,
        contigs: &[usize],
        report: &mut PartitionReport,
        kept: &str,
    ) -> Result<()> {
        let _timer = crate::timing::scope("partition_grid");
        let Some(at) = percentiles(knn) else {
            return Ok(());
        };
        let judge = Judge {
            quality: &self.quality,
            contigs,
            bars: self.settings.bars,
        };
        let bars = crate::quality::Bars {
            completeness: judge.bars.completeness,
            contamination: judge.bars.tier(),
        };
        let kernels = std::iter::once(("fuzzy".to_string(), Kernel::PerContig)).chain(
            BANDWIDTHS
                .iter()
                .map(|scale| (format!("b{scale}"), Kernel::Global(scale * at(50)))),
        );
        let mut cells = Vec::new();
        for (name, kernel) in kernels {
            for percent in CUTOFFS {
                let pass = format!("{kept}_grid_{name}_p{percent}");
                let graph = reshaped_graph(contigs.len(), knn, kernel, at(percent));
                report.pass(&pass);
                let ladder = self.ladder_of(&graph, contigs)?;
                let combined = self.pick_partition(ladder, contigs, Some(report), &pass);
                let over = combined
                    .cluster_map
                    .values()
                    .filter(|members| judge.score(members).clears(bars))
                    .count();
                info!(
                    "Grid cell {name} p{percent}: {} edges, {} combined bins, {over} over the bar.",
                    graph.nnz(),
                    combined.cluster_map.len()
                );
                cells.push(combined);
            }
        }
        report.pass(kept);
        report.add("grid_combined", &combine(&cells, &judge, None));
        Ok(())
    }
}
