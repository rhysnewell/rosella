use std::collections::{HashMap, HashSet};
use std::ops::Range;

use anyhow::Result;
use log::info;

use super::{ONE_GENOME, best_bins, merge};
use crate::embedding::knn::{KnnGraph, nearest_exact};
use crate::recover::floor_walk::{Foreign, Reach};
use crate::recover::recover_engine::RecoverEngine;

pub(super) struct Plan {
    pub(super) taken: Vec<bool>,
    pub(super) searched: usize,
}

impl Plan {
    pub(super) fn given() -> Self {
        Self {
            taken: vec![true],
            searched: 1,
        }
    }
}

impl RecoverEngine {
    // Every stop reads the marker contigs alone, searched exactly, so the search over every
    // contig covers only the bands some bin may still take.
    pub(super) fn plan(
        &mut self,
        order: &[usize],
        spans: &[Range<usize>],
        long_graph: &KnnGraph,
        bin_of: &HashMap<usize, usize>,
        bins: &HashMap<usize, HashSet<usize>>,
    ) -> Result<Plan> {
        let long = (0..long_graph.indices.nrows()).collect::<Vec<_>>();
        let mut reach = Reach::new(
            self.kept,
            self.worth_spread.max(ONE_GENOME),
            self.quality.hit_count(&long),
        );
        let mut foreign = Foreign::default();
        let mut taken = Vec::new();
        let mut ceiling = self.cutoff;
        for span in spans {
            let band = &order[span.clone()];
            let lengths = &self.coverage_table.contig_lengths;
            let low = lengths[band[band.len() - 1]];
            let sequence = band.iter().map(|contig| lengths[*contig]).sum::<usize>();
            let annotation = self.annotator.annotate(low..ceiling)?;
            self.quality
                .fill(annotation, &self.coverage_table.contig_names, band);
            ceiling = low;
            let added = reach.add(self.quality.hit_count(band));
            if added <= reach.bar() {
                info!(
                    "Contigs from {low} bp could add {added:.0} against a bar of {:.0}, so \
                     attach searches no deeper.",
                    reach.bar()
                );
                break;
            }
            let (joins, marked) = self.marked_joins(order, span, long_graph, bin_of);
            let joining = sequence as f64 * joins.len() as f64 / marked.max(1) as f64;
            if joining < self.min_bin_size as f64 {
                info!(
                    "Contigs from {low} bp would join about {joining:.0} bp going by their \
                     marker contigs, short of a {} bp bin, so attach searches no deeper.",
                    self.min_bin_size
                );
                break;
            }
            if taken.last() == Some(&false) {
                taken.push(false);
                continue;
            }
            let seen = self.marker_evidence(bins, &joins);
            let admitted = foreign.admits(
                seen.iter().filter(|seen| seen.in_place).count(),
                seen.iter().map(|seen| seen.complete).sum(),
            );
            if !admitted {
                info!(
                    "Contigs from {low} bp bring the in-place foreign share to {:.2}, so from \
                     there each bin walks on alone while its own markers vouch for its joins.",
                    foreign.share().unwrap_or(f64::NAN),
                );
            }
            taken.push(admitted);
        }
        Ok(Plan {
            searched: taken.len(),
            taken,
        })
    }

    // A contig with a marker is coding and prokaryotic, so it joins at least as often as its band
    // does.
    fn marked_joins(
        &self,
        order: &[usize],
        span: &Range<usize>,
        long_graph: &KnnGraph,
        bin_of: &HashMap<usize, usize>,
    ) -> (Vec<(usize, usize)>, usize) {
        let marked = span
            .clone()
            .filter(|at| self.quality.hit_count(&[order[*at]]) > 0)
            .collect::<Vec<_>>();
        if marked.is_empty() {
            return (Vec::new(), 0);
        }
        let contigs = marked.iter().map(|at| order[*at]).collect::<Vec<_>>();
        let nearest = self.nearest_long(long_graph, &contigs);
        let among = {
            let _timer = crate::timing::scope("attach");
            let prepared = self.features().prepared(&order[..span.end]);
            nearest_exact(span.end, &marked, long_graph.indices.ncols(), &prepared)
        };
        let knn = merge(&nearest, |row| row, &among, marked.len(), |at| order[at]);
        let joins = contigs
            .iter()
            .zip(best_bins(&knn, long_graph.indices.nrows(), bin_of))
            .filter_map(|(contig, best)| {
                let (bin, share) = best?;
                (share > 0.5).then_some((*contig, bin))
            })
            .collect();
        (joins, marked.len())
    }
}
