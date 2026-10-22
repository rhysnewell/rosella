use std::collections::{BTreeMap, HashMap};
use std::ops::Range;

use anyhow::Result;
use log::debug;

use super::{ONE_GENOME, best_bins, merge};
use crate::embedding::knn::{KnnGraph, nearest_exact};
use crate::recover::floor_walk::{Bar, Reach};
use crate::recover::recover_engine::RecoverEngine;
use crate::recover::recover_engine::weight::Worth;

impl RecoverEngine {
    // Every stop reads the marker contigs alone, searched exactly, so the search over every
    // contig covers only the bands attached.
    pub(super) fn plan(
        &mut self,
        order: &[usize],
        spans: &[Range<usize>],
        long_graph: &KnnGraph,
        bin_of: &HashMap<usize, usize>,
        bins: &BTreeMap<usize, Vec<usize>>,
        worth: Worth,
    ) -> Result<Vec<f32>> {
        let long = (0..long_graph.indices.nrows()).collect::<Vec<_>>();
        let mut reach = Reach::new(
            worth.kept,
            worth.spread.max(ONE_GENOME),
            self.quality.hit_count(&long),
        );
        let mut bar = Bar::default();
        let mut bars = Vec::new();
        let mut ceiling = self.settings.lengths.cutoff;
        for span in spans {
            let band = &order[span.clone()];
            let lengths = &self.tables.coverage.contig_lengths;
            let low = lengths[band[band.len() - 1]];
            let sequence = band.iter().map(|contig| lengths[*contig]).sum::<usize>();
            let annotation = self.annotator.annotate(low..ceiling)?;
            self.quality
                .fill(annotation, &self.tables.coverage.contig_names, band);
            ceiling = low;
            let added = reach.add(self.quality.hit_count(band));
            if added <= reach.bar() {
                debug!(
                    "Attach stops at {low} bp: adds {added:.0} against a bar of {:.0}",
                    reach.bar()
                );
                break;
            }
            self.load(band.iter().max().map_or(0, |last| last + 1))?;
            let (joins, marked) = self.marked_joins(order, span, long_graph, bin_of);
            let share_of = joins
                .iter()
                .map(|(contig, _, share)| (*contig, *share))
                .collect::<HashMap<_, _>>();
            let pairs = joins
                .iter()
                .map(|(contig, bin, _)| (*contig, *bin))
                .collect::<Vec<_>>();
            let seen = self.marker_evidence(bins, &pairs);
            let Some(needs) = bar.add(
                seen.iter()
                    .map(|seen| (share_of[&seen.contig], seen.in_place, seen.complete)),
            ) else {
                debug!("Attach stops at {low} bp: no share reads under half foreign");
                break;
            };
            let above = joins.iter().filter(|(_, _, share)| *share > needs).count();
            let joining = sequence as f64 * above as f64 / marked.max(1) as f64;
            if joining < self.settings.bars.min_bin_size as f64 {
                debug!(
                    "Attach stops at {low} bp: joins {joining:.0} bp against a {} bp bin",
                    self.settings.bars.min_bin_size
                );
                break;
            }
            debug!("Attach from {low} bp: neighbour share over {needs:.2}");
            bars.push(needs);
        }
        Ok(bars)
    }

    // A contig with a marker is coding and prokaryotic, so it joins at least as often as its band
    // does.
    fn marked_joins(
        &self,
        order: &[usize],
        span: &Range<usize>,
        long_graph: &KnnGraph,
        bin_of: &HashMap<usize, usize>,
    ) -> (Vec<(usize, usize, f32)>, usize) {
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
        let knn = merge(&nearest, &among, marked.len(), |at| order[at]);
        let joins = contigs
            .iter()
            .zip(best_bins(&knn, long_graph.indices.nrows(), bin_of))
            .filter_map(|(contig, best)| {
                let (bin, share) = best?;
                (share > 0.5).then_some((*contig, bin, share))
            })
            .collect();
        (joins, marked.len())
    }
}
