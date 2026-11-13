use std::collections::{BTreeMap, HashMap};

use rayon::prelude::*;

use super::Proposal;
use crate::embedding::knn::KnnGraph;
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::owners::credit;

pub(super) type Vote = (usize, Option<(usize, f32)>, f32);

// The share attach reads counts short neighbours that cannot vote, so the report also gives the
// best bin's share of the long neighbours alone.
pub(super) fn long_shares(
    knn: &KnnGraph,
    first: usize,
    bin_of: &HashMap<usize, usize>,
) -> Vec<f32> {
    let (sigmas, rhos) = crate::embedding::fuzzy::scales(knn.dists.view());
    (0..knn.indices.nrows())
        .into_par_iter()
        .map(|row| {
            let mut mass = Vec::<(usize, f32)>::new();
            let mut long = 0.0;
            for (neighbour, distance) in knn.indices.row(row).iter().zip(knn.dists.row(row)) {
                if *neighbour == u32::MAX {
                    break;
                }
                let neighbour = *neighbour as usize;
                if neighbour >= first {
                    continue;
                }
                let weight = crate::embedding::fuzzy::membership(*distance, rhos[row], sigmas[row]);
                long += weight;
                if let Some(bin) = bin_of.get(&neighbour) {
                    credit(&mut mass, *bin, weight);
                }
            }
            mass.into_iter().map(|(_, held)| held).fold(0.0, f32::max)
                / f32::max(long, f32::MIN_POSITIVE)
        })
        .collect()
}

pub(super) fn with_long(proposals: Vec<Proposal>, long: Vec<f32>) -> Vec<Vote> {
    proposals
        .into_iter()
        .zip(long)
        .map(|((contig, best), long)| (contig, best, long))
        .collect()
}

impl RecoverEngine {
    pub(super) fn report_votes(
        &self,
        name: &str,
        votes: &[Vote],
        bins: &BTreeMap<usize, Vec<usize>>,
    ) {
        let Some(report) = &self.stages else {
            return;
        };
        let names = &self.tables.coverage.contig_names;
        let lengths = &self.tables.coverage.contig_lengths;
        let rows = votes.iter().map(|(contig, best, long)| {
            let (first, share) = best.map_or(("", 0.0), |(bin, share)| {
                (names[bins[&bin][0]].as_str(), share)
            });
            format!(
                "{}\t{}\t{first}\t{share:.4}\t{long:.4}",
                names[*contig], lengths[*contig]
            )
        });
        report.table(name, "contig\tlength\tbin_first\tshare\tlong_share", rows);
    }
}
