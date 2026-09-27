use std::collections::{HashMap, HashSet};

use log::info;
use rayon::prelude::*;

use crate::embedding::knn::KnnGraph;
use crate::quality::Scorer;
use crate::recover::recover_engine::RecoverEngine;
use crate::recover::recover_engine::attract::Nearest;

struct Evidence {
    bin: usize,
    repeats: bool,
    complete: f64,
}

impl RecoverEngine {
    // Only long contigs vouch for a bin, over a neighbourhood that counts every contig. Letting
    // short ones vouch fills sink bins through chains of short contigs.
    pub(super) fn attach(
        &self,
        bins: &mut HashMap<usize, HashSet<usize>>,
        unbinned: &mut HashSet<usize>,
        parked: &[usize],
        nearest: &Nearest,
    ) {
        let _timer = crate::timing::scope("attach");
        let knn = self.neighbourhoods(parked, nearest);
        let bin_of = bins
            .iter()
            .flat_map(|(bin, members)| members.iter().map(move |contig| (*contig, *bin)))
            .collect::<HashMap<_, _>>();
        let (sigmas, rhos) = crate::embedding::fuzzy::scales(knn.dists.view(), knn.indices.ncols());
        let joins = parked
            .par_iter()
            .enumerate()
            .filter_map(|(row, contig)| {
                let mut mass = HashMap::<usize, f32>::new();
                let mut total = 0.0;
                for (neighbour, distance) in knn.indices.row(row).iter().zip(knn.dists.row(row)) {
                    if *neighbour == u32::MAX {
                        break;
                    }
                    let weight =
                        crate::embedding::fuzzy::membership(*distance, rhos[row], sigmas[row]);
                    total += weight;
                    let neighbour = *neighbour as usize;
                    if neighbour < nearest.first
                        && let Some(bin) = bin_of.get(&neighbour)
                    {
                        *mass.entry(*bin).or_default() += weight;
                    }
                }
                let (bin, held) = mass
                    .into_iter()
                    .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)))?;
                (held > total / 2.0).then_some((*contig, bin))
            })
            .collect::<Vec<_>>();
        let refused = self.refused(bins, &joins);
        let mut taken = 0;
        for (contig, bin) in joins.iter().filter(|(_, bin)| !refused.contains(bin)) {
            unbinned.remove(contig);
            bins.entry(*bin).or_default().insert(*contig);
            taken += 1;
        }
        info!(
            "{} of {} parked short contigs sit in one bin's neighbourhood. {} bins refuse theirs \
             on marker evidence and {taken} join.",
            joins.len(),
            parked.len(),
            refused.len(),
        );
    }

    // A foreign contig repeats a marker as often as its bin is complete and an own one almost never,
    // so each bin weighs its own fills against its repeats at the run's contamination weight.
    fn refused(
        &self,
        bins: &HashMap<usize, HashSet<usize>>,
        joins: &[(usize, usize)],
    ) -> HashSet<usize> {
        let evidence = self.marker_evidence(bins, joins);
        let repeats = evidence.iter().filter(|seen| seen.repeats).count() as f64;
        let complete = evidence.iter().map(|seen| seen.complete).sum::<f64>();
        let foreign = (repeats / complete).min(1.0);
        let mut trade = HashMap::<usize, f64>::new();
        for seen in &evidence {
            *trade.entry(seen.bin).or_default() += match seen.repeats {
                true => -self.worth,
                false => 1.0 - foreign * (1.0 - seen.complete),
            };
        }
        info!(
            "Short contigs' markers repeat {repeats} times against {complete:.1} if all were \
             foreign, a foreign share of {foreign:.2}."
        );
        trade
            .into_iter()
            .filter(|(_, gain)| *gain <= 0.0)
            .map(|(bin, _)| bin)
            .collect()
    }

    fn marker_evidence(
        &self,
        bins: &HashMap<usize, HashSet<usize>>,
        joins: &[(usize, usize)],
    ) -> Vec<Evidence> {
        let members = bins
            .iter()
            .map(|(bin, contigs)| {
                let mut contigs = contigs.iter().copied().collect::<Vec<_>>();
                contigs.sort_unstable();
                (*bin, contigs)
            })
            .collect::<HashMap<_, _>>();
        joins
            .par_iter()
            .filter(|(contig, _)| self.quality.hit_count(&[*contig]) > 0)
            .filter_map(|(contig, bin)| {
                let rest = &members[bin];
                let mut with = rest.clone();
                with.insert(with.partition_point(|at| at < contig), *contig);
                Some(Evidence {
                    bin: *bin,
                    repeats: self.quality.repeats(&with, *contig)?,
                    complete: self.quality.score(rest).completeness / 100.0,
                })
            })
            .collect()
    }

    // The long half of each neighbourhood is already known from the nearest search, so only the
    // parked contigs are searched among themselves and the cost follows their count.
    fn neighbourhoods(&self, parked: &[usize], nearest: &Nearest) -> KnnGraph {
        let among = self.features().knn_of(
            parked,
            self.n_neighbours,
            self.seeds,
            self.knn_candidates,
            crate::embedding::KNN_ATTACH,
        );
        let width = nearest.knn.indices.ncols();
        let mut merged = KnnGraph {
            indices: ndarray::Array2::from_elem((parked.len(), width), u32::MAX),
            dists: ndarray::Array2::from_elem((parked.len(), width), f32::INFINITY),
        };
        for (row, contig) in parked.iter().enumerate() {
            let own = contig - nearest.first;
            let mut both = nearest
                .knn
                .indices
                .row(own)
                .iter()
                .zip(nearest.knn.dists.row(own))
                .map(|(at, distance)| (*at, *distance))
                .chain(
                    among
                        .indices
                        .row(row)
                        .iter()
                        .zip(among.dists.row(row))
                        .filter(|(at, _)| **at != u32::MAX)
                        .map(|(at, distance)| (parked[*at as usize] as u32, *distance)),
                )
                .filter(|(at, _)| *at != u32::MAX)
                .collect::<Vec<_>>();
            both.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
            for (slot, (at, distance)) in both.into_iter().take(width).enumerate() {
                merged.indices[[row, slot]] = at;
                merged.dists[[row, slot]] = distance;
            }
        }
        merged
    }
}
