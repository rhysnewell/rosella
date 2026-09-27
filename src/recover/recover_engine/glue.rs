use std::collections::{HashMap, HashSet};

use log::info;
use rayon::prelude::*;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::knn::{KnnGraph, nearest_in};
use crate::quality::Scorer;
use crate::recover::recover_engine::RecoverEngine;

const UNBINNED: u32 = u32::MAX;

pub(super) struct Glue {
    pub(super) chosen: Vec<usize>,
    pub(super) nearest: Nearest,
}

pub(super) struct Nearest {
    pub(super) first: usize,
    pub(super) knn: KnnGraph,
}

impl Nearest {
    pub(super) fn rows(&self, contigs: &[usize]) -> ndarray::Array2<u32> {
        let rows = contigs
            .iter()
            .map(|contig| contig - self.first)
            .collect::<Vec<_>>();
        self.knn.indices.select(ndarray::Axis(0), &rows)
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum Class {
    Glue,
    Lone,
    Shared,
    Loose,
}

impl RecoverEngine {
    // Held by two long contigs of one bin, a short contig ties that bin together. Held across
    // bins, gold says it fuses two genomes rather than rejoining a split one.
    pub(super) fn glue(
        &self,
        knn: &KnnGraph,
        settled: &Partitioning,
        contigs: &[usize],
        long: usize,
    ) -> Glue {
        let _timer = crate::timing::scope("glue");
        let mut bin_of = vec![UNBINNED; long];
        for (label, members) in &settled.cluster_map {
            for at in members {
                bin_of[*at] = *label as u32;
            }
        }
        let reach = (0..long)
            .into_par_iter()
            .map(|at| own_reach(knn, &bin_of, at))
            .collect::<Vec<_>>();
        let prepared = self.features().prepared(contigs);
        let nearest = nearest_in(
            knn,
            contigs.len() - long,
            knn.indices.ncols(),
            self.knn_candidates,
            self.seeds.knn,
            |short, base| prepared.distance(long + short, base),
        );
        let classes = (0..contigs.len() - long)
            .into_par_iter()
            .map(|short| held_by(&nearest, short, &reach, &bin_of))
            .collect::<Vec<_>>();
        let mut counts = HashMap::<Class, usize>::new();
        for class in &classes {
            *counts.entry(*class).or_default() += 1;
        }
        info!(
            "Short contigs: {} are held by two or more long contigs of one bin, {} by one, {} \
             across bins, {} by none.",
            counts.get(&Class::Glue).unwrap_or(&0),
            counts.get(&Class::Lone).unwrap_or(&0),
            counts.get(&Class::Shared).unwrap_or(&0),
            counts.get(&Class::Loose).unwrap_or(&0),
        );
        let chosen = classes
            .iter()
            .enumerate()
            .filter(|(_, class)| **class == Class::Glue)
            .map(|(short, _)| contigs[long + short])
            .collect();
        Glue {
            chosen,
            nearest: Nearest {
                first: contigs[long],
                knn: nearest,
            },
        }
    }

    // Only long contigs vouch for a bin, over a neighbourhood that counts every contig. Letting
    // short ones vouch fills sink bins through chains of short contigs.
    pub(super) fn attach(
        &self,
        bins: &mut HashMap<usize, HashSet<usize>>,
        unbinned: &mut HashSet<usize>,
        parked: &[usize],
        nearest: &Nearest,
        held: &KnnGraph,
    ) {
        let _timer = crate::timing::scope("attach");
        let knn = self.neighbourhoods(parked, nearest, held);
        let bin_of = bins
            .iter()
            .flat_map(|(bin, members)| members.iter().map(move |contig| (*contig, *bin)))
            .collect::<HashMap<_, _>>();
        let (sigmas, rhos) = crate::embedding::fuzzy::scales(knn.dists.view(), knn.indices.ncols());
        let joins = parked
            .par_iter()
            .filter_map(|contig| {
                let mut mass = HashMap::<usize, f32>::new();
                let mut total = 0.0;
                for (neighbour, distance) in
                    knn.indices.row(*contig).iter().zip(knn.dists.row(*contig))
                {
                    if *neighbour == u32::MAX {
                        break;
                    }
                    let weight = crate::embedding::fuzzy::membership(
                        *distance,
                        rhos[*contig],
                        sigmas[*contig],
                    );
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
        let (repeats, expected) = self.foreign_evidence(bins, &joins);
        let own = (repeats as f64) < expected / 2.0;
        info!(
            "{} of {} parked short contigs sit in one bin's neighbourhood. Their markers repeat \
             {repeats} times against {expected:.1} if all were foreign, so they {}.",
            joins.len(),
            parked.len(),
            if own { "join" } else { "stay out" }
        );
        if !own {
            return;
        }
        for (contig, bin) in joins {
            unbinned.remove(&contig);
            bins.entry(bin).or_default().insert(contig);
        }
    }

    // A foreign contig repeats a marker its bin already has as often as the bin is complete, and
    // a contig of the bin's own genome almost never does, so repeats over that sum is the share foreign.
    fn foreign_evidence(
        &self,
        bins: &HashMap<usize, HashSet<usize>>,
        joins: &[(usize, usize)],
    ) -> (usize, f64) {
        let members = bins
            .iter()
            .map(|(bin, contigs)| {
                let mut contigs = contigs.iter().copied().collect::<Vec<_>>();
                contigs.sort_unstable();
                (*bin, contigs)
            })
            .collect::<HashMap<_, _>>();
        let evidence = joins
            .par_iter()
            .filter(|(contig, _)| self.quality.hit_count(&[*contig]) > 0)
            .filter_map(|(contig, bin)| {
                let rest = &members[bin];
                let mut with = rest.clone();
                with.insert(with.partition_point(|at| at < contig), *contig);
                let repeats = self.quality.repeats(&with, *contig)?;
                let complete = self.quality.score(rest).completeness / 100.0;
                Some((usize::from(repeats), complete))
            })
            .collect::<Vec<_>>();
        (
            evidence.iter().map(|(repeats, _)| repeats).sum(),
            evidence.iter().map(|(_, complete)| complete).sum(),
        )
    }

    fn neighbourhoods(&self, parked: &[usize], nearest: &Nearest, held: &KnnGraph) -> KnnGraph {
        let width = held.indices.ncols().max(nearest.knn.indices.ncols());
        let mut start = ndarray::Array2::from_elem((self.n_contigs, width), u32::MAX);
        start
            .slice_mut(ndarray::s![..held.n_points(), ..held.indices.ncols()])
            .assign(&held.indices);
        for contig in parked {
            start
                .slice_mut(ndarray::s![*contig, ..nearest.knn.indices.ncols()])
                .assign(&nearest.knn.indices.row(contig - nearest.first));
        }
        let everyone = (0..self.n_contigs).collect::<Vec<_>>();
        self.features().knn_from(
            &everyone,
            &start,
            self.n_neighbours,
            self.seeds,
            self.knn_candidates,
        )
    }
}

fn own_reach(knn: &KnnGraph, bin_of: &[u32], at: usize) -> f32 {
    let mut last = 0.0;
    for (neighbour, distance) in knn.indices.row(at).iter().zip(knn.dists.row(at)) {
        if *neighbour == u32::MAX {
            break;
        }
        if bin_of[*neighbour as usize] != bin_of[at] {
            return *distance;
        }
        last = *distance;
    }
    last
}

fn held_by(nearest: &KnnGraph, short: usize, reach: &[f32], bin_of: &[u32]) -> Class {
    let mut bin = None;
    let mut holders = 0;
    for (base, distance) in nearest
        .indices
        .row(short)
        .iter()
        .zip(nearest.dists.row(short))
    {
        if *base == u32::MAX {
            break;
        }
        let base = *base as usize;
        if bin_of[base] == UNBINNED || *distance >= reach[base] {
            continue;
        }
        holders += 1;
        match bin {
            None => bin = Some(bin_of[base]),
            Some(held) if held != bin_of[base] => return Class::Shared,
            Some(_) => {}
        }
    }
    match holders {
        0 => Class::Loose,
        1 => Class::Lone,
        _ => Class::Glue,
    }
}
