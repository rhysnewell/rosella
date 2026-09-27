use std::collections::{HashMap, HashSet};

use log::info;
use rayon::prelude::*;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::knn::{KnnGraph, nearest_in};
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

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum Class {
    Marker,
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
            .map(
                |short| match self.quality.hit_count(&[contigs[long + short]]) > 0 {
                    true => Class::Marker,
                    false => held_by(&nearest, short, &reach, &bin_of),
                },
            )
            .collect::<Vec<_>>();
        let mut counts = HashMap::<Class, usize>::new();
        for class in &classes {
            *counts.entry(*class).or_default() += 1;
        }
        info!(
            "Short contigs: {} carry markers, {} are held by two or more long contigs of one bin, \
             {} by one, {} across bins, {} by none.",
            counts.get(&Class::Marker).unwrap_or(&0),
            counts.get(&Class::Glue).unwrap_or(&0),
            counts.get(&Class::Lone).unwrap_or(&0),
            counts.get(&Class::Shared).unwrap_or(&0),
            counts.get(&Class::Loose).unwrap_or(&0),
        );
        let chosen = classes
            .iter()
            .enumerate()
            .filter(|(_, class)| matches!(class, Class::Marker | Class::Glue))
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

    // Membership is weighed the way the graph weighs an edge, so the join needs no scale of its own.
    pub(super) fn attach(
        &self,
        bins: &mut HashMap<usize, HashSet<usize>>,
        unbinned: &mut HashSet<usize>,
        parked: &[usize],
        nearest: &Nearest,
    ) {
        let _timer = crate::timing::scope("attach");
        let bin_of = bins
            .iter()
            .flat_map(|(bin, members)| members.iter().map(move |contig| (*contig, *bin)))
            .collect::<HashMap<_, _>>();
        let width = nearest.knn.indices.ncols();
        let (sigmas, rhos) = crate::embedding::fuzzy::scales(nearest.knn.dists.view(), width);
        let joins = parked
            .par_iter()
            .filter_map(|contig| {
                let row = contig - nearest.first;
                let mut mass = HashMap::<usize, f32>::new();
                let mut total = 0.0;
                for (base, distance) in nearest
                    .knn
                    .indices
                    .row(row)
                    .iter()
                    .zip(nearest.knn.dists.row(row))
                {
                    if *base == u32::MAX {
                        break;
                    }
                    let weight =
                        crate::embedding::fuzzy::membership(*distance, rhos[row], sigmas[row]);
                    total += weight;
                    if let Some(bin) = bin_of.get(&(*base as usize)) {
                        *mass.entry(*bin).or_default() += weight;
                    }
                }
                let (bin, held) = mass
                    .into_iter()
                    .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)))?;
                (held > total / 2.0).then_some((*contig, bin))
            })
            .collect::<Vec<_>>();
        info!(
            "{} of {} parked short contigs joined a bin by their own neighbours.",
            joins.len(),
            parked.len()
        );
        for (contig, bin) in joins {
            unbinned.remove(&contig);
            bins.entry(bin).or_default().insert(contig);
        }
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
