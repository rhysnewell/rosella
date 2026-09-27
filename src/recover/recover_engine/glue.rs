use std::collections::HashMap;

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
