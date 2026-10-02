use std::collections::{BTreeMap, HashSet};

use anyhow::{Result, bail};
use log::debug;
use rayon::prelude::*;

use crate::clustering::codelength::codelength_saving;
use crate::clustering::graph_partition::{Partition, label_propagation, node_degrees, sized};
use crate::clustering::leiden::{base_level, leiden_from, resolutions};
use crate::embedding::Graph;

// The ladder ranks every rung on codelength, so a caller with a better judge can have the
// whole ladder for what the winner cost.
fn ladder(mut scored: Vec<(Vec<i32>, Option<f64>, Partition)>) -> Vec<Partitioning> {
    if scored.iter().all(|held| held.1.is_some()) {
        let saving = |held: &(Vec<i32>, Option<f64>, Partition)| {
            held.1
                .filter(|saving| !saving.is_nan())
                .unwrap_or(f64::NEG_INFINITY)
        };
        scored.sort_by(|a, b| saving(b).total_cmp(&saving(a)));
        debug!("Best validity {:?}", scored[0].1);
    }
    scored
        .into_iter()
        .map(|(labels, validity, arm)| Partitioning::from_labels(&labels, validity, arm))
        .collect()
}

pub fn find_partitions(
    graph: &Graph,
    lengths: &[usize],
    partition_seed: u64,
    kind: Partition,
    rank_rungs: bool,
) -> Result<Vec<Partitioning>> {
    let sized = sized(graph, lengths);
    let graph = &sized.graph;
    let sizes = Some(sized.sizes.as_slice());
    let degrees = rank_rungs.then(|| node_degrees(graph));
    let rank = |labels: &[i32]| {
        degrees
            .as_deref()
            .map(|degrees| codelength_saving(graph, degrees, labels))
    };

    // Label propagation is one thread's work, so it runs beside the rungs instead of before them.
    let (propagated, rungs) = rayon::join(
        || {
            kind.runs_labelprop().then(|| {
                let _timer = crate::timing::scope("partition_labelprop");
                let labels = label_propagation(graph, partition_seed);
                let validity = rank(&labels);
                debug!("label propagation validity {validity:?}");
                (labels, validity, Partition::LabelProp)
            })
        },
        || {
            kind.runs_leiden()
                .then(|| leiden_rungs(graph, sizes, partition_seed, &rank))
        },
    );
    let mut scored = propagated.into_iter().collect::<Vec<_>>();
    scored.extend(rungs.into_iter().flatten());

    if scored.is_empty() {
        anyhow::bail!("the resolution ladder produced no labelling");
    }

    let mut rungs = ladder(scored);
    for rung in rungs.iter_mut() {
        rung.seed = partition_seed;
    }
    Ok(rungs)
}

fn leiden_rungs(
    graph: &Graph,
    sizes: Option<&[f64]>,
    partition_seed: u64,
    rank: &(impl Fn(&[i32]) -> Option<f64> + Sync),
) -> Vec<(Vec<i32>, Option<f64>, Partition)> {
    let _timer = crate::timing::scope("partition_leiden");
    let rungs = resolutions(graph, sizes, crate::tuning::SWEEP_WIDTH);
    let base = base_level(graph, sizes);
    let progress =
        crate::progress::counted(crate::progress::Stage::Partitioning, rungs.len() as u64);
    let scored = rungs
        .par_iter()
        .map(|resolution| {
            let labels = leiden_from(&base, *resolution, partition_seed);
            let validity = rank(&labels);
            progress.inc(1);
            debug!(
                "resolution {resolution:.3e} communities {} validity {validity:?}",
                labels.iter().collect::<HashSet<_>>().len()
            );
            (labels, validity, Partition::Leiden)
        })
        .collect();
    progress.finish_and_clear();
    scored
}

pub fn find_best_partition(
    graph: &Graph,
    lengths: &[usize],
    partition_seed: u64,
    kind: Partition,
) -> Result<Partitioning> {
    Ok(find_partitions(graph, lengths, partition_seed, kind, true)?.swap_remove(0))
}

// Members and outliers are held sorted, so every consumer reads them in one order on every run.
pub struct Partitioning {
    pub cluster_map: BTreeMap<usize, Vec<usize>>,
    pub outliers: Vec<usize>,
    pub score: Option<f64>,
    pub arm: Partition,
    pub seed: u64,
}

impl Partitioning {
    fn from_labels(labels: &[i32], score: Option<f64>, arm: Partition) -> Self {
        let mut cluster_map: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        let mut outliers = Vec::new();
        for (index, label) in labels.iter().enumerate() {
            match usize::try_from(*label) {
                Ok(label) => cluster_map.entry(label).or_default().push(index),
                Err(_) => outliers.push(index),
            }
        }
        Self {
            cluster_map,
            outliers,
            score,
            arm,
            seed: 0,
        }
    }

    // Renumbered by lowest member, so one partition gets the same bin names on every run.
    pub fn merge(&mut self, other: Partitioning) {
        let base = self.cluster_map.keys().max().map_or(0, |id| id + 1);
        let mut incoming = other.cluster_map.into_values().collect::<Vec<_>>();
        incoming.sort_by_key(|indices| indices.first().copied().unwrap_or(usize::MAX));
        for (offset, indices) in incoming.into_iter().enumerate() {
            self.cluster_map.insert(base + offset, indices);
        }
        self.outliers = other.outliers;
        self.score = None;
    }

    pub fn reindex_clusters(&mut self, subset: &[usize]) {
        let map = |points: &mut Vec<usize>| {
            for point in points.iter_mut() {
                *point = subset[*point];
            }
            points.sort_unstable();
        };
        self.cluster_map.values_mut().for_each(map);
        map(&mut self.outliers);
    }
}

// Every stage happens to partition exactly, and nothing checked it. A contig in two bins
// survives to the writer, where the label map is keyed on name and one insert silently wins.
pub fn placed_once(
    placements: impl IntoIterator<Item = usize>,
    allowed: &HashSet<usize>,
) -> Result<HashSet<usize>> {
    let mut placed = HashSet::with_capacity(allowed.len());
    for index in placements {
        if !allowed.contains(&index) {
            bail!("contig {index} was placed but is not among the contigs handed in");
        }
        if !placed.insert(index) {
            bail!("contig {index} was placed more than once");
        }
    }
    Ok(placed)
}

pub fn conserved(
    placements: impl IntoIterator<Item = usize>,
    expected: &HashSet<usize>,
) -> Result<()> {
    let placed = placed_once(placements, expected)?;
    if placed.len() != expected.len() {
        bail!(
            "{} of {} contigs came out of binning",
            placed.len(),
            expected.len()
        );
    }
    Ok(())
}
