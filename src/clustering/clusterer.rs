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
fn ladder(mut scored: Vec<(Vec<u32>, Option<f64>, Partition)>, seed: u64) -> Vec<Partitioning> {
    if scored.iter().all(|held| held.1.is_some()) {
        let saving = |held: &(Vec<u32>, Option<f64>, Partition)| {
            held.1
                .filter(|saving| !saving.is_nan())
                .unwrap_or(f64::NEG_INFINITY)
        };
        scored.sort_by(|a, b| saving(b).total_cmp(&saving(a)));
        debug!("Best validity: {:?}", scored[0].1);
    }
    scored
        .into_iter()
        .map(|(labels, validity, arm)| Partitioning::from_labels(&labels, validity, arm, seed))
        .collect()
}

// Only the labelling hangs on the seed, so every seed's rungs share one weighted graph and one
// pool, and no seed waits on another's slowest rung. Each seed keeps its own ranked ladder.
pub fn find_partitions(
    graph: &Graph,
    lengths: &[usize],
    seeds: &[u64],
    kind: Partition,
    rank_rungs: bool,
) -> Result<Vec<Partitioning>> {
    let sized = sized(graph, lengths);
    let graph = &sized.graph;
    let sizes = sized.sizes.as_slice();
    let degrees = rank_rungs.then(|| node_degrees(graph));
    let rank = |labels: &[u32]| {
        degrees
            .as_deref()
            .map(|degrees| codelength_saving(graph, degrees, labels))
    };
    let rungs = match kind.runs_leiden() {
        true => resolutions(graph, sizes, crate::tuning::SWEEP_WIDTH),
        false => Vec::new(),
    };
    let base = kind.runs_leiden().then(|| base_level(graph, sizes));
    let jobs = seeds
        .iter()
        .enumerate()
        .flat_map(|(at, seed)| {
            let propagate = kind.runs_labelprop().then_some((at, *seed, None));
            propagate
                .into_iter()
                .chain(rungs.iter().map(move |rung| (at, *seed, Some(*rung))))
        })
        .collect::<Vec<_>>();

    let _timer = crate::timing::scope("partition_leiden");
    let progress =
        crate::progress::counted(crate::progress::Stage::Partitioning, jobs.len() as u64);
    let scored = jobs
        .par_iter()
        .map(|(at, seed, rung)| {
            let held = match (rung, &base) {
                (Some(resolution), Some(base)) => {
                    let labels = leiden_from(base, *resolution, *seed);
                    let validity = rank(&labels);
                    debug!(
                        "Resolution {resolution:.3e}: {} communities, validity {validity:?}",
                        labels.iter().collect::<HashSet<_>>().len()
                    );
                    (labels, validity, Partition::Leiden)
                }
                _ => {
                    let _timer = crate::timing::scope("partition_labelprop");
                    let labels = label_propagation(graph, *seed);
                    let validity = rank(&labels);
                    debug!("Label propagation validity: {validity:?}");
                    (labels, validity, Partition::LabelProp)
                }
            };
            progress.inc(1);
            (*at, held)
        })
        .collect::<Vec<_>>();
    crate::progress::done(&progress);

    let mut ladders = seeds.iter().map(|_| Vec::new()).collect::<Vec<_>>();
    for (at, held) in scored {
        ladders[at].push(held);
    }
    let mut partitions = Vec::new();
    for (ladder_of, seed) in ladders.into_iter().zip(seeds) {
        if ladder_of.is_empty() {
            anyhow::bail!("the resolution ladder produced no labelling");
        }
        partitions.extend(ladder(ladder_of, *seed));
    }
    Ok(partitions)
}

pub fn find_best_partition(
    graph: &Graph,
    lengths: &[usize],
    partition_seed: u64,
    kind: Partition,
) -> Result<Partitioning> {
    Ok(find_partitions(graph, lengths, &[partition_seed], kind, true)?.swap_remove(0))
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
    fn from_labels(labels: &[u32], score: Option<f64>, arm: Partition, seed: u64) -> Self {
        let mut cluster_map: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (index, label) in labels.iter().enumerate() {
            cluster_map.entry(*label as usize).or_default().push(index);
        }
        Self {
            cluster_map,
            outliers: Vec::new(),
            score,
            arm,
            seed,
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
