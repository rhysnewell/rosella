use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
};

use anyhow::{Result, bail};
use log::{debug, trace};
use rayon::prelude::*;

use crate::clustering::graph_partition::{NodeSize, Partition, label_propagation};
use crate::clustering::leiden::{leiden, resolutions};
use crate::clustering::objective::ClusterObjective;
use crate::embedding::Graph;

const RESOLUTION_STEPS: usize = SWEEP_WIDTH;

pub const SWEEP_WIDTH: usize = 10;

/// The ladder ranks every rung on the objective, so a caller with a better judge can have
/// the whole ladder for what the winner cost.
fn ladder(mut scored: Vec<(Vec<i32>, f64)>) -> Vec<Partitioning> {
    scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
    debug!("Best validity {}", scored[0].1);
    scored
        .into_iter()
        .map(|(labels, validity)| Partitioning::from_labels(&labels, validity))
        .collect()
}

pub fn find_partitions(
    graph: &Graph,
    lengths: &[usize],
    node_size: NodeSize,
    objective: &dyn ClusterObjective,
    partition_seed: u64,
    kind: Partition,
    resolution: Option<f64>,
    theta: Option<f64>,
) -> Result<Vec<Partitioning>> {
    let _timer = crate::timing::scope("partition");
    let sized = node_size.apply(graph, lengths);
    let graph = sized.graph.as_ref();
    let sizes = sized.sizes.as_deref();
    let rank = |labels: &[i32]| objective.score_graph(graph, labels);

    if kind == Partition::LabelProp {
        let labels = label_propagation(graph, partition_seed);
        let validity = rank(&labels);
        debug!("label propagation validity {validity}");
        return Ok(vec![Partitioning::from_labels(&labels, validity)]);
    }

    let rungs = resolution.map_or_else(
        || resolutions(graph, sizes, RESOLUTION_STEPS),
        |one| vec![one],
    );
    let scored = rungs
        .par_iter()
        .map(|resolution| {
            let labels = leiden(graph, sizes, *resolution, theta, partition_seed);
            let validity = rank(&labels);
            trace!("resolution {} validity {}", resolution, validity);
            (labels, validity)
        })
        .collect::<Vec<_>>();

    if scored.is_empty() {
        anyhow::bail!("the resolution ladder produced no labelling");
    }

    Ok(ladder(scored))
}

pub fn find_best_partition(
    graph: &Graph,
    lengths: &[usize],
    node_size: NodeSize,
    objective: &dyn ClusterObjective,
    partition_seed: u64,
    kind: Partition,
    resolution: Option<f64>,
    theta: Option<f64>,
) -> Result<Partitioning> {
    Ok(find_partitions(
        graph,
        lengths,
        node_size,
        objective,
        partition_seed,
        kind,
        resolution,
        theta,
    )?
    .swap_remove(0))
}

pub struct Partitioning {
    pub cluster_map: HashMap<usize, HashSet<usize>>,
    pub outliers: HashSet<usize>,
    pub score: f64,
}

impl Partitioning {
    pub fn from_labels(labels: &[i32], score: f64) -> Self {
        let mut cluster_map: HashMap<usize, HashSet<usize>> = HashMap::new();
        let mut outliers = HashSet::new();

        for (index, label) in labels.iter().enumerate() {
            if *label < 0 {
                outliers.insert(index);
            } else {
                cluster_map
                    .entry(*label as usize)
                    .or_default()
                    .insert(index);
            }
        }

        Self {
            cluster_map,
            outliers,
            score,
        }
    }

    /// Fold another result in, renumbering its clusters so nothing collides. Ordered by
    /// lowest member, because hash order would give the same partition different bin names
    /// on every run.
    pub fn merge(&mut self, other: Partitioning) {
        let mut next_cluster_id = self.cluster_map.keys().max().map_or(0, |id| id + 1);
        let mut incoming = other.cluster_map.into_values().collect::<Vec<_>>();
        incoming
            .sort_unstable_by_key(|indices| indices.iter().min().copied().unwrap_or(usize::MAX));
        for indices in incoming {
            self.cluster_map.insert(next_cluster_id, indices);
            next_cluster_id += 1;
        }
        self.outliers = other.outliers;
        self.score = f64::NAN;
    }

    /// Map positions within a subset back to their original contig indices.
    pub fn reindex_clusters(&mut self, contig_map: HashMap<usize, usize>) {
        self.cluster_map = self
            .cluster_map
            .par_iter()
            .map(|(cluster, points)| {
                let indices = points.par_iter().map(|point| contig_map[point]).collect();
                (*cluster, indices)
            })
            .collect();

        self.outliers = self
            .outliers
            .par_iter()
            .map(|point| contig_map[point])
            .collect();
    }
}

/// Every stage happens to partition exactly, and nothing checked it. A contig in two bins
/// survives to the writer, where the label map is keyed on name and one insert silently wins.
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
