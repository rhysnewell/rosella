use std::{fs::File, io::Write, path::Path};

use anyhow::Result;
use ndarray::{ArrayBase, Data, Ix2};
use rayon::prelude::*;

use crate::clustering::{
    codelength::codelength_and_saving,
    graph_partition::{NodeSize, label_propagation},
    leiden::{leiden, resolutions},
    modularity::modularity,
    objective::{ClusterObjective, Dbcv, EVALUATION_GAMMA, EmbeddingSample},
    stability::mean_pairwise,
};
use crate::embedding::Graph;

pub struct LadderRow {
    source: &'static str,
    rung: usize,
    resolution: f64,
    nodes: usize,
    communities: usize,
    modularity: f64,
    codelength: f64,
    saving: f64,
    dbcv: f64,
    mean_ari: f64,
}

fn measure(
    source: &'static str,
    rung: usize,
    resolution: f64,
    labellings: &[Vec<i32>],
    graph: &Graph,
    contigs: &[usize],
    sample: &EmbeddingSample,
    dbcv: &Dbcv<'_>,
) -> LadderRow {
    let labels = &labellings[0];
    let communities = labels
        .iter()
        .copied()
        .filter(|label| *label >= 0)
        .max()
        .map_or(0, |highest| highest as usize + 1);

    let (codelength, saving) = codelength_and_saving(graph, labels);

    LadderRow {
        source,
        rung,
        resolution,
        nodes: graph.rows(),
        communities,
        modularity: modularity(graph, labels, EVALUATION_GAMMA),
        codelength,
        saving,
        dbcv: dbcv.score(sample, contigs, labels),
        mean_ari: mean_pairwise(labellings),
    }
}

/// Every rung scored by every candidate criterion, so the one whose peak matches the assembly
/// can be picked without scoring a bin. Labelprop joins on the same scale as the reference.
pub fn rows<S: Data<Elem = f64> + Sync>(
    graph: &Graph,
    embeddings: &ArrayBase<S, Ix2>,
    contigs: &[usize],
    lengths: &[usize],
    node_size: NodeSize,
    dbcv: &Dbcv<'_>,
    sample_seed: u64,
    partition_seed: u64,
    seeds: usize,
    steps: usize,
    theta: Option<f64>,
) -> Vec<LadderRow> {
    let sample = EmbeddingSample::new(embeddings.view(), sample_seed);
    let sized = node_size.apply(graph, lengths);
    let graph = sized.graph.as_ref();
    let sizes = sized.sizes.as_deref();
    let seeds = seeds.max(1);
    let at = |offset: usize| partition_seed.wrapping_add(offset as u64);

    let mut rows = resolutions(graph, sizes, steps)
        .par_iter()
        .enumerate()
        .map(|(index, resolution)| {
            let labellings = (0..seeds)
                .map(|offset| leiden(graph, sizes, *resolution, theta, at(offset)))
                .collect::<Vec<_>>();
            measure(
                "leiden",
                index + 1,
                *resolution,
                &labellings,
                graph,
                contigs,
                &sample,
                dbcv,
            )
        })
        .collect::<Vec<_>>();

    let propagated = (0..seeds)
        .map(|offset| label_propagation(graph, at(offset)))
        .collect::<Vec<_>>();
    rows.push(measure(
        "labelprop",
        0,
        f64::NAN,
        &propagated,
        graph,
        contigs,
        &sample,
        dbcv,
    ));

    rows
}

pub fn write(rows: &[LadderRow], path: &Path) -> Result<()> {
    let mut out = File::create(path)?;
    writeln!(
        out,
        "source\trung\tresolution\tnodes\tcommunities\tmodularity\tcodelength\tsaving\tdbcv\tmean_ari"
    )?;
    for row in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.source,
            row.rung,
            row.resolution,
            row.nodes,
            row.communities,
            row.modularity,
            row.codelength,
            row.saving,
            row.dbcv,
            row.mean_ari
        )?;
    }
    Ok(())
}
