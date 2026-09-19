use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use rayon::prelude::*;

use crate::embedding::features::ContigFeatures;
use crate::embedding::knn::KnnGraph;

struct Home {
    rank: usize,
    distance: f64,
    within_graph: bool,
}

pub fn write(
    path: &Path,
    groups: &[Vec<usize>],
    features: &ContigFeatures,
    knn: &KnnGraph,
    share: f64,
    lengths: &[usize],
    names: &[String],
) -> Result<()> {
    let n_contigs = names.len();
    let mut genome_of = vec![usize::MAX; n_contigs];
    for (genome, members) in groups.iter().enumerate() {
        for contig in members {
            genome_of[*contig] = genome;
        }
    }
    let all = (0..n_contigs).collect::<Vec<_>>();
    let metric = features.prepared(&all);
    let width = knn.indices.ncols();
    let (sigmas, rhos) = crate::embedding::fuzzy::scales(knn.dists.view(), width);
    let edge =
        crate::embedding::selective::edge_membership(knn.dists.view(), width, &sigmas, &rhos);
    let chosen = crate::embedding::selective::needy(&edge, share);

    let rows = all
        .par_iter()
        .filter(|contig| {
            let genome = genome_of[**contig];
            genome != usize::MAX && groups[genome].len() > 1
        })
        .map(|contig| {
            let contig = *contig;
            let genome = genome_of[contig];
            let home = home_of(contig, genome, &groups[genome], &genome_of, knn, &metric);
            let depth = features
                .coverage_row(contig)
                .chunks_exact(2)
                .map(|sample| sample[0])
                .sum::<f64>()
                / (features.n_samples().max(1) as f64);
            format!(
                "{}\t{}\t{:.4}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\n",
                names[contig],
                lengths[contig],
                depth,
                groups[genome].len(),
                if home.within_graph { 1 } else { 0 },
                home.rank,
                home.distance,
                knn.dists[[contig, width - 1]],
                edge[contig],
                if chosen[contig] { 1 } else { 0 }
            )
        })
        .collect::<Vec<_>>();

    let mut out = BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        out,
        "contig\tlength\tdepth\tgenome_contigs\tin_graph\thome_rank\thome_distance\tedge_distance\tedge_membership\tselected"
    )?;
    for row in rows {
        out.write_all(row.as_bytes())?;
    }
    out.flush()?;
    Ok(())
}

fn home_of(
    contig: usize,
    genome: usize,
    members: &[usize],
    genome_of: &[usize],
    knn: &KnnGraph,
    metric: &crate::embedding::metrics::prepared::PreparedAggregate,
) -> Home {
    for (rank, neighbour) in knn.indices.row(contig).iter().enumerate() {
        if genome_of[*neighbour as usize] == genome {
            return Home {
                rank,
                distance: knn.dists[[contig, rank]] as f64,
                within_graph: true,
            };
        }
    }
    let nearest = members
        .iter()
        .filter(|member| **member != contig)
        .map(|member| metric.distance(contig, *member))
        .fold(f64::INFINITY, f64::min);
    let closer = (0..genome_of.len())
        .into_par_iter()
        .filter(|other| *other != contig && metric.distance(contig, *other) < nearest)
        .count();
    Home {
        rank: closer,
        distance: nearest,
        within_graph: false,
    }
}
