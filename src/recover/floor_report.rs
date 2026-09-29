use std::collections::HashSet;
use std::io::Write;
use std::path::Path;

use anyhow::Result;
use log::info;
use ndarray::Array2;
use rayon::prelude::*;

use crate::cli::RecoverArgs;
use crate::defaults::SEED_STRIDE;
use crate::embedding::fuzzy::{membership, scales};
use crate::embedding::metrics::{MIN_VAR, prepared::PreparedAggregate};
use crate::seeds::sample_positions;

const EDGES: [usize; 7] = [250, 375, 500, 625, 750, 1000, 1250];
const PER_BAND: usize = 400;

struct Query {
    band: usize,
    length: usize,
    contig: usize,
    parent: Option<usize>,
}

type Neighbours = Vec<(u32, f32)>;

// A piece of a long contig has a known home, so it shows what a contig of that length that does
// belong to a binned genome looks like. The real contigs of the band are read against it.
pub fn write(args: &RecoverArgs, path: &Path) -> Result<()> {
    let tables = crate::tables::Tables::build(&super::inputs::sources(args, 1, 1))?;
    let (coverage, tnf) = (&tables.coverage, &tables.tnf);
    let lengths = &coverage.contig_lengths;
    let cutoff = args.binning.cutoff();
    let k = args.graph.n_neighbours;
    let seed = super::settings::seeds(&args.seeds).knn;
    let long = (0..lengths.len())
        .filter(|row| lengths[*row] >= cutoff)
        .collect::<Vec<_>>();

    let mut edges = EDGES
        .iter()
        .copied()
        .filter(|edge| *edge < cutoff)
        .collect::<Vec<_>>();
    edges.push(cutoff);
    let mut pieces = Vec::new();
    let mut real = Vec::new();
    for (step, window) in edges.windows(2).enumerate() {
        let band = (0..lengths.len())
            .filter(|row| (window[0]..window[1]).contains(&lengths[*row]))
            .collect::<Vec<_>>();
        let draw = |salt: u64, n: usize, k: usize| {
            sample_positions(n, k, seed.wrapping_add(salt.wrapping_mul(SEED_STRIDE)))
        };
        let chosen = draw(2 * step as u64 + 1, band.len(), PER_BAND.min(band.len()));
        let parents = draw(
            2 * step as u64 + 2,
            long.len(),
            chosen.len().min(long.len()),
        );
        for (at, position) in chosen.iter().enumerate() {
            let contig = band[*position];
            real.push(Query {
                band: window[0],
                length: lengths[contig],
                contig,
                parent: None,
            });
            let parent = long[parents[at % parents.len()]];
            pieces.push(Query {
                band: window[0],
                length: lengths[contig],
                contig: parent,
                parent: Some(parent),
            });
        }
    }

    let names = &coverage.contig_names;
    let piece_tnf = crate::kmers::kmer_counting::prefixes(
        &args.assembly,
        &pieces
            .iter()
            .map(|piece| (names[piece.contig].as_str(), piece.length))
            .collect::<Vec<_>>(),
        &tnf.kmer_sizes(),
    )?;
    let rows = long
        .iter()
        .chain(pieces.iter().map(|piece| &piece.contig))
        .chain(real.iter().map(|query| &query.contig))
        .copied()
        .collect::<Vec<_>>();
    let coverage_rows = coverage.table.select(ndarray::Axis(0), &rows);
    let mut tnf_rows = tnf.kmer_table.select(ndarray::Axis(0), &rows);
    tnf_rows
        .slice_mut(ndarray::s![long.len()..long.len() + pieces.len(), ..])
        .assign(&piece_tnf);
    let row_lengths = long
        .iter()
        .map(|row| lengths[*row])
        .chain(pieces.iter().chain(&real).map(|query| query.length))
        .collect::<Vec<_>>();
    let all = (0..rows.len()).collect::<Vec<_>>();
    let metric = PreparedAggregate::new(
        &coverage_rows,
        &tnf_rows,
        &all,
        &vec![MIN_VAR; rows.len()],
        &row_lengths,
        tables.distance,
    );
    let composition = PreparedAggregate::new(
        &coverage_rows,
        &tnf_rows,
        &all,
        &vec![MIN_VAR; rows.len()],
        &row_lengths,
        tables.distance.composition_only(),
    );
    let nearest = |row: usize, skip: Option<usize>| closest(&metric, long.len(), row, skip, k);

    let long_row = long
        .iter()
        .enumerate()
        .map(|(at, row)| (*row, at))
        .collect::<std::collections::HashMap<_, _>>();
    let parent_rows = pieces
        .iter()
        .filter_map(|piece| piece.parent.map(|parent| long_row[&parent]))
        .collect::<HashSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let homes = parent_rows
        .par_iter()
        .map(|row| {
            let set = nearest(*row, Some(*row))
                .into_iter()
                .map(|(at, _)| at)
                .collect::<HashSet<_>>();
            (*row, set)
        })
        .collect::<std::collections::HashMap<_, _>>();
    let queried = pieces
        .par_iter()
        .enumerate()
        .map(|(at, piece)| {
            let parent = long_row[&piece.parent.expect("a piece has a parent")];
            nearest(long.len() + at, Some(parent))
        })
        .chain(
            real.par_iter()
                .enumerate()
                .map(|(at, _)| nearest(long.len() + pieces.len() + at, None)),
        )
        .collect::<Vec<_>>();

    let width = queried.iter().map(Vec::len).max().unwrap_or(0);
    let mut dists = Array2::from_elem((queried.len(), width), f32::INFINITY);
    for (row, held) in queried.iter().enumerate() {
        for (slot, (_, distance)) in held.iter().enumerate() {
            dists[[row, slot]] = *distance;
        }
    }
    let (sigmas, rhos) = scales(dists.view(), width);

    let references = sample_positions(long.len(), PER_BAND.min(long.len()), seed);
    let mut by_composition = std::io::BufWriter::new(std::fs::File::create(
        path.with_extension("composition.tsv"),
    )?);
    writeln!(by_composition, "kind\tband\tcontig\tneighbours")?;
    let composed = real
        .par_iter()
        .enumerate()
        .map(|(at, query)| {
            let row = long.len() + pieces.len() + at;
            (
                "real",
                query.band,
                query.contig,
                closest(&composition, long.len(), row, None, 10),
            )
        })
        .chain(references.par_iter().map(|at| {
            (
                "long",
                cutoff,
                long[*at],
                closest(&composition, long.len(), *at, Some(*at), 10),
            )
        }))
        .collect::<Vec<_>>();
    for (kind, band, contig, held) in composed {
        let listed = held
            .iter()
            .map(|(other, _)| names[long[*other as usize]].as_str())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            by_composition,
            "{kind}\t{band}\t{}\t{listed}",
            names[contig]
        )?;
    }

    let mut out = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(out, "kind\tband\tlength\tcontig\td1\tdk\thome\tneighbours")?;
    let mut summary = std::collections::BTreeMap::<usize, [Vec<f64>; 3]>::new();
    for (at, (query, held)) in pieces.iter().chain(&real).zip(&queried).enumerate() {
        let weights = held
            .iter()
            .map(|(_, distance)| membership(*distance, rhos[at], sigmas[at]))
            .collect::<Vec<_>>();
        let total = weights.iter().sum::<f32>().max(f32::MIN_POSITIVE);
        let home = query.parent.map(|parent| {
            let set = &homes[&long_row[&parent]];
            held.iter()
                .zip(&weights)
                .filter(|((other, _), _)| set.contains(other))
                .map(|(_, weight)| weight)
                .sum::<f32>()
                / total
        });
        let (d1, dk) = (held[0].1, held[held.len() - 1].1);
        let entry = summary.entry(query.band).or_default();
        match home {
            Some(home) => {
                entry[0].push(home as f64);
                entry[1].push(d1 as f64);
            }
            None => entry[2].push(d1 as f64),
        }
        let listed = held
            .iter()
            .zip(&weights)
            .map(|((other, _), weight)| {
                format!("{}:{:.3}", names[long[*other as usize]], weight / total)
            })
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{d1:.5}\t{dk:.5}\t{}\t{listed}",
            if home.is_some() { "piece" } else { "real" },
            query.band,
            query.length,
            names[query.contig],
            home.map_or("NA".to_string(), |home| format!("{home:.3}")),
        )?;
    }
    for (band, [home, piece_d1, real_d1]) in &summary {
        let median = median(piece_d1);
        let within =
            real_d1.iter().filter(|d1| **d1 <= median).count() as f64 / real_d1.len().max(1) as f64;
        info!(
            "Floor band {band}: pieces place {:.3} of their mass home, {:.3} of pieces hold more \
             than half, real contigs within the pieces' median nearest distance {within:.3}.",
            home.iter().sum::<f64>() / home.len().max(1) as f64,
            home.iter().filter(|home| **home > 0.5).count() as f64 / home.len().max(1) as f64,
        );
    }
    Ok(())
}

fn closest(
    metric: &PreparedAggregate,
    n_long: usize,
    row: usize,
    skip: Option<usize>,
    k: usize,
) -> Neighbours {
    let mut held = (0..n_long)
        .filter(|other| Some(*other) != skip)
        .map(|other| (other as u32, metric.distance(row, other) as f32))
        .collect::<Vec<_>>();
    let k = k.min(held.len());
    held.select_nth_unstable_by(k.saturating_sub(1), |a, b| a.1.total_cmp(&b.1));
    held.truncate(k);
    held.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    held
}

fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    sorted.get(sorted.len() / 2).copied().unwrap_or(f64::NAN)
}
