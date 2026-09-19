use std::collections::{BTreeMap, HashMap};
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::embedding::features::ContigFeatures;
use crate::embedding::knn::KnnGraph;
use crate::embedding::metrics::{AggregateMetric, metabat_with, rho};
use crate::markers::ContigMarkers;
use crate::quality::Scorer;
use crate::refine::audit::{neighbour_weight, share};
use crate::refine::recruit::{Profile, claim, needed, wanted};

pub struct Inputs<'a> {
    pub features: &'a ContigFeatures<'a>,
    pub markers: &'a ContigMarkers,
    pub knn: &'a KnnGraph,
    pub lengths: &'a [usize],
    pub names: &'a [String],
}

/// Every bin a shed contig could be offered instead of the unbinned, with what the neighbours,
/// the claim and the markers each make of the move.
pub fn write(path: &Path, bins: &BTreeMap<usize, Vec<usize>>, inputs: Inputs<'_>) -> Result<()> {
    let metric = AggregateMetric::new(
        inputs.features.n_samples() * 2,
        inputs.features.distance_settings(),
    );
    let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut owner: HashMap<usize, usize> = HashMap::new();
    for (label, contigs) in bins {
        let mut held = contigs.clone();
        held.sort_unstable();
        for contig in &held {
            owner.insert(*contig, *label);
        }
        members.insert(*label, held);
    }
    let profiles = members
        .iter()
        .filter_map(|(label, held)| {
            Profile::of(inputs.features, &metric, held).map(|profile| (*label, profile))
        })
        .collect::<HashMap<_, _>>();
    let families = members
        .iter()
        .map(|(label, held)| (*label, inputs.markers.features(held)))
        .collect::<HashMap<_, _>>();

    let mut labels = members.keys().copied().collect::<Vec<_>>();
    labels.sort_unstable();

    let mut out = BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        out,
        "contig\tlength\town_bin\town_share\tmarkers\ttwin\tshared\t\
         twin_distance\ttwin_coverage\ttwin_composition\ttwin_rank\tbin_median\t\
         twin_comp_rank\tbin_comp_median\tbin_members\t\
         rival_bin\trival_share\tclaim\tcompletes"
    )?;
    for label in labels {
        let held = &members[&label];
        for entry in inputs.markers.redundant_traced(held) {
            let contig = entry.contig;
            let twin = entry.twin.map_or("-", |other| inputs.names[other].as_str());
            let mut weights = neighbour_weight(contig, &owner, inputs.knn, inputs.lengths);
            let Some(own) = share(&mut weights, label) else {
                continue;
            };
            let total = weights.iter().map(|(_, weight)| *weight).sum::<f64>();
            let row = inputs.features.rows(&[contig]);
            let floor = inputs.features.floors(&[contig]);
            let carried = inputs.markers.features(&[contig]);
            let leaving = profiles
                .get(&label)
                .map_or(1.0, |profile| profile.to(&metric, &row[0], floor[0]));
            let leaving_odds = needed(inputs.markers, held, contig);
            let pair = carrier_pair(&inputs, &metric, held, contig, entry.twin);
            for (rival, weight) in weights.iter().filter(|(bin, _)| *bin != label) {
                let taking = profiles
                    .get(rival)
                    .map_or(1.0, |profile| profile.to(&metric, &row[0], floor[0]));
                let odds = families
                    .get(rival)
                    .map_or(1.0, |families| wanted(&carried, families));
                let completes = members.get(rival).is_some_and(|into| {
                    let mut offered = into.clone();
                    offered.push(contig);
                    offered.sort_unstable();
                    inputs.markers.completes(&offered, contig)
                });
                writeln!(
                    out,
                    "{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t\
                     {:.6}\t{:.6}\t{:.6}\t{}\t{:.6}\t{}\t{:.6}\t{}\t\
                     {}\t{:.6}\t{:.6}\t{}",
                    inputs.names[contig],
                    inputs.lengths[contig],
                    label,
                    own,
                    entry.markers,
                    twin,
                    entry.shared,
                    pair.distance,
                    pair.coverage,
                    pair.composition,
                    pair.rank,
                    pair.median,
                    pair.comp_rank,
                    pair.comp_median,
                    pair.members,
                    rival,
                    weight / total,
                    claim(taking, odds, leaving, leaving_odds),
                    u8::from(completes)
                )?;
            }
        }
    }
    out.flush()?;
    Ok(())
}

struct Pair {
    distance: f64,
    coverage: f64,
    composition: f64,
    rank: usize,
    median: f64,
    comp_rank: usize,
    comp_median: f64,
    members: usize,
}

impl Pair {
    fn unknown(members: usize) -> Self {
        Self {
            distance: f64::NAN,
            coverage: f64::NAN,
            composition: f64::NAN,
            rank: 0,
            median: f64::NAN,
            comp_rank: 0,
            comp_median: f64::NAN,
            members,
        }
    }
}

fn place(values: &mut [f64], of: f64) -> (usize, f64) {
    let rank = values.iter().filter(|other| **other < of).count();
    values.sort_unstable_by(f64::total_cmp);
    (rank, values[values.len() / 2])
}

/// The carrier pair rather than the bin, because a carrier the candidate cannot be told apart
/// from is what a second genome looks like and what the bin's own sequence rarely does.
fn carrier_pair(
    inputs: &Inputs<'_>,
    metric: &AggregateMetric,
    held: &[usize],
    contig: usize,
    twin: Option<usize>,
) -> Pair {
    let places = |target: usize| held.iter().position(|member| *member == target);
    let (Some(twin), Some(mine)) = (twin.and_then(places), places(contig)) else {
        return Pair::unknown(held.len());
    };
    let rows = inputs.features.rows(held);
    let floor = inputs.features.floors(held)[0];
    let split = inputs.features.n_samples() * 2;
    let presence = inputs.features.distance_settings().presence_fraction;
    let (mine_coverage, mine_tnf) = rows[mine].split_at(split);
    let (twin_coverage, twin_tnf) = rows[twin].split_at(split);
    let mut distances = Vec::with_capacity(held.len());
    let mut compositions = Vec::with_capacity(held.len());
    for (position, row) in rows.iter().enumerate() {
        if position != mine {
            distances.push(metric.distance(&rows[mine], row, floor, floor));
            compositions.push(rho(mine_tnf, row.split_at(split).1));
        }
    }
    let distance = metric.distance(&rows[mine], &rows[twin], floor, floor);
    let composition = rho(mine_tnf, twin_tnf);
    let (rank, median) = place(&mut distances, distance);
    let (comp_rank, comp_median) = place(&mut compositions, composition);
    Pair {
        distance,
        coverage: metabat_with(mine_coverage, twin_coverage, floor, floor, presence).0,
        composition,
        rank,
        median,
        comp_rank,
        comp_median,
        members: held.len(),
    }
}
