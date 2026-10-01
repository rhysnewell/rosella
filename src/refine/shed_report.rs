use std::collections::BTreeMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::embedding::metrics::{AggregateMetric, abundance_distance, rho_between};
use crate::markers::ContigMarkers;
use crate::refine::report_context::{Context, Inputs};

/// Every bin a shed contig could be offered instead of the unbinned, with what the neighbours,
/// the claim and the markers each make of the move.
pub fn write(
    path: &Path,
    bins: &BTreeMap<usize, Vec<usize>>,
    inputs: &Inputs<'_>,
    markers: &ContigMarkers,
) -> Result<()> {
    let context = Context::of(bins, inputs.features, inputs.quality);
    let mut out = BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        out,
        "contig\tlength\town_bin\town_share\tmarkers\ttwin\tshared\t\
         twin_distance\ttwin_coverage\ttwin_composition\ttwin_rank\tbin_median\t\
         twin_comp_rank\tbin_comp_median\tbin_members\t\
         rival_bin\trival_share\tclaim\tcompletes"
    )?;
    for label in context.labels() {
        let held = &context.members[&label];
        for entry in markers.redundant_traced(held) {
            let contig = entry.contig;
            let twin = entry.twin.map_or("-", |other| inputs.names[other].as_str());
            let Some((own, rivals)) = context.rivals(contig, label, inputs) else {
                continue;
            };
            let pair = carrier_pair(inputs, &context.metric, held, contig, entry.twin);
            for rival in rivals {
                let completes = context.members.get(&rival.bin).is_some_and(|into| {
                    let mut offered = into.clone();
                    offered.push(contig);
                    offered.sort_unstable();
                    markers.completes(&offered, contig)
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
                    rival.bin,
                    rival.share,
                    rival.claim,
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
    let points = inputs.features.points(held);
    let mut distances = Vec::with_capacity(held.len());
    let mut compositions = Vec::with_capacity(held.len());
    for (position, point) in points.iter().enumerate() {
        if position != mine {
            distances.push(metric.distance(&points[mine], point));
            compositions.push(rho_between(&points[mine].composition, &point.composition));
        }
    }
    let (mine, twin) = (&points[mine], &points[twin]);
    let distance = metric.distance(mine, twin);
    let composition = rho_between(&mine.composition, &twin.composition);
    let (rank, median) = place(&mut distances, distance);
    let (comp_rank, comp_median) = place(&mut compositions, composition);
    Pair {
        distance,
        coverage: abundance_distance(&mine.abundance, &twin.abundance).0,
        composition,
        rank,
        median,
        comp_rank,
        comp_median,
        members: held.len(),
    }
}
