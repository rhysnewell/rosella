use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::embedding::features::ContigFeatures;
use crate::embedding::knn::KnnGraph;
use crate::embedding::metrics::AggregateMetric;
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
pub fn write(path: &Path, bins: &HashMap<usize, HashSet<usize>>, inputs: Inputs<'_>) -> Result<()> {
    let metric = AggregateMetric::new(
        inputs.features.n_samples() * 2,
        inputs.features.distance_settings(),
    );
    let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut owner: HashMap<usize, usize> = HashMap::new();
    for (label, contigs) in bins {
        let mut held = contigs.iter().copied().collect::<Vec<_>>();
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
                    "{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}",
                    inputs.names[contig],
                    inputs.lengths[contig],
                    label,
                    own,
                    entry.markers,
                    twin,
                    entry.shared,
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
