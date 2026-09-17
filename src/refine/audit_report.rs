use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::embedding::features::ContigFeatures;
use crate::embedding::knn::KnnGraph;
use crate::embedding::metrics::AggregateMetric;
use crate::quality::Scorer;
use crate::refine::audit::{neighbour_weight, share};
use crate::refine::recruit::{Profile, claim, needed, wanted};

pub fn write(
    path: &Path,
    bins: &HashMap<usize, HashSet<usize>>,
    features: &ContigFeatures,
    quality: &dyn Scorer,
    knn: &KnnGraph,
    lengths: &[usize],
    names: &[String],
) -> Result<()> {
    let metric = AggregateMetric::new(features.n_samples() * 2, features.distance_settings());
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
            Profile::of(features, &metric, held).map(|profile| (*label, profile))
        })
        .collect::<HashMap<_, _>>();
    let families = members
        .iter()
        .map(|(label, held)| (*label, quality.features(held)))
        .collect::<HashMap<_, _>>();

    let mut audited = owner
        .iter()
        .map(|(contig, label)| (*contig, *label))
        .collect::<Vec<_>>();
    audited.sort_unstable();

    let mut out = BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        out,
        "contig\tlength\town_bin\town_share\trival_bin\trival_share\tclaim"
    )?;
    for (contig, label) in audited {
        let mut weights = neighbour_weight(contig, &owner, knn, lengths);
        let Some(own) = share(&mut weights, label) else {
            continue;
        };
        let total = weights.iter().map(|(_, weight)| *weight).sum::<f64>();
        let row = features.rows(&[contig]);
        let floor = features.floors(&[contig]);
        let held = quality.features(&[contig]);
        let leaving = profiles
            .get(&label)
            .map_or(1.0, |profile| profile.to(&metric, &row[0], floor[0]));
        let leaving_odds = members
            .get(&label)
            .map_or(1.0, |donor| needed(quality, donor, contig));
        for (rival, weight) in weights.iter().filter(|(bin, _)| *bin != label) {
            let taking = profiles
                .get(rival)
                .map_or(1.0, |profile| profile.to(&metric, &row[0], floor[0]));
            let odds = families
                .get(rival)
                .map_or(1.0, |families| wanted(&held, families));
            writeln!(
                out,
                "{}\t{}\t{}\t{:.6}\t{}\t{:.6}\t{:.6}",
                names[contig],
                lengths[contig],
                label,
                own,
                rival,
                weight / total,
                claim(taking, odds, leaving, leaving_odds)
            )?;
        }
    }
    out.flush()?;
    Ok(())
}
