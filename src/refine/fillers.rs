use std::collections::{HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::quality::Scorer;

// A contig from outside the parent that brings a marker the parent lacked raises completeness
// whatever genome it is from, so only the parent's own depth range can vouch for it.
pub fn eject_fillers(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    dissolved: &[(usize, Vec<usize>)],
    promoted: &mut [Vec<usize>],
) -> Vec<usize> {
    let origin = dissolved
        .iter()
        .flat_map(|(bin, contigs)| contigs.iter().map(|contig| (*contig, *bin)))
        .collect::<HashMap<_, _>>();
    let mut ejected = Vec::new();
    for claim in promoted.iter_mut() {
        let mut taken = HashMap::<usize, usize>::new();
        for contig in claim.iter() {
            if let Some(bin) = origin.get(contig) {
                *taken.entry(*bin).or_default() += features.length(*contig);
            }
        }
        let Some((parent, held)) = taken
            .into_iter()
            .max_by_key(|(bin, held)| (*held, std::cmp::Reverse(*bin)))
        else {
            continue;
        };
        if 2 * held <= features.bin_size(claim) {
            continue;
        }
        let (core, outside): (Vec<usize>, Vec<usize>) = claim
            .iter()
            .partition(|contig| origin.get(contig) == Some(&parent));
        let held = quality.features(&core);
        let fillers = outside
            .into_iter()
            .filter(|contig| {
                quality
                    .features(&[*contig])
                    .iter()
                    .any(|marker| !held.contains(marker))
            })
            .collect::<Vec<_>>();
        if fillers.is_empty() {
            continue;
        }
        let centre = centre(features, &core);
        let reach = core
            .iter()
            .map(|contig| departure(features, *contig, &centre))
            .fold(0.0, f64::max);
        let out = fillers
            .into_iter()
            .filter(|contig| departure(features, *contig, &centre) > reach)
            .collect::<HashSet<_>>();
        if out.is_empty() {
            continue;
        }
        claim.retain(|contig| !out.contains(contig));
        ejected.extend(out);
    }
    ejected.sort_unstable();
    ejected
}

fn depth(features: &ContigFeatures, contig: usize, sample: usize) -> f64 {
    features.coverage_row(contig)[2 * sample]
}

fn centre(features: &ContigFeatures, core: &[usize]) -> Vec<Option<f64>> {
    let bases = features.bin_size(core) as f64;
    (0..features.n_samples())
        .map(|sample| {
            let depth = core
                .iter()
                .map(|contig| features.length(*contig) as f64 * depth(features, *contig, sample))
                .sum::<f64>()
                / bases;
            (depth > 0.0).then(|| depth.ln())
        })
        .collect()
}

fn departure(features: &ContigFeatures, contig: usize, centre: &[Option<f64>]) -> f64 {
    centre
        .iter()
        .enumerate()
        .filter_map(|(sample, centre)| centre.map(|centre| (sample, centre)))
        .map(|(sample, centre)| match depth(features, contig, sample) {
            depth if depth > 0.0 => (depth.ln() - centre).abs(),
            _ => f64::INFINITY,
        })
        .fold(0.0, f64::max)
}
