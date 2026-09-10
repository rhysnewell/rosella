use std::collections::{BTreeMap, HashMap};

use rayon::prelude::*;

use crate::embedding::features::ContigFeatures;
use crate::quality::Scorer;

#[derive(Debug, Default, Clone, Copy)]
pub struct ConflictLedger {
    pub examined: usize,
    pub offered: usize,
    pub ejected: usize,
    pub bases: usize,
    pub cleared: usize,
}

impl std::fmt::Display for ConflictLedger {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "{} bins over the bar, {} contigs offered, {} ejected holding {} bp, {} bins cleared",
            self.examined, self.offered, self.ejected, self.bases, self.cleared
        )
    }
}

/// Two contigs holding the same gene family are two genomes or one gene cut in half, and
/// nothing in the sequence says which. The scorer says: a real second copy is the only one
/// whose removal takes contamination down without taking completeness with it.
fn carriers(quality: &dyn Scorer, contigs: &[usize]) -> Vec<usize> {
    let mut held = HashMap::<u32, usize>::new();
    let per_contig = contigs
        .iter()
        .map(|contig| quality.features(std::slice::from_ref(contig)))
        .collect::<Vec<_>>();
    for features in &per_contig {
        for feature in features {
            *held.entry(*feature).or_default() += 1;
        }
    }
    contigs
        .iter()
        .zip(&per_contig)
        .filter(|(_, features)| {
            features
                .iter()
                .any(|feature| held.get(feature).is_some_and(|count| *count > 1))
        })
        .map(|(contig, _)| *contig)
        .collect()
}

fn peel(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    contigs: &[usize],
    bar: f64,
    min_bin_size: usize,
) -> Option<(Vec<usize>, Vec<usize>)> {
    let mut held = quality.score(contigs);
    if held.contamination <= bar {
        return None;
    }

    let mut offered = carriers(quality, contigs);
    offered.sort_unstable_by_key(|contig| (features.length(*contig), *contig));

    let mut kept = contigs.to_vec();
    let mut ejected = Vec::new();
    for contig in &offered {
        let trial = kept
            .iter()
            .copied()
            .filter(|other| other != contig)
            .collect::<Vec<_>>();
        if features.bin_size(&trial) < min_bin_size {
            break;
        }
        let after = quality.score(&trial);
        if after.completeness < held.completeness || after.contamination >= held.contamination {
            continue;
        }
        kept = trial;
        held = after;
        ejected.push(*contig);
        if held.contamination <= bar {
            break;
        }
    }
    Some((offered, ejected))
}

/// The sketch eject sees a bin holding one organism's sequence twice. A passenger from a
/// different genome shares no k-mers with anything here, and only its genes give it away.
pub fn eject_conflicts(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    bar: f64,
    min_bin_size: usize,
) -> (Vec<usize>, ConflictLedger) {
    let proposals = bins
        .par_iter()
        .filter_map(|(bin_id, contigs)| {
            peel(features, quality, contigs, bar, min_bin_size).map(|found| (*bin_id, found))
        })
        .collect::<Vec<_>>();

    let mut ledger = ConflictLedger::default();
    let mut leaving = Vec::new();
    for (bin_id, (offered, ejected)) in proposals {
        ledger.examined += 1;
        ledger.offered += offered.len();
        if ejected.is_empty() {
            continue;
        }
        let bin = bins.get_mut(&bin_id).expect("bin was read from this map");
        bin.retain(|contig| !ejected.contains(contig));
        ledger.ejected += ejected.len();
        ledger.bases += features.bin_size(&ejected);
        ledger.cleared += usize::from(quality.score(bin).contamination <= bar);
        leaving.extend(ejected);
    }

    leaving.sort_unstable();
    (leaving, ledger)
}
