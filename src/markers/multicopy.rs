use log::debug;

use super::ContigMarkers;
use crate::quality::{FAMILY_ERROR, Scorer, upper_tail};

const BINNED_COMPLETENESS: f64 = 50.0;

impl ContigMarkers {
    // A second genome duplicates every marker it carries at once, while a family that is
    // multi-copy in this community duplicates alone, far above the rate the others share.
    pub fn learn_multicopy<'a>(&self, bins: impl Iterator<Item = &'a [usize]>) {
        if self.multicopy.get().is_some() {
            return;
        }
        let mut doubled = vec![0u32; self.set.len()];
        let mut held = vec![0u32; self.set.len()];
        for contigs in bins {
            if self.score(contigs).completeness < BINNED_COMPLETENESS {
                continue;
            }
            let Some((chosen, counts)) = self.chosen(contigs) else {
                continue;
            };
            for (marker, tally) in counts.iter().enumerate() {
                if self.set.sets.holds(chosen, marker) && tally.any > 0 {
                    held[marker] += 1;
                    doubled[marker] += u32::from(tally.complete >= 2);
                }
            }
        }
        let trials = held.iter().map(|n| f64::from(*n)).sum::<f64>();
        let rate = doubled.iter().map(|k| f64::from(*k)).sum::<f64>() / trials.max(1.0);
        let tested = held.iter().filter(|n| **n > 0).count().max(1) as f64;
        let mask = doubled
            .iter()
            .zip(&held)
            .map(|(k, n)| *k > 0 && upper_tail(*k, *n, rate) < FAMILY_ERROR / tested)
            .collect::<Vec<_>>();
        debug!(
            "Markers multi-copy in this community: {}",
            mask.iter()
                .enumerate()
                .filter(|(_, flagged)| **flagged)
                .map(|(marker, _)| self.set.name(marker as u16))
                .collect::<Vec<_>>()
                .join(",")
        );
        let _ = self.multicopy.set(mask);
    }
}
