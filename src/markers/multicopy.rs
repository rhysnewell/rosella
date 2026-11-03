use log::debug;

use super::ContigMarkers;
use crate::quality::{FAMILY_ERROR, Scorer, upper_tail};

pub(super) const BINNED_COMPLETENESS: f64 = 50.0;

struct Tallies {
    doubled: Vec<u32>,
    held: Vec<u32>,
}

impl Tallies {
    fn new(len: usize) -> Self {
        Self {
            doubled: vec![0; len],
            held: vec![0; len],
        }
    }

    fn add(&mut self, at: usize, copies: u32) {
        self.held[at] += 1;
        self.doubled[at] += u32::from(copies >= 2);
    }

    fn flagged(&self) -> Vec<bool> {
        let trials = self.held.iter().map(|n| f64::from(*n)).sum::<f64>();
        let rate = self.doubled.iter().map(|k| f64::from(*k)).sum::<f64>() / trials.max(1.0);
        let tested = self.held.iter().filter(|n| **n > 0).count().max(1) as f64;
        self.doubled
            .iter()
            .zip(&self.held)
            .map(|(k, n)| *k > 0 && upper_tail(*k, *n, rate) < FAMILY_ERROR / tested)
            .collect()
    }
}

impl ContigMarkers {
    // A second genome duplicates every marker it carries at once, while a family that is
    // multi-copy in this community duplicates alone, far above the rate the others share.
    pub fn learn_multicopy<'a>(&self, bins: impl Iterator<Item = &'a [usize]>) {
        if self.multicopy.get().is_some() {
            return;
        }
        let mut gtdb = Tallies::new(self.set.len());
        let mut checkm = Tallies::new(self.set.checkm.len());
        for contigs in bins {
            if self.score(contigs).completeness < BINNED_COMPLETENESS {
                continue;
            }
            let Some((chosen, counts)) = self.chosen(contigs) else {
                continue;
            };
            for (marker, tally) in counts.iter().enumerate() {
                if self.set.sets.holds(chosen, marker) && tally.any > 0 {
                    gtdb.add(marker, tally.complete);
                }
            }
            let panel = &self.set.checkm;
            let Some(at) = panel
                .lineage(self.set.sets.name(chosen))
                .filter(|at| !panel.reads_gtdb(*at))
            else {
                continue;
            };
            let copies = self.lineage_copies(contigs, &counts, at, true);
            for model in panel.models(at) {
                let held = copies[*model as usize];
                if held > 0 {
                    checkm.add(*model as usize, held);
                }
            }
        }
        let mask = gtdb.flagged();
        let checkm_mask = checkm.flagged();
        debug!(
            "Markers multi-copy in this community: {} (CheckM {})",
            named(&mask, |marker| self.set.name(marker)),
            named(&checkm_mask, |model| self.set.checkm.name(model))
        );
        let _ = self.multicopy.set(mask);
        let _ = self.checkm_multicopy.set(checkm_mask);
    }

    pub(super) fn is_checkm_multicopy(&self, model: u16) -> bool {
        self.checkm_multicopy
            .get()
            .and_then(|mask| mask.get(model as usize))
            .copied()
            .unwrap_or(false)
    }
}

fn named<'a>(mask: &[bool], name: impl Fn(u16) -> &'a str) -> String {
    mask.iter()
        .enumerate()
        .filter(|(_, flagged)| **flagged)
        .map(|(at, _)| name(at as u16))
        .collect::<Vec<_>>()
        .join(",")
}
