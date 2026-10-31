use std::collections::{BTreeMap, HashSet};

use super::{ContigMarkers, Hit, Tally};
use crate::quality::{Quality, Scorer};

pub struct Fair<'a>(pub &'a ContigMarkers);

impl Scorer for Fair<'_> {
    fn score(&self, contigs: &[usize]) -> Quality {
        self.0.score_by(contigs, true)
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        self.0.features(contigs)
    }

    fn points(&self, contigs: &[usize], weight: f64) -> Vec<(usize, f64)> {
        self.0.points_by(contigs, weight, true)
    }
}

impl ContigMarkers {
    pub fn count_fair(&mut self) {
        self.fair = true;
    }

    pub(super) fn score_by(&self, contigs: &[usize], fair: bool) -> Quality {
        let Some((chosen, counts)) = self.chosen_by(contigs, fair) else {
            return Quality::default();
        };
        match self.gtdb_on(&counts, chosen) {
            reading if reading.markers == 0 => Quality::default(),
            reading => reading.quality(chosen),
        }
    }

    // A gene cut by a contig end is a copy of its own, unless another cut gene covers the other
    // end of the model, when the two can be one gene broken across two contigs.
    pub(super) fn pair_cuts(&self, contigs: &[usize], counts: &mut [Tally]) {
        let mut held = BTreeMap::<u16, Vec<&Hit>>::new();
        for hit in contigs
            .iter()
            .filter(|contig| !self.presence_only(**contig))
            .flat_map(|contig| &self.per_contig[*contig])
        {
            if hit.partial && counts[hit.marker as usize].any >= 2 {
                held.entry(hit.marker).or_default().push(hit);
            }
        }
        for (marker, mut cuts) in held {
            cuts.sort_by_key(|hit| hit.place.reach.model_from);
            let mut genes: Vec<(&Hit, bool)> = Vec::new();
            for cut in cuts {
                match genes
                    .iter_mut()
                    .find(|(first, paired)| !*paired && !cut.same_part(first))
                {
                    Some(gene) => gene.1 = true,
                    None => genes.push((cut, false)),
                }
            }
            counts[marker as usize].complete += genes.len() as u32;
        }
    }
}
