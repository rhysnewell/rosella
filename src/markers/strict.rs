use std::collections::HashSet;

use super::ContigMarkers;
use crate::quality::{Quality, Scorer};

// GTDB picked tree markers that rarely duplicate, so its contamination reads under CheckM's on
// the same bin. Taking the higher of the two holds a bin to the bar a CheckM scorer applies.
pub struct Strict<'a> {
    pub inner: &'a dyn Scorer,
    pub markers: &'a ContigMarkers,
}

impl Strict<'_> {
    fn stricter(&self, quality: Quality, contigs: &[usize]) -> Quality {
        match self.markers.checkm_contamination(contigs, quality.set) {
            Some(contamination) if contamination > quality.contamination => Quality {
                contamination,
                ..quality
            },
            _ => quality,
        }
    }
}

impl Scorer for Strict<'_> {
    fn score(&self, contigs: &[usize]) -> Quality {
        self.stricter(self.inner.score(contigs), contigs)
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        self.inner.features(contigs)
    }

    fn points(&self, contigs: &[usize], weight: f64) -> Vec<(usize, f64)> {
        self.inner.points(contigs, weight)
    }

    fn score_on(&self, contigs: &[usize], set: u16) -> Quality {
        self.stricter(self.inner.score_on(contigs, set), contigs)
    }

    fn points_on(&self, contigs: &[usize], weight: f64, set: u16) -> Vec<(usize, f64)> {
        self.inner.points_on(contigs, weight, set)
    }
}
