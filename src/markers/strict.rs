use std::collections::HashSet;

use ndarray::Array2;

use super::ContigMarkers;
use crate::quality::{Quality, Scorer};

// GTDB picked tree markers that rarely duplicate, so its contamination reads under CheckM's on
// the same bin. CheckM's sets are bacterial and archaeal, so a bin that reads otherwise keeps GTDB's.
pub struct Strict<'a> {
    pub inner: &'a dyn Scorer,
    pub markers: &'a ContigMarkers,
    pub fence: Option<f64>,
    pub points: Option<&'a Array2<f64>>,
}

impl Strict<'_> {
    pub fn judged(&self, contigs: &[usize], quality: Quality) -> Option<Quality> {
        if self.markers.reads_eukaryotic(contigs, quality, self.fence) {
            return None;
        }
        Some(
            match self.markers.checkm_contamination(
                contigs,
                quality.set,
                self.points
                    .map(|points| {
                        super::core::outside(contigs, &self.markers.lengths, points, |contig| {
                            self.markers.carries(contig)
                        })
                    })
                    .as_ref(),
            ) {
                Some(contamination) if contamination > quality.contamination => Quality {
                    contamination,
                    ..quality
                },
                _ => quality,
            },
        )
    }

    fn stricter(&self, quality: Quality, contigs: &[usize]) -> Quality {
        self.judged(contigs, quality).unwrap_or(quality)
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
