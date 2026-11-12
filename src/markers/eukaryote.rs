use std::collections::HashSet;

use log::debug;

use super::ContigMarkers;
use super::replicon::quantile;
use crate::quality::{Bars, Quality, Scorer};

pub(super) const LEAST_READINGS: usize = 4;
// Tukey's far-out fence (Tukey, 1977) and MIMAG medium quality (Bowers et al., 2017).
const FAR_OUT: f64 = 3.0;
pub(super) const MEDIUM_QUALITY: f64 = 50.0;

#[derive(Clone, Copy, Debug, Default)]
pub struct Reading {
    pub density: Option<f64>,
    pub quality: Quality,
}

pub(super) fn far_out_low(mut values: Vec<f64>) -> Option<f64> {
    if values.len() < LEAST_READINGS {
        return None;
    }
    values.sort_unstable_by(f64::total_cmp);
    let (low, high) = (quantile(&values, 0.25), quantile(&values, 0.75));
    Some(low - FAR_OUT * (high - low))
}

// Finished bins are prokaryotes by construction, so they set how densely a prokaryote codes here.
fn fence(readings: &[Reading], bars: Bars) -> Option<f64> {
    let finished = readings
        .iter()
        .filter(|reading| reading.quality.clears(bars))
        .filter_map(|reading| reading.density)
        .collect::<Vec<_>>();
    let count = finished.len();
    let fence = far_out_low(finished)?;
    debug!("Eukaryote fence: coding density {fence:.3} ({count} finished bins)");
    Some(fence)
}

fn below(reading: &Reading, fence: f64) -> bool {
    reading.quality.completeness < MEDIUM_QUALITY
        && reading.density.is_some_and(|density| density < fence)
}

pub fn eukaryotic(readings: &[Reading], bars: Bars) -> Vec<bool> {
    let Some(fence) = fence(readings, bars) else {
        return vec![false; readings.len()];
    };
    readings
        .iter()
        .map(|reading| below(reading, fence))
        .collect()
}

impl ContigMarkers {
    fn reading(&self, contigs: &[usize], quality: Quality) -> (bool, Reading) {
        let panel = self.set.sets.eukaryotic(usize::from(quality.set));
        let density = self.coding(contigs).and_then(|coding| coding.density());
        let density = density.filter(|_| !panel);
        (panel, Reading { density, quality })
    }

    pub fn eukaryotic<'a>(
        &self,
        bins: impl IntoIterator<Item = &'a [usize]>,
        bars: Bars,
    ) -> Vec<bool> {
        let (panel, readings): (Vec<_>, Vec<_>) = bins
            .into_iter()
            .map(|contigs| self.reading(contigs, self.score(contigs)))
            .unzip();
        eukaryotic(&readings, bars)
            .into_iter()
            .zip(panel)
            .map(|(called, panel)| called || panel)
            .collect()
    }

    pub fn eukaryote_fence<'a>(
        &self,
        bins: impl IntoIterator<Item = &'a [usize]>,
        bars: Bars,
    ) -> Option<f64> {
        let readings = bins
            .into_iter()
            .map(|contigs| self.reading(contigs, self.score(contigs)).1)
            .collect::<Vec<_>>();
        fence(&readings, bars)
    }

    pub fn reads_eukaryotic(
        &self,
        contigs: &[usize],
        quality: Quality,
        fence: Option<f64>,
    ) -> bool {
        let (panel, reading) = self.reading(contigs, quality);
        panel || fence.is_some_and(|fence| below(&reading, fence))
    }

    // A fungal paralog of a bacterial marker says the marker is there, not that a second genome
    // is, so a whole fungus would otherwise score below any half of it.
    pub fn count_sparse_as_presence(&mut self) {
        let carriers = (0..self.per_contig.len())
            .filter(|contig| self.carries(*contig))
            .collect::<Vec<_>>();
        let mut sparse = vec![false; self.per_contig.len()];
        for (contig, kept) in carriers.iter().zip(self.prokaryotic(&carriers)) {
            sparse[*contig] = !kept;
        }
        debug!(
            "Sparse carriers: {} of {} count as presence only",
            sparse.iter().filter(|sparse| **sparse).count(),
            carriers.len()
        );
        self.sparse = sparse;
    }

    // Fungi carry homologs of the same few bacterial markers, so two whole fungi share most of the
    // smaller one's families where two pieces of one fungus split them.
    pub fn shared(&self, one: &[usize], other: &[usize]) -> Option<f64> {
        let (set, _) = self.chosen(&[one, other].concat())?;
        let present = |contigs: &[usize]| {
            self.counts(contigs)
                .iter()
                .enumerate()
                .filter(|(marker, tally)| tally.any > 0 && self.set.sets.holds(set, *marker))
                .map(|(marker, _)| marker)
                .collect::<HashSet<_>>()
        };
        let (one, other) = (present(one), present(other));
        let fewer = one.len().min(other.len());
        (fewer > 0).then(|| one.intersection(&other).count() as f64 / fewer as f64)
    }

    pub(super) fn presence_only(&self, contig: usize) -> bool {
        self.sparse.get(contig).copied().unwrap_or(false)
    }

    // Introns leave a fungal contig coding sparsely, and its paralog of a bacterial marker would
    // part one genome from itself. Most carriers are prokaryotes, so they set the fence.
    pub fn prokaryotic(&self, contigs: &[usize]) -> Vec<bool> {
        let densities = contigs
            .iter()
            .map(|contig| self.coding(&[*contig]).and_then(|coding| coding.density()))
            .collect::<Vec<_>>();
        let fence = far_out_low(densities.iter().flatten().copied().collect());
        densities
            .iter()
            .map(|density| match (density, fence) {
                (Some(density), Some(fence)) => *density >= fence,
                _ => true,
            })
            .collect()
    }
}
