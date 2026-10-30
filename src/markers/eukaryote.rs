use log::debug;

use super::ContigMarkers;
use super::replicon::quantile;
use crate::quality::{Bars, Quality, Scorer};

const LEAST_READINGS: usize = 4;
// Tukey's far-out fence (Tukey, 1977) and MIMAG medium quality (Bowers et al., 2017).
const FAR_OUT: f64 = 3.0;
const MEDIUM_QUALITY: f64 = 50.0;

#[derive(Clone, Copy, Debug, Default)]
pub struct Reading {
    pub density: Option<f64>,
    pub quality: Quality,
}

fn far_out_low(mut values: Vec<f64>) -> Option<f64> {
    if values.len() < LEAST_READINGS {
        return None;
    }
    values.sort_unstable_by(f64::total_cmp);
    let (low, high) = (quantile(&values, 0.25), quantile(&values, 0.75));
    Some(low - FAR_OUT * (high - low))
}

// Finished bins are prokaryotes by construction, so they set how densely a prokaryote codes here.
pub fn eukaryotic(readings: &[Reading], bars: Bars) -> Vec<bool> {
    let finished = readings
        .iter()
        .filter(|reading| reading.quality.clears(bars))
        .filter_map(|reading| reading.density)
        .collect::<Vec<_>>();
    let count = finished.len();
    let Some(fence) = far_out_low(finished) else {
        return vec![false; readings.len()];
    };
    debug!("Eukaryote fence: coding density {fence:.3} ({count} finished bins)");
    readings
        .iter()
        .map(|reading| {
            reading.quality.completeness < MEDIUM_QUALITY
                && reading.density.is_some_and(|density| density < fence)
        })
        .collect()
}

impl ContigMarkers {
    pub fn eukaryotic<'a>(
        &self,
        bins: impl IntoIterator<Item = &'a [usize]>,
        bars: Bars,
    ) -> Vec<bool> {
        let readings = bins
            .into_iter()
            .map(|contigs| Reading {
                density: self.coding(contigs).and_then(|coding| coding.density()),
                quality: self.score(contigs),
            })
            .collect::<Vec<_>>();
        eukaryotic(&readings, bars)
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
