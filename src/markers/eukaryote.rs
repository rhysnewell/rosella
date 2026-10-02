use log::debug;

use super::ContigMarkers;
use super::replicon::quantile;
use crate::quality::{Bars, Quality, Scorer};

const LEAST_FINISHED: usize = 4;
// Tukey's far-out fence (Tukey, 1977) and MIMAG medium quality (Bowers et al., 2017).
const FAR_OUT: f64 = 3.0;
const MEDIUM_QUALITY: f64 = 50.0;

#[derive(Clone, Copy, Debug, Default)]
pub struct Reading {
    pub density: Option<f64>,
    pub quality: Quality,
}

// Finished bins are prokaryotes by construction, so they set how densely a prokaryote codes here.
pub fn eukaryotic(readings: &[Reading], bars: Bars) -> Vec<bool> {
    let mut finished = readings
        .iter()
        .filter(|reading| reading.quality.clears(bars))
        .filter_map(|reading| reading.density)
        .collect::<Vec<_>>();
    if finished.len() < LEAST_FINISHED {
        return vec![false; readings.len()];
    }
    finished.sort_unstable_by(f64::total_cmp);
    let (low, high) = (quantile(&finished, 0.25), quantile(&finished, 0.75));
    let fence = low - FAR_OUT * (high - low);
    debug!(
        "Eukaryotic below coding density {fence:.3}, from {} finished bins",
        finished.len()
    );
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
}
