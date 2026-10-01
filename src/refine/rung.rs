use crate::embedding::features::ContigFeatures;
use crate::quality::{Quality, Scorer};

pub const DEFAULT_COMPLETENESS: f64 = 90.0;
pub const DEFAULT_CONTAMINATION: f64 = 5.0;
pub const DEFAULT_WORTH_CONTAMINATION: f64 = 2.0;

pub const DEFAULT_RUNG_FLOOR: f64 = 0.56;

pub const RUNGS: usize = 5;

const TIER_MULTIPLE: f64 = 2.0;

/// The rungs the loop walks to take the genomes it is sure of first and only then the ones it
/// is not. Completeness falls evenly from the full bar to the floor; the other two do not.
fn contamination_multiple(rung: usize) -> f64 {
    1.0 + (rung / 2) as f64
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    Adopt,
    TooSmall,
    Incomplete,
    Contaminated,
    Consumed,
}

impl Verdict {
    pub fn label(&self) -> &'static str {
        match self {
            Self::Adopt => "adopt",
            Self::TooSmall => "small",
            Self::Incomplete => "incomplete",
            Self::Contaminated => "contaminated",
            Self::Consumed => "consumed",
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Rung {
    pub completeness: f64,
    pub contamination: f64,
    pub min_bin_size: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct Bars {
    pub min_bin_size: usize,
    pub completeness: f64,
    pub contamination: f64,
    pub worth: f64,
    pub rung_floor: f64,
}

impl Bars {
    /// The tier a recovered genome is counted at, not the accept bar, because a rung is judged
    /// on how many genomes it would yield rather than on what the pool will take.
    pub fn tier(&self) -> f64 {
        self.contamination * TIER_MULTIPLE
    }

    /// What the run would report rather than what the pool would adopt. No contamination
    /// ceiling: a genome whose markers duplicate reads over any ceiling whole or in pieces.
    pub fn reported(&self) -> Rung {
        Rung {
            completeness: self.completeness * self.rung_floor,
            contamination: f64::INFINITY,
            ..self.at(RUNGS - 1)
        }
    }

    pub fn at(&self, rung: usize) -> Rung {
        let steps = (RUNGS - 1) as f64;
        let complete = 1.0 - (1.0 - self.rung_floor) * rung as f64 / steps;
        Rung {
            completeness: self.completeness * complete,
            contamination: self.contamination * contamination_multiple(rung),
            min_bin_size: self.min_bin_size,
        }
    }
}

pub fn judge(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    contigs: &[usize],
    rung: Rung,
) -> Verdict {
    // Checking genome scale here only let the pool adopt a small whole genome once padded past it.
    let bases = features.bin_size(contigs);
    if bases < rung.min_bin_size {
        return Verdict::TooSmall;
    }
    verdict(bases, quality.score(contigs), rung)
}

pub fn verdict(bases: usize, held: Quality, rung: Rung) -> Verdict {
    if bases < rung.min_bin_size {
        return Verdict::TooSmall;
    }
    if held.contamination > rung.contamination {
        return Verdict::Contaminated;
    }
    match held.completeness >= rung.completeness {
        true => Verdict::Adopt,
        false => Verdict::Incomplete,
    }
}
