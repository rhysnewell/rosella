use crate::embedding::features::ContigFeatures;
use crate::quality::Scorer;

pub const DEFAULT_COMPLETENESS: f64 = 90.0;
pub const DEFAULT_CONTAMINATION: f64 = 5.0;
pub const DEFAULT_WORTH_CONTAMINATION: f64 = 2.0;

pub const DEFAULT_RUNG_FLOOR: f64 = 0.56;

/// The rungs the loop walks to take the genomes it is sure of first and only then the ones it
/// is not. Completeness falls evenly from the full bar to the floor; contamination does not.
const CONTAMINATION_LADDER: [f64; 5] = [1.0, 1.0, 2.0, 2.0, 3.0];

/// Floor as a share of the gap between the bin floor and genome scale. Rung zero is the fixed
/// bar the single pass always used.
const LADDER: [f64; 5] = [1.0, 0.75, 0.5, 0.5, 0.5];

pub const RUNGS: usize = LADDER.len();

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

#[derive(Debug, Clone, Copy)]
pub struct Rung {
    pub floor: usize,
    pub completeness: f64,
    pub contamination: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Bars {
    pub min_bin_size: usize,
    pub completeness: f64,
    pub contamination: f64,
    pub worth: crate::quality::Worth,
    pub rung_floor: f64,
}

impl Bars {
    pub fn at(&self, top: usize, rung: usize) -> Rung {
        let share = LADDER[rung];
        let floor = self.min_bin_size
            + (share * top.saturating_sub(self.min_bin_size) as f64).round() as usize;
        let steps = (RUNGS - 1) as f64;
        let complete = 1.0 - (1.0 - self.rung_floor) * rung as f64 / steps;
        let contaminated = CONTAMINATION_LADDER[rung];
        Rung {
            floor,
            completeness: self.completeness * complete,
            contamination: self.contamination * contaminated,
        }
    }
}

pub fn judge(
    features: &ContigFeatures,
    quality: &dyn Scorer,
    contigs: &[usize],
    rung: Rung,
) -> Verdict {
    if features.bin_size(contigs) < rung.floor {
        return Verdict::TooSmall;
    }
    let held = quality.score(contigs);
    if held.contamination > rung.contamination {
        return Verdict::Contaminated;
    }
    match held.completeness >= rung.completeness {
        true => Verdict::Adopt,
        false => Verdict::Incomplete,
    }
}
