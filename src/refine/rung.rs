use crate::embedding::features::ContigFeatures;
use crate::quality::Scorer;

pub const DEFAULT_COMPLETENESS: f64 = 90.0;
pub const DEFAULT_CONTAMINATION: f64 = 5.0;
pub const DEFAULT_DUPLICATION_BAR: f64 = 0.05;
pub const DEFAULT_WORTH_CONTAMINATION: f64 = 2.0;

pub const DEFAULT_RUNG_FLOOR: f64 = 0.56;
pub const DEFAULT_RUNG_CEILING: f64 = 3.0;

/// The rungs the loop walks to take the genomes it is sure of first and only then the ones it
/// is not. Completeness falls evenly from the full bar to the floor; contamination does not.
const CONTAMINATION_LADDER: [f64; 5] = [1.0, 1.0, 2.0, 2.0, 3.0];

/// Floor as a share of the gap between the bin floor and genome scale, and the duplication bar
/// as a multiple of its setting. Rung zero is the fixed bar the single pass always used.
const LADDER: [(f64, f64); 5] = [(1.0, 1.0), (0.75, 1.0), (0.5, 1.0), (0.5, 2.0), (0.5, 3.0)];

pub const RUNGS: usize = LADDER.len();

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    Adopt,
    TooSmall,
    Incomplete,
    Contaminated,
    Duplicated,
    Consumed,
}

impl Verdict {
    pub fn label(&self) -> &'static str {
        match self {
            Self::Adopt => "adopt",
            Self::TooSmall => "small",
            Self::Incomplete => "incomplete",
            Self::Contaminated => "contaminated",
            Self::Duplicated => "duplicated",
            Self::Consumed => "consumed",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Rung {
    pub floor: usize,
    pub bar: f64,
    pub completeness: f64,
    pub contamination: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Bars {
    pub min_bin_size: usize,
    pub duplication_bar: f64,
    pub completeness: f64,
    pub contamination: f64,
    pub worth_contamination: f64,
    pub rung_floor: f64,
    pub rung_ceiling: f64,
}

impl Bars {
    /// A scorer that reads gene counts can call a small genome whole on its own, so the run's
    /// genome scale is only a test for one that cannot.
    pub fn at(&self, top: usize, rung: usize, sees_scale: bool) -> Rung {
        let (share, multiple) = LADDER[rung];
        let floor = self.min_bin_size
            + (share * top.saturating_sub(self.min_bin_size) as f64).round() as usize;
        let steps = (RUNGS - 1) as f64;
        let complete = 1.0 - (1.0 - self.rung_floor) * rung as f64 / steps;
        let contaminated = CONTAMINATION_LADDER[rung].min(self.rung_ceiling);
        Rung {
            floor: match sees_scale {
                true => self.min_bin_size,
                false => floor,
            },
            bar: self.duplication_bar * multiple,
            completeness: self.completeness * complete,
            contamination: self.contamination * contaminated,
        }
    }
}

pub fn over_bar(features: &ContigFeatures, contigs: &[usize], bar: f64) -> bool {
    features
        .sketches()
        .and_then(|sketches| sketches.duplication(contigs))
        .is_some_and(|duplication| duplication > bar)
}

pub fn judge(
    features: &ContigFeatures,
    quality: Option<&dyn Scorer>,
    contigs: &[usize],
    rung: Rung,
) -> Verdict {
    if features.bin_size(contigs) < rung.floor {
        return Verdict::TooSmall;
    }
    match quality {
        Some(quality) => {
            let held = quality.score(contigs);
            if held.contamination > rung.contamination {
                return Verdict::Contaminated;
            }
            match held.completeness >= rung.completeness {
                true => Verdict::Adopt,
                false => Verdict::Incomplete,
            }
        }
        None => match over_bar(features, contigs, rung.bar) {
            true => Verdict::Duplicated,
            false => Verdict::Adopt,
        },
    }
}
