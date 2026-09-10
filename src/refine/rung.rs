use crate::embedding::features::ContigFeatures;
use crate::quality::Scorer;

pub const DEFAULT_COMPLETENESS: f64 = 90.0;
pub const DEFAULT_CONTAMINATION: f64 = 5.0;

/// The rungs the loop walks to take the genomes it is sure of first and only then the ones it
/// is not.
const QUALITY_LADDER: [(f64, f64); 5] = [
    (1.0, 1.0),
    (0.89, 1.0),
    (0.78, 2.0),
    (0.67, 2.0),
    (0.56, 3.0),
];

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
}

impl Bars {
    /// A scorer that reads gene counts can call a small genome whole on its own, so the run's
    /// genome scale is only a test for one that cannot.
    pub fn at(&self, top: usize, rung: usize, sees_scale: bool) -> Rung {
        let (share, multiple) = LADDER[rung];
        let floor = self.min_bin_size
            + (share * top.saturating_sub(self.min_bin_size) as f64).round() as usize;
        let (complete, contaminated) = QUALITY_LADDER[rung];
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
