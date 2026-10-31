pub mod bases;
pub mod bins;
pub mod orfs;
pub mod report;

#[derive(Debug, Clone, Copy, Default)]
pub struct Quality {
    pub completeness: f64,
    pub contamination: f64,
    pub set: u16,
}

impl Quality {
    pub fn score(&self, weight: f64) -> f64 {
        self.completeness - weight * self.contamination
    }

    pub fn clears(&self, bars: Bars) -> bool {
        self.completeness >= bars.completeness && self.contamination <= bars.contamination
    }
}

// The standard error of a squared worth sum under resampling the marker catalogue, by the
// infinitesimal jackknife (Jaeckel 1972), so it needs no replicates and no seed.
pub fn edge_spread(sides: impl IntoIterator<Item = (f64, Vec<(usize, f64)>)>) -> f64 {
    let mut gradient = std::collections::BTreeMap::<usize, f64>::new();
    for (sign, points) in sides {
        let count = points.len() as f64;
        let worth = points.iter().map(|(_, point)| point).sum::<f64>() / count;
        if points.is_empty() || worth <= 0.0 {
            continue;
        }
        for (marker, point) in points {
            *gradient.entry(marker).or_default() += sign * 2.0 * worth * (point - worth) / count;
        }
    }
    gradient
        .values()
        .map(|slope| slope * slope)
        .sum::<f64>()
        .sqrt()
}

pub const FAMILY_ERROR: f64 = 0.05;

pub fn upper_tail(k: u32, n: u32, p: f64) -> f64 {
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return 1.0;
    }
    let ln_choose = |i: u32| {
        libm::lgamma(f64::from(n) + 1.0)
            - libm::lgamma(f64::from(i) + 1.0)
            - libm::lgamma(f64::from(n - i) + 1.0)
    };
    (k..=n)
        .map(|i| (ln_choose(i) + f64::from(i) * p.ln() + f64::from(n - i) * (1.0 - p).ln()).exp())
        .sum()
}

#[derive(Debug, Clone, Copy)]
pub struct Bars {
    pub completeness: f64,
    pub contamination: f64,
}

pub trait Scorer: Sync {
    fn score(&self, contigs: &[usize]) -> Quality;

    // A partner that brings no feature the bin lacks cannot raise its completeness, which is
    // most pairs, so this keeps the join off the scorer.
    fn features(&self, contigs: &[usize]) -> std::collections::HashSet<u32>;

    fn points(&self, _contigs: &[usize], _weight: f64) -> Vec<(usize, f64)> {
        Vec::new()
    }
}
