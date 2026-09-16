const SAMPLE_PAIRS: usize = 200_000;
const BANDS: usize = 12;
const FLOOR_QUANTILE: f64 = 0.05;
const SAMPLE_SEED: u64 = 0x9E3779B97F4A7C15;
const STEP: u64 = 6364136223846793005;
const OFFSET: u64 = 1442695040888963407;

// Two short contigs sit further apart than two long ones out of the same genome, because a
// kmer vector is a multinomial sample of its contig's length.
#[derive(Debug, Clone, Copy, Default)]
pub struct LengthCalibration {
    slope: f64,
    floor: f64,
    centre: f64,
}

impl LengthCalibration {
    pub fn fit(pairs: &[(f64, f64)]) -> Option<Self> {
        if pairs.len() < BANDS * 8 {
            return None;
        }
        let mut sorted = pairs.to_vec();
        sorted.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));

        let width = (sorted.len() / BANDS).max(1);
        let mut points = Vec::with_capacity(BANDS);
        for band in sorted.chunks(width) {
            // A short trailing band takes its quantile from too few points and reads high,
            // which tilts the whole fit.
            if band.len() < width {
                continue;
            }
            let mut distances = band.iter().map(|pair| pair.1).collect::<Vec<_>>();
            distances.sort_unstable_by(f64::total_cmp);
            let at = ((distances.len() - 1) as f64 * FLOOR_QUANTILE).round() as usize;
            points.push((band[band.len() / 2].0, distances[at]));
        }
        if points.len() < 3 {
            return None;
        }

        let n = points.len() as f64;
        let mean_x = points.iter().map(|point| point.0).sum::<f64>() / n;
        let mean_y = points.iter().map(|point| point.1).sum::<f64>() / n;
        let spread = points
            .iter()
            .map(|point| (point.0 - mean_x) * (point.0 - mean_x))
            .sum::<f64>();
        if spread <= 0.0 {
            return None;
        }
        let covariance = points
            .iter()
            .map(|point| (point.0 - mean_x) * (point.1 - mean_y))
            .sum::<f64>();

        let slope = (covariance / spread).max(0.0);
        let floor = mean_y - slope * mean_x;
        let centre =
            pairs.iter().map(|pair| slope * pair.0 + floor).sum::<f64>() / pairs.len() as f64;
        Some(Self {
            slope,
            floor,
            centre,
        })
    }

    // Centred, not rescaled, so the range every downstream threshold was set against survives.
    pub fn apply(&self, distance: f64, reciprocal_sum: f64) -> f64 {
        (distance - self.slope * reciprocal_sum - self.floor + self.centre).max(0.0)
    }
}

// Deterministic on purpose: no second source of seed spread on a graph that already has one.
pub fn sample_pairs(rows: usize, mut take: impl FnMut(usize, usize)) {
    if rows < 2 {
        return;
    }
    let wanted = SAMPLE_PAIRS.min(rows * (rows - 1) / 2);
    let mut state = SAMPLE_SEED;
    let mut drawn = 0;
    while drawn < wanted {
        state = state.wrapping_mul(STEP).wrapping_add(OFFSET);
        let a = (state >> 33) as usize % rows;
        state = state.wrapping_mul(STEP).wrapping_add(OFFSET);
        let b = (state >> 33) as usize % rows;
        if a == b {
            continue;
        }
        take(a.min(b), a.max(b));
        drawn += 1;
    }
}
