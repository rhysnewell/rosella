use ndarray::linalg::general_mat_mul;
use ndarray::{Array2, ArrayView2, s};
use rand::Rng;
use rayon::prelude::*;

use crate::learned::corpus::centred_unit;

pub const RELATIVES: usize = 10;
pub const NONE: u32 = u32::MAX;
pub const CROSS: f64 = 0.5;
const CROSS_SCALE: f64 = 1.0;
const BLOCK: usize = 512;
const SWEEP: usize = 8192;

// A relative is a near composition under another label, so the depth it lends a view is another
// genome's, as a cross-mapped read would be.
pub fn relatives(composition: ArrayView2<f32>, labels: &[u32]) -> Vec<[u32; RELATIVES]> {
    let n = composition.nrows();
    (0..n.div_ceil(BLOCK))
        .into_par_iter()
        .flat_map_iter(|block| {
            let rows = block * BLOCK..((block + 1) * BLOCK).min(n);
            let query = composition.slice(s![rows.clone(), ..]);
            let mut best = vec![[(f32::NEG_INFINITY, NONE); RELATIVES]; rows.len()];
            let mut similarity = Array2::<f32>::zeros((rows.len(), SWEEP));
            for start in (0..n).step_by(SWEEP) {
                let end = (start + SWEEP).min(n);
                let mut view = similarity.slice_mut(s![.., ..end - start]);
                let others = composition.slice(s![start..end, ..]);
                general_mat_mul(1.0, &query, &others.t(), 0.0, &mut view);
                for (at, row) in view.outer_iter().enumerate() {
                    let own = labels[rows.start + at];
                    let kept = &mut best[at];
                    for (offset, value) in row.iter().enumerate() {
                        let other = start + offset;
                        if *value <= kept[RELATIVES - 1].0 || labels[other] == own {
                            continue;
                        }
                        let mut place = RELATIVES - 1;
                        while place > 0 && kept[place - 1].0 < *value {
                            kept[place] = kept[place - 1];
                            place -= 1;
                        }
                        kept[place] = (*value, other as u32);
                    }
                }
            }
            best.into_iter().map(|kept| kept.map(|(_, other)| other))
        })
        .collect()
}

pub struct Lender {
    relatives: Vec<[u32; RELATIVES]>,
    depth: Vec<f64>,
    peaks: Vec<f64>,
    samples: usize,
}

impl Lender {
    pub fn new(composition: Array2<f32>, labels: &[u32], depth: Vec<f64>, samples: usize) -> Self {
        let peaks = depth
            .chunks(samples.max(1))
            .map(|depth| depth.iter().copied().fold(1e-6, f64::max))
            .collect();
        let relatives = {
            let _timer = crate::timing::scope("learned_relatives");
            relatives(centred_unit(composition).view(), labels)
        };
        Self {
            relatives,
            depth,
            peaks,
            samples,
        }
    }

    pub fn depth(&self, row: usize) -> &[f64] {
        &self.depth[row * self.samples..(row + 1) * self.samples]
    }

    // A share of a relative's depth, scaled to the row's own peak, as reads cross-mapped from it.
    pub fn lend(&self, rng: &mut impl Rng, chance: f64, row: usize, depth: &mut [f64]) {
        if rng.random::<f64>() >= chance {
            return;
        }
        let relative = self.relatives[row][rng.random_range(0..RELATIVES)];
        let share = rng.random::<f64>() * CROSS_SCALE;
        if relative == NONE {
            return;
        }
        let relative = relative as usize;
        let scale = share * self.peaks[row] / self.peaks[relative];
        for (depth, lent) in depth.iter_mut().zip(self.depth(relative)) {
            *depth += scale * lent;
        }
    }
}
