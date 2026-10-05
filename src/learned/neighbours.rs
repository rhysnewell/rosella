use ndarray::linalg::general_mat_mul;
use ndarray::{Array2, ArrayView2, s};
use rayon::prelude::*;

pub const RELATIVES: usize = 10;
pub const NONE: u32 = u32::MAX;
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
