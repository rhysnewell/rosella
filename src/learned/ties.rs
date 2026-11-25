use std::sync::LazyLock;

use ndarray::ArrayView2;
use rand::seq::SliceRandom;
use rand::{Rng, rngs::StdRng};

pub static TIES: LazyLock<bool> =
    LazyLock::new(|| std::env::var_os("ROSELLA_TIE_BATCHES").is_some());

pub fn keys(coverage: ArrayView2<f64>) -> Vec<(u32, f64)> {
    coverage
        .rows()
        .into_iter()
        .map(|row| {
            row.iter()
                .step_by(2)
                .enumerate()
                .fold((0, f64::MIN), |best, (sample, depth)| {
                    match *depth > best.1 {
                        true => (sample as u32, *depth),
                        false => best,
                    }
                })
        })
        .collect()
}

// Contigs held by one sample at one depth share a batch, so depth cannot part them and the loss
// has to part them by sequence. The batches then run in a random order.
pub fn arrange(order: &mut Vec<u32>, keys: &[(u32, f64)], batch: usize, rng: &mut StdRng) {
    if order.is_empty() || batch == 0 {
        return;
    }
    order.sort_by(|a, b| {
        let (a, b) = (keys[*a as usize], keys[*b as usize]);
        a.0.cmp(&b.0).then(a.1.total_cmp(&b.1))
    });
    let shift = rng.random_range(0..order.len());
    order.rotate_left(shift);
    let mut chunks = order.chunks(batch).map(<[u32]>::to_vec).collect::<Vec<_>>();
    let short = chunks.pop_if(|chunk| chunk.len() < batch);
    chunks.shuffle(rng);
    *order = chunks.into_iter().chain(short).flatten().collect();
}
