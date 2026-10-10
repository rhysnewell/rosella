use ndarray::linalg::general_mat_mul;
use ndarray::{Array2, ArrayView2, s};
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::gemm::{gpu, product};

fn table(rng: &mut StdRng, rows: usize, columns: usize) -> Array2<f32> {
    Array2::from_shape_simple_fn((rows, columns), || rng.random_range(-1.0..1.0))
}

#[test]
fn gpu_products_match_the_cpu_on_every_layout() {
    let Some(_gpu) = gpu(true) else {
        return;
    };
    let mut rng = StdRng::seed_from_u64(3);
    let wide = table(&mut rng, 1100, 1030);
    let tall = table(&mut rng, 1030, 1100);
    let wider = table(&mut rng, 1100, 1047);
    let cases: [(ArrayView2<f32>, ArrayView2<f32>); 4] = [
        (wide.view(), tall.slice(s![.., ..1070])),
        (tall.t(), tall.slice(s![.., 5..1075])),
        (wider.slice(s![.., 17..]), tall.slice(s![.., ..1029])),
        (wide.view(), wide.t()),
    ];
    for (at, (a, b)) in cases.into_iter().enumerate() {
        let start = table(&mut rng, a.nrows(), b.ncols());
        let mut expected = start.clone();
        general_mat_mul(0.7, &a, &b, 0.3, &mut expected);
        let mut given = start.clone();
        product(0.7, a, b, 0.3, given.view_mut());
        let worst = (&given - &expected)
            .iter()
            .fold(0.0f32, |m, d| m.max(d.abs()));
        assert!(worst < 1e-2, "case {at}: off by {worst}");
    }
}
