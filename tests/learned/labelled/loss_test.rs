use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::{labelled::loss::contrastive, net::Net};

const TEMPERATURE: f32 = 0.5;

fn loss_of(net: &mut Net, input: &Array2<f32>, labels: &[u32]) -> f64 {
    let pass = net.forward(input.clone());
    contrastive(pass.output.view(), labels, TEMPERATURE).0
}

#[test]
fn backward_matches_finite_differences() {
    let mut rng = StdRng::seed_from_u64(7);
    let mut net = Net::new(6, &[512, 512], 64, &mut rng);
    let input = Array2::from_shape_simple_fn((8, 6), || rng.random_range(-1.0..1.0));
    let labels = [0, 0, 1, 1, 2, 2, 3, 3];
    let pass = net.forward(input.clone());
    let (_, gradient, _) = contrastive(pass.output.view(), &labels, TEMPERATURE);
    let analytic = net
        .backward(&pass, gradient)
        .slices()
        .into_iter()
        .map(<[f32]>::to_vec)
        .collect::<Vec<_>>();
    let step = 1e-3;
    for (tensor, values) in analytic.iter().enumerate() {
        for at in [0, values.len() / 2, values.len() - 1] {
            let original = net.params.slices_mut()[tensor][at];
            net.params.slices_mut()[tensor][at] = original + step;
            let up = loss_of(&mut net, &input, &labels);
            net.params.slices_mut()[tensor][at] = original - step;
            let down = loss_of(&mut net, &input, &labels);
            net.params.slices_mut()[tensor][at] = original;
            let numeric = (up - down) / f64::from(2.0 * step);
            let given = f64::from(values[at]);
            assert!(
                (numeric - given).abs() < 1e-3 + 0.05 * numeric.abs(),
                "tensor {tensor} entry {at}: backward {given:.6} against {numeric:.6}"
            );
        }
    }
}
