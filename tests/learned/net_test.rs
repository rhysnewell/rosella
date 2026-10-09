use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::{
    loss::per_positive,
    net::{HIDDEN, Net},
};

const TEMPERATURE: f32 = 0.5;
const DROP: Option<(f32, u64)> = Some((0.2, 77));
const LABELS: [u32; 8] = [0, 0, 1, 1, 2, 2, 3, 3];

fn loss_of(net: &mut Net, input: &Array2<f32>, drop: Option<(f32, u64)>) -> f64 {
    let pass = net.forward_dropping(input.clone(), drop);
    per_positive(pass.output.view(), &LABELS, TEMPERATURE).0
}

fn close(given: f32, numeric: f64) -> bool {
    (numeric - f64::from(given)).abs() < 1e-3 + 0.05 * numeric.abs()
}

#[test]
fn backward_matches_finite_differences() {
    let mut rng = StdRng::seed_from_u64(9);
    let mut net = Net::new(6, HIDDEN, &mut rng);
    let input = Array2::from_shape_simple_fn((8, 6), || rng.random_range(-1.0..1.0));
    let pass = net.forward_dropping(input.clone(), DROP);
    let (_, gradient) = per_positive(pass.output.view(), &LABELS, TEMPERATURE);
    let analytic = net.backward(&pass, gradient).slices().map(<[f32]>::to_vec);
    for (tensor, values) in analytic.iter().enumerate() {
        for at in [0, values.len() / 2, values.len() - 1] {
            let original = net.params.slices_mut()[tensor][at];
            net.params.slices_mut()[tensor][at] = original + 1e-3;
            let up = loss_of(&mut net, &input, DROP);
            net.params.slices_mut()[tensor][at] = original - 1e-3;
            let down = loss_of(&mut net, &input, DROP);
            net.params.slices_mut()[tensor][at] = original;
            let numeric = (up - down) / 2e-3;
            assert!(
                close(values[at], numeric),
                "tensor {tensor} entry {at}: {:.6} against {numeric:.6}",
                values[at]
            );
        }
    }
}

#[test]
fn input_gradient_matches_finite_differences() {
    let mut rng = StdRng::seed_from_u64(11);
    let mut net = Net::new(6, HIDDEN, &mut rng);
    let input = Array2::from_shape_simple_fn((8, 6), || rng.random_range(-1.0..1.0));
    let pass = net.forward_dropping(input.clone(), None);
    let (_, gradient) = per_positive(pass.output.view(), &LABELS, TEMPERATURE);
    let (_, analytic) = net.backward_input(&pass, gradient);
    for (row, column) in [(0, 0), (3, 2), (7, 5)] {
        let mut moved = input.clone();
        moved[[row, column]] += 1e-3;
        let up = loss_of(&mut net, &moved, None);
        moved[[row, column]] -= 2e-3;
        let down = loss_of(&mut net, &moved, None);
        let numeric = (up - down) / 2e-3;
        assert!(
            close(analytic[[row, column]], numeric),
            "input {row},{column}: {:.6} against {numeric:.6}",
            analytic[[row, column]]
        );
    }
}
