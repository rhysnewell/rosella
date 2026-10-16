use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::{loss::per_positive, net::Net};

const TEMPERATURE: f32 = 0.5;
const DROP: Option<(f32, u64)> = Some((0.2, 77));
const LABELS: [u32; 8] = [0, 0, 1, 1, 2, 2, 3, 3];
const LAYERS: [usize; 4] = [1024, 1024, 512, 256];

fn loss_of(net: &mut Net, input: &Array2<f32>, drop: Option<(f32, u64)>) -> f64 {
    let pass = net.forward_dropping(input.clone(), drop);
    per_positive(pass.output.view(), &LABELS, TEMPERATURE).0
}

fn close(given: f32, numeric: f64) -> bool {
    (numeric - f64::from(given)).abs() < 1e-3 + 0.05 * numeric.abs()
}

// A step that crosses a ReLU kink on one side only leaves the central difference meaningless,
// so a kinked entry passes when its gradient matches the slope on either side.
fn matches(given: f32, (down, base, up): (f64, f64, f64), step: f64) -> bool {
    let (left, right) = ((base - down) / step, (up - base) / step);
    let kinked = !close(left as f32, right);
    close(given, (up - down) / (2.0 * step))
        || kinked && (close(given, left) || close(given, right))
}

#[test]
fn backward_matches_finite_differences() {
    let mut rng = StdRng::seed_from_u64(1);
    let mut net = Net::new(6, &LAYERS, 128, &mut rng);
    let input = Array2::from_shape_simple_fn((8, 6), || rng.random_range(-1.0..1.0));
    let pass = net.forward_dropping(input.clone(), DROP);
    let (_, gradient, _) = per_positive(pass.output.view(), &LABELS, TEMPERATURE);
    let analytic = net
        .backward(&pass, gradient)
        .slices()
        .into_iter()
        .map(<[f32]>::to_vec)
        .collect::<Vec<_>>();
    for (tensor, values) in analytic.iter().enumerate() {
        for at in [0, values.len() / 2, values.len() - 1] {
            let original = net.params.slices_mut()[tensor][at];
            net.params.slices_mut()[tensor][at] = original + 1e-3;
            let up = loss_of(&mut net, &input, DROP);
            net.params.slices_mut()[tensor][at] = original - 1e-3;
            let down = loss_of(&mut net, &input, DROP);
            net.params.slices_mut()[tensor][at] = original;
            let base = loss_of(&mut net, &input, DROP);
            assert!(
                matches(values[at], (down, base, up), 1e-3),
                "tensor {tensor} entry {at}: {:.6} against {:.6}",
                values[at],
                (up - down) / 2e-3
            );
        }
    }
}

#[test]
fn input_gradient_matches_finite_differences() {
    let mut rng = StdRng::seed_from_u64(11);
    let mut net = Net::new(6, &LAYERS, 128, &mut rng);
    let input = Array2::from_shape_simple_fn((8, 6), || rng.random_range(-1.0..1.0));
    let pass = net.forward_dropping(input.clone(), None);
    let (_, gradient, _) = per_positive(pass.output.view(), &LABELS, TEMPERATURE);
    let (_, analytic) = net.backward_input(&pass, gradient, 0);
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

#[test]
fn a_column_slice_of_the_input_gradient_matches_the_full_one() {
    let mut rng = StdRng::seed_from_u64(4);
    let mut net = Net::new(264, &LAYERS, 128, &mut rng);
    let input = Array2::from_shape_simple_fn((512, 264), || rng.random_range(-1.0..1.0));
    let gradient = Array2::from_shape_simple_fn((512, 128), || rng.random_range(-1.0..1.0));
    let pass = net.forward_dropping(input, DROP);
    let (_, full) = net.backward_input(&pass, gradient.clone(), 0);
    let (_, slice) = net.backward_input(&pass, gradient, 136);
    let full = full.slice(ndarray::s![.., 136..]);
    assert_eq!(slice.dim(), full.dim());
    assert!(
        slice
            .iter()
            .zip(full)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}
