use anyhow::Result;
use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::device::{self, Found, Trainer};
use rosella::learned::loss::per_positive;
use rosella::learned::towers::Towers;

const TEMPERATURE: f32 = 0.07;

fn norm(values: impl Iterator<Item = f32>) -> f64 {
    values
        .map(|value| f64::from(value).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn norm_gap(host: &[f32], device: &[f32]) -> f64 {
    norm(host.iter().zip(device).map(|(host, device)| host - device)) / norm(host.iter().copied())
}

// The products round in another order, so a unit at the ReLU edge can fall the other way and
// move every gradient below it a little. The gradients agree in norm, not entry by entry.
#[test]
fn a_device_step_matches_the_host_gradients() -> Result<()> {
    let Found::Device(_) = device::open()? else {
        return Ok(());
    };
    let mut rng = StdRng::seed_from_u64(5);
    let (rows, width, split) = (602, 20, 14);
    let mut towers = Towers::new(width, split, &mut rng);
    let Some(mut trainer) = Trainer::new(&towers, 1e-4)? else {
        return Ok(());
    };
    let features = Array2::from_shape_simple_fn((rows, width), || rng.random_range(-2.0..2.0));
    let labels = (0..rows as u32).map(|row| row / 4).collect::<Vec<_>>();
    let drop = Some((0.2, 99));
    let pass = towers.forward(features.clone(), drop);
    let (loss, gradient, losses) = per_positive(pass.output(), &labels, TEMPERATURE);
    let (joint, depth) = towers.gradients(&pass, gradient);
    let (device_loss, device_losses) = trainer.gradients(&features, &labels, drop, TEMPERATURE)?;
    assert!(
        (loss - device_loss).abs() < 1e-5 * loss,
        "loss {loss} against {device_loss}"
    );
    for (at, (host, device)) in losses.iter().zip(&device_losses).enumerate() {
        assert_eq!(host.1, device.1, "row {at} pairs");
        assert!(
            (host.0 - device.0).abs() < 1e-4 * host.0.max(1.0),
            "row {at} loss"
        );
    }
    let [device_joint, device_depth] = trainer.grads()?;
    for (tower, (host, device)) in [(joint, device_joint), (depth, device_depth)]
        .into_iter()
        .enumerate()
    {
        for (tensor, (host, device)) in host.slices().into_iter().zip(&device).enumerate() {
            let gap = norm_gap(host, device);
            assert!(
                gap < 2e-2,
                "tower {tower} tensor {tensor}: off by {gap:.2e} in norm"
            );
        }
    }
    Ok(())
}
