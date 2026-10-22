use anyhow::Result;
use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::device::{self, Found, Single, Trainer};
use rosella::learned::labelled::loss::contrastive;
use rosella::learned::loss::per_positive;
use rosella::learned::net::{Net, Params};
use rosella::learned::towers::Towers;

const TEMPERATURE: f32 = 0.07;

type Losses = (f64, Vec<(f64, usize)>);

fn norm(values: impl Iterator<Item = f32>) -> f64 {
    values
        .map(|value| f64::from(value).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn norm_gap(host: &[f32], device: &[f32]) -> f64 {
    norm(host.iter().zip(device).map(|(host, device)| host - device)) / norm(host.iter().copied())
}

fn features(rng: &mut StdRng, rows: usize, width: usize) -> Array2<f32> {
    Array2::from_shape_simple_fn((rows, width), || rng.random_range(-2.0..2.0))
}

fn assert_losses((loss, losses): Losses, (device_loss, device_losses): Losses) {
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
}

// The products round in another order, so a unit at the ReLU edge can fall the other way and
// move every gradient below it a little. The gradients agree in norm, not entry by entry.
fn assert_grads(net: &str, host: &Params, device: &[Vec<f32>]) {
    for (tensor, (host, device)) in host.slices().into_iter().zip(device).enumerate() {
        let gap = norm_gap(host, device);
        assert!(
            gap < 2e-2,
            "{net} tensor {tensor}: off by {gap:.2e} in norm"
        );
    }
}

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
    let features = features(&mut rng, rows, width);
    let labels = (0..rows as u32).map(|row| row / 4).collect::<Vec<_>>();
    let drop = Some((0.2, 99));
    let pass = towers.forward(features.clone(), drop);
    let (loss, gradient, losses) = per_positive(pass.output(), &labels, TEMPERATURE);
    let (joint, depth) = towers.gradients(&pass, gradient);
    let device = trainer.gradients(&features, &labels, drop, TEMPERATURE)?;
    assert_losses((loss, losses), device);
    let [device_joint, device_depth] = trainer.grads()?;
    assert_grads("joint", &joint, &device_joint);
    assert_grads("depth", &depth, &device_depth);
    Ok(())
}

// Bins far past the per-positive loss's mate cap, crop pairs and rows with no positive at all.
#[test]
fn a_labelled_device_step_matches_the_host_gradients() -> Result<()> {
    let Found::Device(_) = device::open()? else {
        return Ok(());
    };
    let mut rng = StdRng::seed_from_u64(9);
    let (rows, width) = (600, 20);
    let mut net = Net::new(width, &[512, 512], 64, &mut rng);
    let Some(mut single) = Single::new(&net, 700, 1e-4)? else {
        return Ok(());
    };
    let features = features(&mut rng, rows, width);
    let labels = (0..rows as u32)
        .map(|row| match row {
            0..200 => row / 50,
            200..590 => row / 2,
            _ => row + 1000,
        })
        .collect::<Vec<_>>();
    let pass = net.forward(features.clone());
    let (loss, gradient, losses) = contrastive(pass.output.view(), &labels, TEMPERATURE);
    let grads = net.backward(&pass, gradient);
    let device = single.gradients(&features, &labels, TEMPERATURE)?;
    assert_losses((loss, losses), device);
    assert_grads("labelled", &grads, &single.grads()?);
    Ok(())
}
