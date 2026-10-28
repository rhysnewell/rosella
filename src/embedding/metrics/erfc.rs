use std::sync::LazyLock;

const STEPS: f64 = 32.0;
const LOW: f64 = -6.0;
const HIGH: f64 = 27.25;
const DEGREE: usize = 7;
const NODES: usize = ((HIGH - LOW) * STEPS) as usize + 1;

#[repr(align(64))]
struct Row([f64; DEGREE + 1]);

// The distance calls erfc four times per sample per pair, so a Taylor table that needs no exp
// and no branch replaces the rational fits. Past LOW erfc rounds to 2 and past HIGH it is 0.
static TABLE: LazyLock<Vec<Row>> = LazyLock::new(|| {
    (0..NODES)
        .map(|at| taylor(LOW + at as f64 / STEPS))
        .collect()
});

fn taylor(z: f64) -> Row {
    let scale = 2.0 / std::f64::consts::PI.sqrt();
    let gauss = (-z * z).exp();
    let mut hermite = [0.0; DEGREE];
    hermite[0] = 1.0;
    hermite[1] = 2.0 * z;
    for n in 2..DEGREE {
        hermite[n] = 2.0 * z * hermite[n - 1] - 2.0 * (n - 1) as f64 * hermite[n - 2];
    }
    let mut row = [0.0; DEGREE + 1];
    row[0] = libm::erfc(z);
    let mut factorial = 1.0;
    for k in 1..=DEGREE {
        factorial *= k as f64;
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        row[k] = sign * scale * hermite[k - 1] * gauss / factorial;
    }
    Row(row)
}

pub fn erfc(z: f64) -> f64 {
    let z = z.clamp(LOW, HIGH);
    let at = ((z - LOW) * STEPS).round();
    let offset = z - (LOW + at / STEPS);
    let row = &TABLE[at as usize].0;
    row[..DEGREE]
        .iter()
        .rev()
        .fold(row[DEGREE], |sum, coefficient| {
            fused(sum, offset, *coefficient)
        })
}

// Without a hardware fused multiply add, mul_add is a slow library call.
#[cfg(any(target_arch = "aarch64", target_feature = "fma"))]
fn fused(sum: f64, offset: f64, coefficient: f64) -> f64 {
    sum.mul_add(offset, coefficient)
}

#[cfg(not(any(target_arch = "aarch64", target_feature = "fma")))]
fn fused(sum: f64, offset: f64, coefficient: f64) -> f64 {
    sum * offset + coefficient
}
