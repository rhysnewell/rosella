use rand::Rng;

// Past a mean of 30 a rounded normal sits within a read of the Poisson, and the product method
// would take a draw per read.
const PRODUCT_LIMIT: f64 = 30.0;

pub fn poisson(rng: &mut impl Rng, mean: f64) -> f64 {
    if mean <= 0.0 {
        return 0.0;
    }
    if mean >= PRODUCT_LIMIT {
        return (mean + mean.sqrt() * normal(rng)).round().max(0.0);
    }
    let limit = (-mean).exp();
    let mut product = rng.random::<f64>();
    let mut count = 0.0;
    while product > limit {
        product *= rng.random::<f64>();
        count += 1.0;
    }
    count
}

pub fn normal(rng: &mut impl Rng) -> f64 {
    let radius = (-2.0 * (1.0 - rng.random::<f64>()).ln()).sqrt();
    radius * (std::f64::consts::TAU * rng.random::<f64>()).cos()
}

// Marsaglia & Tsang (2000), for a shape of one or more.
pub fn gamma(rng: &mut impl Rng, shape: f64) -> f64 {
    let d = shape - 1.0 / 3.0;
    let c = (9.0 * d).sqrt().recip();
    loop {
        let x = normal(rng);
        let v = (1.0 + c * x).powi(3);
        if v <= 0.0 {
            continue;
        }
        let u = 1.0 - rng.random::<f64>();
        if u.ln() < 0.5 * x * x + d - d * v + d * v.ln() {
            return d * v;
        }
    }
}
