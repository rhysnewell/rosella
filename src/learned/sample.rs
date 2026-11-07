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
