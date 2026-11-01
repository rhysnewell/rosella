use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::labelled::corpus::read_spans;

#[test]
fn each_sample_fits_its_own_read_length() {
    let reads = [150.0, 8000.0];
    let (bins, contigs) = (6, 300);
    let mut rng = StdRng::seed_from_u64(7);
    let lengths = (0..bins * contigs)
        .map(|_| rng.random_range(2000..60000))
        .collect::<Vec<usize>>();
    let mut coverage = Array2::zeros((lengths.len(), 2 * reads.len()));
    for row in 0..lengths.len() {
        let depth = 5.0 + (row / contigs) as f64 * 7.0;
        for (sample, read) in reads.iter().enumerate() {
            let gaussian = (-2.0 * rng.random::<f64>().ln()).sqrt()
                * (std::f64::consts::TAU * rng.random::<f64>()).cos();
            coverage[[row, 2 * sample]] =
                depth + gaussian * (depth * read / lengths[row] as f64).sqrt();
        }
    }
    let labels = (0..bins)
        .map(|bin| (bin * contigs..(bin + 1) * contigs).collect())
        .collect::<Vec<Vec<usize>>>();
    let fitted = read_spans(&labels, coverage.view(), &lengths);
    for (fit, read) in fitted.iter().zip(reads) {
        assert!((fit / read - 1.0).abs() < 0.2, "fitted {fit:.0} for {read}");
    }
}
