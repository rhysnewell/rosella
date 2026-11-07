use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::{Learned, Space, corpus::Scales, towers::Towers};

// A half that holds the whole contig at near exact depth must land on the contig's own row, or
// the weight fit compares rows from two different spaces.
#[test]
fn a_windowed_half_lands_beside_its_own_contig() {
    let mut rng = StdRng::seed_from_u64(5);
    let (contigs, split, samples) = (6, 10, 3);
    let composition =
        Array2::from_shape_simple_fn((contigs, split), || rng.random_range(-1.0..1.0));
    let coverage =
        Array2::from_shape_simple_fn((contigs, 2 * samples), || rng.random_range(5.0..80.0));
    let net = Towers::new(split + samples, split, &mut rng);
    let space = Space::Windowed(Box::new(Learned::new(
        net,
        Scales::fit(composition.view(), coverage.view()),
    )));
    let whole = space.embed(composition.view(), coverage.view());
    let rows = (0..contigs)
        .map(|row| coverage.row(row).to_slice().unwrap())
        .collect::<Vec<_>>();
    let [first, _] = space.halves(
        [composition.clone(), composition.clone()],
        &rows,
        &[1_000_000_000; 6],
        42,
    );
    for contig in 0..contigs {
        let nearest = (0..contigs)
            .max_by(|a, b| {
                let dot = |other: usize| first.row(contig).dot(&whole.row(other));
                dot(*a).total_cmp(&dot(*b))
            })
            .unwrap();
        assert_eq!(nearest, contig);
        assert!(first.row(contig).dot(&whole.row(contig)) > 0.99);
    }
}
