use ndarray::Array2;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rosella::learned::{Learned, corpus::Scales, loss::per_positive, store, towers::Towers};

const LABELS: [u32; 8] = [0, 0, 1, 1, 2, 2, 3, 3];

// One step moves every weight, norm parameter and running statistic off its start, so a part read
// back into another's slot changes the embedding.
#[test]
fn a_kept_space_embeds_as_the_trained_one_and_only_under_its_key() {
    let mut rng = StdRng::seed_from_u64(3);
    let (contigs, split, samples) = (16, 10, 2);
    let composition =
        Array2::from_shape_simple_fn((contigs, split), || rng.random_range(-1.0..1.0));
    let coverage =
        Array2::from_shape_simple_fn((contigs, 2 * samples), || rng.random_range(0.5..50.0));
    let mut net = Towers::new(split + samples, split, &mut rng);
    let mut optimisers = net.optimisers(1e-4);
    let input = Array2::from_shape_simple_fn((LABELS.len(), split + samples), || {
        rng.random_range(-1.0f32..1.0)
    });
    let pass = net.forward(input, Some((0.2, 7)));
    let (_, gradient, _) = per_positive(pass.output(), &LABELS, 0.5);
    net.step(&pass, gradient, &mut optimisers, (1e-2, 0.9));
    let learned = Learned::new(net, Scales::fit(composition.view(), coverage.view()));

    let directory = tempfile::tempdir().unwrap();
    let path = store::path(directory.path(), "key");
    learned.save(&path, "key").unwrap();
    let loaded = Learned::load(&path, "key")
        .unwrap()
        .expect("the key it was kept under");
    assert_eq!(
        learned.embed(composition.view(), coverage.view()),
        loaded.embed(composition.view(), coverage.view())
    );
    assert!(Learned::load(&path, "other").unwrap().is_none());
}
