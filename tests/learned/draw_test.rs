use rand::{SeedableRng, rngs::StdRng};
use rosella::learned::draw::Draw;

#[test]
fn a_plain_run_draws_every_contig_each_epoch() {
    let mut draw = Draw::new(50, false);
    let mut rng = StdRng::seed_from_u64(7);
    for _ in 0..3 {
        let mut drawn = draw.epoch(&mut rng).to_vec();
        drawn.sort_unstable();
        assert_eq!(drawn, (0..50).collect::<Vec<u32>>());
    }
}

#[test]
fn a_fast_run_skips_placed_contigs_and_keeps_hard_and_unseen_ones() {
    let mut draw = Draw::new(4, true);
    draw.record(&[0, 2, 3, 4], &[(9.0, 1), (0.0, 1), (1.0, 1), (3.0, 1)], 2);
    let mut rng = StdRng::seed_from_u64(7);
    for _ in 0..20 {
        let drawn = draw.epoch(&mut rng).to_vec();
        assert!(!drawn.contains(&0), "{drawn:?}");
        assert!(drawn.contains(&2) && drawn.contains(&3), "{drawn:?}");
    }
}
