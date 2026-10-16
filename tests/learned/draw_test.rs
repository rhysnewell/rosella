use rand::{SeedableRng, rngs::StdRng};
use rosella::learned::draw::Draw;

#[test]
fn a_plain_run_draws_every_contig_each_epoch() {
    let mut draw = Draw::new(&[2000; 50], false);
    let mut rng = StdRng::seed_from_u64(7);
    for _ in 0..3 {
        let mut drawn = draw.epoch(&mut rng).to_vec();
        drawn.sort_unstable();
        assert_eq!(drawn, (0..50).collect::<Vec<u32>>());
    }
}

#[test]
fn a_fast_run_skips_placed_contigs_and_keeps_hard_and_unseen_ones() {
    let mut draw = Draw::new(&[2000; 4], true);
    draw.record(&[0, 2, 3, 4], &[(9.0, 1), (0.0, 1), (1.0, 1), (3.0, 1)], 2);
    let mut rng = StdRng::seed_from_u64(7);
    for _ in 0..20 {
        let drawn = draw.epoch(&mut rng).to_vec();
        assert!(!drawn.contains(&0), "{drawn:?}");
        assert!(drawn.contains(&2) && drawn.contains(&3), "{drawn:?}");
    }
}

#[test]
fn a_short_contig_is_judged_against_short_peers() {
    let lengths = (0..20).map(|at| 1000 + 1000 * at).collect::<Vec<_>>();
    let mut draw = Draw::new(&lengths, true);
    let labels = (0..20).collect::<Vec<u32>>();
    let rows = (0..20)
        .map(|at| (if at < 10 { 9.0 } else { 1.0 }, 1))
        .collect::<Vec<_>>();
    draw.record(&labels, &rows, 0);
    let mut rng = StdRng::seed_from_u64(7);
    let drawn = (0..200).map(|_| draw.epoch(&mut rng).len()).sum::<usize>();
    assert!(drawn > 3000, "{drawn}");
}
