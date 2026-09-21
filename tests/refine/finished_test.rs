//! Where the majority reading of the pool's own hold-back count falls, including the tie and
//! the empty pool the stage order can hand it.

use rosella::refine::finished::Finished;

#[test]
fn a_tie_counts_as_mostly_finished() {
    assert!(Finished::new(21, 42).mostly());
    assert!(!Finished::new(20, 42).mostly());
}

#[test]
fn an_empty_pool_is_never_mostly_finished() {
    let none = Finished::default();
    assert!(!none.mostly());
    assert_eq!(none.share(), 0.0);
}

#[test]
fn the_share_is_the_held_back_fraction() {
    assert!((Finished::new(434, 600).share() - 0.7233).abs() < 1e-4);
}
