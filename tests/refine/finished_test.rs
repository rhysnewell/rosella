use rosella::refine::finished::Finished;

// A tie counts as mostly finished, and the empty pool the stage order can hand it never does.
#[test]
fn a_tie_is_mostly_finished_and_an_empty_pool_is_not() {
    assert!(Finished::new(21, 42).mostly());
    assert!(!Finished::new(20, 42).mostly());
    assert!(!Finished::default().mostly());
    assert_eq!(Finished::default().share(), 0.0);
}
