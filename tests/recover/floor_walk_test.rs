use rosella::recover::floor_walk::{Foreign, Reach, bands};

#[test]
fn bands_double_what_is_in_play_and_never_split_a_length() {
    let lengths = [1400, 1300, 1300, 1300, 1200, 1100, 1000, 900, 800, 700];
    assert_eq!(bands(&lengths, 2), vec![0..4, 4..10]);
}

#[test]
fn the_foreign_share_pools_from_the_cutoff_down() {
    let mut foreign = Foreign::default();
    assert!(foreign.admits(10, 40.0));
    assert!(
        foreign.admits(6, 8.0),
        "a band over half alone stays under half pooled"
    );
    assert!(!foreign.admits(30, 10.0));
    assert!(
        !Foreign::default().admits(0, 0.0),
        "no marker evidence attaches nothing"
    );
}

#[test]
fn a_band_is_worth_what_it_adds_on_top_of_the_bands_taken() {
    let mut reach = Reach::new(100.0, 0.0, 100);
    let first = reach.add(10);
    let second = reach.add(10);
    assert!((first - 21.0).abs() < 1e-9 && (second - 23.0).abs() < 1e-9);
}
