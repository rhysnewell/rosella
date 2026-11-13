use rosella::recover::floor_walk::{Bar, Reach, bands};

#[test]
fn bands_double_what_is_in_play_and_never_split_a_length() {
    let lengths = [1400, 1300, 1300, 1300, 1200, 1100, 1000, 900, 800, 700];
    assert_eq!(bands(&lengths, 2), vec![0..4, 4..10]);
}

#[test]
fn the_share_bar_rises_until_the_joins_above_it_read_mostly_home() {
    let mut bar = Bar::default();
    let home = [(0.95, false, 1.0), (0.9, false, 1.0)];
    assert_eq!(bar.add(home), Some(0.5), "all home keeps every join");
    let foreign = [(0.8, true, 1.0), (0.7, true, 1.0), (0.6, true, 1.0)];
    assert_eq!(
        bar.add(foreign),
        Some(0.7),
        "pooled with the bands above, joins over 0.7 read one repeat in three"
    );
    assert_eq!(
        Bar::default().add([(0.9, true, 1.0)]),
        None,
        "no share reads under half foreign"
    );
    assert_eq!(
        Bar::default().add([]),
        None,
        "no marker evidence attaches nothing"
    );
}

#[test]
fn the_marginal_bar_stops_where_joins_at_that_share_turn_foreign() {
    let home = |share| (share, false, 1.0);
    let foreign = |share| (share, true, 1.0);
    assert_eq!(
        Bar::marginal().add([home(0.9), home(0.8), home(0.7), foreign(0.4), foreign(0.3)]),
        Some(0.4),
        "home joins above do not carry the foreign ones under them"
    );
    assert_eq!(
        Bar::marginal().add([home(0.45), home(0.35), foreign(0.2)]),
        Some(0.2),
        "the bar goes under a majority where the markers read home"
    );
    assert_eq!(
        Bar::marginal().add([foreign(0.9), home(0.8), home(0.7)]),
        Some(0.0),
        "a lone foreign join on top is pooled with the home joins under it"
    );
}

#[test]
fn a_band_is_worth_what_it_adds_on_top_of_the_bands_taken() {
    let mut reach = Reach::new(100.0, 0.0, 100);
    let first = reach.add(10);
    let second = reach.add(10);
    assert!((first - 21.0).abs() < 1e-9 && (second - 23.0).abs() < 1e-9);
}
