use std::collections::HashMap;

use rosella::recover::offered_votes::{Ballot, ballots};
use rosella::recover::origin::BinSet;

fn set(name: &str, bins: &[&[usize]]) -> BinSet {
    BinSet {
        name: name.to_string(),
        bins: bins
            .iter()
            .enumerate()
            .map(|(at, members)| (format!("{name}{at}"), members.to_vec()))
            .collect(),
    }
}

fn ballot(bins: &[(usize, f32)], say: f32) -> Ballot {
    Ballot {
        bins: bins.to_vec(),
        say,
    }
}

#[test]
fn a_set_speaks_in_a_band_as_often_as_it_bins_that_band() {
    let lengths = [2000, 2000, 2000, 2000, 1000, 900, 800, 600];
    let bin_of = HashMap::from([(0, 10), (1, 10), (2, 11), (3, 11)]);
    let sets = [
        set("a", &[&[0, 4], &[1, 2, 6], &[7]]),
        set("b", &[&[2, 3, 4]]),
    ];
    let order = [4, 5, 6, 7];
    let spans = [0..2, 2..4];
    assert_eq!(
        ballots(&order, &spans, &lengths, 1500, &sets, &bin_of),
        vec![
            ballot(&[(10, 0.5), (11, 1.0)], 1.5),
            ballot(&[], 0.5),
            ballot(&[], 1.0),
            ballot(&[], 1.0),
        ],
        "a bins half the first band and all the second, though its split bin and its bin of short \
         contigs vote for nothing. b bins nothing under 1000 bp, so it has no say below that"
    );
    assert!(ballots(&order, &spans, &lengths, 1500, &[], &bin_of).is_empty());
}

#[test]
fn a_contig_joins_only_where_some_bin_can_still_win_a_majority() {
    assert!(
        ballot(&[], 0.0).can_join(),
        "no say leaves the neighbourhood alone"
    );
    assert!(
        !ballot(&[], 1.0).can_join(),
        "a set binning the whole band caps an abstained share at a half"
    );
    assert!(
        ballot(&[], 0.02).can_join(),
        "a set that barely bins the band barely counts"
    );
    assert!(ballot(&[(7, 1.0)], 1.0).can_join());
    assert!(
        !ballot(&[(7, 1.0), (8, 1.0)], 3.0).can_join(),
        "a split vote leaves no bin past a half"
    );
    assert!(ballot(&[(7, 1.0), (7, 1.0)], 3.0).can_join());
}
