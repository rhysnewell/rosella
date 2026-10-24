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

#[test]
fn an_offered_bin_votes_for_the_final_bin_holding_most_of_its_long_bases() {
    let lengths = [2000, 2000, 2000, 2000, 1000, 800, 600];
    let bin_of = HashMap::from([(0, 10), (1, 10), (2, 11), (3, 11)]);
    let sets = [
        set("a", &[&[0, 4], &[1, 2, 5], &[6]]),
        set("b", &[&[2, 3, 4]]),
    ];
    let order = [4, 5, 6];
    let ballot = |bins: &[usize], voters| Ballot {
        bins: bins.to_vec(),
        voters,
    };
    assert_eq!(
        ballots(&order, &lengths, 1500, &sets, &bin_of),
        vec![ballot(&[10, 11], 2), ballot(&[], 1), ballot(&[], 1)],
        "a split long half gives no vote, a bin of short contigs has no home, and b bins nothing \
         under 1000 bp so it has no say on the shorter two"
    );
    assert!(ballots(&order, &lengths, 1500, &[], &bin_of).is_empty());
}
