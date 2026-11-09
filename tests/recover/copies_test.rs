use std::collections::HashMap;

use rosella::recover::copies::Copies;

#[test]
fn copies_need_shared_reads_and_one_length() {
    let lengths = vec![
        2000, 2010, 2020, 2005, 1200, 2003, 5000, 3000, 1750, 1752, 1749,
    ];
    let reads = (0..lengths.len())
        .map(|contig| (contig, 100.0))
        .collect::<HashMap<_, _>>();
    let mut ties = HashMap::new();
    let mut tie = |first: usize, second: usize| {
        ties.insert((first.min(second), first.max(second)), 900);
    };
    for first in 0..4 {
        for second in first + 1..5 {
            tie(first, second);
        }
    }
    tie(6, 7);
    tie(8, 9);
    tie(9, 10);
    let copies = Copies {
        lengths: &lengths,
        reads: &reads,
        ties: &ties,
        read_length: 150.0,
    };
    assert_eq!(
        copies.of(&(0..lengths.len()).collect::<Vec<_>>()),
        vec![vec![0, 1, 2, 3, 4]],
        "contig 5 shares no reads, 6 and 7 differ in length, 8 and 10 share none"
    );
}
