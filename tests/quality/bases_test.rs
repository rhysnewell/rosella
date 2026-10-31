use std::collections::HashMap;

use rosella::quality::bases::{Bases, composition, n50};

#[test]
fn gc_spread_is_about_the_pooled_gc_over_contigs_longer_than_a_kilobase() {
    let contigs = [
        format!("{}{}", "G".repeat(1500), "A".repeat(500)),
        format!("{}NN", "a".repeat(2000)),
        "G".repeat(1000),
    ];
    let bases = contigs
        .iter()
        .enumerate()
        .map(|(at, contig)| (at, Bases::count(contig.as_bytes())))
        .collect::<HashMap<_, _>>();
    let lengths = contigs.iter().map(String::len).collect::<Vec<_>>();

    let held = composition(&bases, &lengths, &[0, 1, 2]).unwrap();

    assert!((held.gc - 50.0).abs() < 1e-9);
    assert!((held.gc_spread - 100.0 * 0.15625f64.sqrt()).abs() < 1e-9);
    assert_eq!(held.ambiguous, 2);
    assert!(composition(&bases, &lengths, &[0, 3]).is_none());
}

#[test]
fn n50_is_the_contig_that_carries_the_running_total_past_half() {
    assert_eq!(n50(&mut [100, 400, 200, 300]), 300);
    assert_eq!(n50(&mut [200, 500, 300]), 500);
}
