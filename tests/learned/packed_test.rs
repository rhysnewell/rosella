use std::collections::HashSet;

use rosella::learned::{packed, tied_families};

// Carriers shared between families, and a family past the step cap whose last chunk is one
// contig with no family mate to be told apart from.
#[test]
fn packed_steps_hold_every_carrier_once_and_keep_small_families_whole() {
    let per_step = 8;
    let families = vec![
        vec![0, 1, 2, 3],
        vec![2, 3, 4, 5],
        vec![6, 7, 8],
        (10..27).collect(),
    ];
    let steps = packed(&families, per_step);
    for step in &steps {
        let unique = step.iter().collect::<HashSet<_>>();
        assert_eq!(unique.len(), step.len(), "a contig repeats in {step:?}");
        assert!(step.len() <= per_step && step.len() > 1, "{step:?}");
    }
    let held = steps.iter().flatten().copied().collect::<HashSet<u32>>();
    let carriers = families[..3]
        .iter()
        .flatten()
        .copied()
        .chain(10..26)
        .collect::<HashSet<_>>();
    assert_eq!(held, carriers);
    for family in &families[..3] {
        assert!(
            steps
                .iter()
                .any(|step| family.iter().all(|contig| step.contains(contig))),
            "{family:?} is split across {steps:?}"
        );
    }
}

// A carrier under the gate leaves, and a family left with one carrier is dropped.
#[test]
fn tied_families_keep_carriers_through_the_gate() {
    let families = vec![vec![0, 1, 2], vec![3, 4], vec![1, 5, 6]];
    let tied = HashSet::from([0, 1, 3, 5, 6]);
    let kept = tied_families(&families, |carrier| tied.contains(&carrier));
    assert_eq!(kept, vec![vec![0, 1], vec![1, 5, 6]]);
}
