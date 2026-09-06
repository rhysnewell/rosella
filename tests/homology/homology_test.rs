use rosella::external::skani_engine::parse_table;
use rosella::homology::{Homology, HomologySettings, Pair};

fn names() -> Vec<String> {
    ["c0", "c1", "c2", "c3"]
        .iter()
        .map(|name| (*name).to_string())
        .collect()
}

fn pair(one: usize, other: usize, identity: f64, aligned: f64) -> Pair {
    Pair {
        one,
        other,
        identity,
        aligned_one: aligned,
        aligned_other: aligned,
    }
}

#[test]
fn the_log_skani_writes_before_the_table_is_not_the_header() {
    let table = "\
Sketching files...
Wrote 4 sketches.
Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name
x.fna\tx.fna\t98.7\t71.2\t64.0\tc0 length=100\tc2 length=90
x.fna\tx.fna\t91.0\t12.0\t80.0\tc1\tc3
";
    let pairs = parse_table(table, &names());
    assert_eq!(pairs.len(), 2);
    assert_eq!((pairs[0].one, pairs[0].other), (0, 2));
    assert!((pairs[0].identity - 98.7).abs() < 1e-9);
    assert!((pairs[1].aligned_fraction() - 12.0).abs() < 1e-9);
}

#[test]
fn a_pair_covered_on_one_side_only_is_a_repeat_not_an_organism() {
    let lengths = [400_000, 3_000, 400_000, 3_000];
    let settings = HomologySettings::default();
    let repeat = Homology::from_pairs([pair(0, 2, 99.0, 4.0)], &lengths, settings);
    assert!(!repeat.cannot_link(0, 2));

    let both_ends = Homology::from_pairs([pair(0, 2, 99.0, 60.0)], &lengths, settings);
    assert!(both_ends.cannot_link(0, 2) && both_ends.cannot_link(2, 0));
}

#[test]
fn two_short_contigs_are_held_back_by_the_pair_length_bar() {
    let lengths = [400_000, 3_000, 4_000, 3_000];
    let settings = HomologySettings {
        min_pair_length: 20_000,
        ..HomologySettings::default()
    };
    let table = Homology::from_pairs(
        [pair(1, 3, 99.0, 90.0), pair(0, 2, 99.0, 90.0)],
        &lengths,
        settings,
    );
    assert!(!table.cannot_link(1, 3));
    assert!(table.cannot_link(0, 2));
}

/// The behaviour solo gets wrong: two genome-sized contigs that do not align are one genome.
#[test]
fn contigs_that_never_align_stay_in_one_group() {
    let lengths = [3_000_000, 2_500_000, 2_800_000, 1_000];
    let settings = HomologySettings::default();

    let silent = Homology::from_pairs([], &lengths, settings);
    assert_eq!(silent.groups(&[0, 1, 2]), vec![vec![0, 1, 2]]);

    let one_pair = Homology::from_pairs([pair(0, 1, 99.0, 80.0)], &lengths, settings);
    assert_eq!(one_pair.groups(&[0, 1, 2]), vec![vec![0, 2], vec![1]]);

    let all = Homology::from_pairs(
        [
            pair(0, 1, 99.0, 80.0),
            pair(0, 2, 99.0, 80.0),
            pair(1, 2, 99.0, 80.0),
        ],
        &lengths,
        settings,
    );
    assert_eq!(all.groups(&[0, 1, 2]), vec![vec![0], vec![1], vec![2]]);
}

/// A bin holding one side of a homologous pair proves nothing. Both sides prove two organisms.
#[test]
fn only_both_sides_in_one_bin_are_evidence() {
    let lengths = [400_000, 400_000, 400_000, 400_000];
    let table = Homology::from_pairs(
        [pair(0, 2, 98.0, 70.0)],
        &lengths,
        HomologySettings::default(),
    );
    assert!(table.holds_pair(&[3, 0, 1, 2]));
    assert!(!table.holds_pair(&[0, 1, 3]));
    assert!(!table.holds_pair(&[1, 3]));
}
