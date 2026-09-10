use rosella::markers::{ContigMarkers, Hit, MarkerRules, MarkerSet};
use rosella::quality::Scorer;

const TABLE: &str = "model_name\tdomain\tubiquity_percent\n\
                     alpha\tbac120\t100.0\n\
                     beta\tbac120\t50.0\n";

fn markers(per_contig: Vec<Vec<Hit>>, rules: MarkerRules) -> ContigMarkers {
    ContigMarkers::new(per_contig, MarkerSet::parse(TABLE), rules)
}

fn hit(marker: u16, partial: bool) -> Hit {
    Hit { marker, partial }
}

#[test]
fn a_gene_cut_by_two_contig_ends_is_one_marker_present_and_no_second_copy() {
    let split = vec![vec![hit(0, true)], vec![hit(0, true)]];
    let rules = MarkerRules {
        partial_counts: true,
        ..MarkerRules::default()
    };

    let shipped = markers(split.clone(), MarkerRules::default()).score(&[0, 1]);
    let counted = markers(split, rules).score(&[0, 1]);

    assert_eq!(shipped.completeness, 0.0);
    assert_eq!(counted.completeness, 50.0);
    assert_eq!(counted.contamination, 0.0);
}

#[test]
fn two_whole_copies_are_contamination_under_either_rule() {
    let doubled = vec![vec![hit(0, false), hit(0, false)]];
    let rules = MarkerRules {
        partial_counts: true,
        ..MarkerRules::default()
    };

    assert_eq!(
        markers(doubled.clone(), MarkerRules::default())
            .score(&[0])
            .contamination,
        50.0
    );
    assert_eq!(markers(doubled, rules).score(&[0]).contamination, 50.0);
}

#[test]
fn ubiquity_judges_presence_against_what_a_genome_is_expected_to_carry() {
    let one_of_two = vec![vec![hit(0, false)]];
    let rules = MarkerRules {
        ubiquity: true,
        ..MarkerRules::default()
    };

    assert_eq!(
        markers(one_of_two.clone(), MarkerRules::default())
            .score(&[0])
            .completeness,
        50.0
    );
    // The set expects 1.5 models, so holding the one every genome carries is two thirds whole.
    let weighted = markers(one_of_two, rules).score(&[0]).completeness;
    assert!((weighted - 100.0 / 1.5).abs() < 1e-9, "{weighted}");
}
