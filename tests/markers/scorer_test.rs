use rosella::markers::{ContigMarkers, Hit, MarkerSet, Partials};
use rosella::quality::Scorer;

const TABLE: &str = "model_name\tdomain\n\
                     alpha\tbac120\n\
                     beta\tbac120\n";

fn markers(per_contig: Vec<Vec<Hit>>) -> ContigMarkers {
    ContigMarkers::new(per_contig, MarkerSet::parse(TABLE))
}

fn hit(marker: u16, partial: bool) -> Hit {
    Hit { marker, partial }
}

#[test]
fn a_gene_cut_by_two_contig_ends_is_one_marker_present_and_no_second_copy() {
    let split = vec![vec![hit(0, true)], vec![hit(0, true)]];
    let scored = markers(split).score(&[0, 1]);

    assert_eq!(scored.completeness, 50.0);
    assert_eq!(scored.contamination, 0.0);
}

#[test]
fn counting_partials_charges_a_fragment_beside_a_whole_copy() {
    let whole_and_fragment = vec![vec![hit(0, false)], vec![hit(0, true)]];
    let counted = markers(whole_and_fragment.clone())
        .with_partials(Partials::Count)
        .score(&[0, 1]);
    let ignored = markers(whole_and_fragment).score(&[0, 1]);

    assert_eq!(ignored.contamination, 0.0);
    assert_eq!(counted.contamination, 50.0);
    assert_eq!(counted.completeness, ignored.completeness);
}

#[test]
fn two_whole_copies_are_contamination() {
    let doubled = vec![vec![hit(0, false), hit(0, false)]];

    assert_eq!(markers(doubled).score(&[0]).contamination, 50.0);
}

#[test]
fn a_duplicate_counts_by_how_single_copy_the_marker_is() {
    const RATED: &str = "model_name\tdomain\tsets\tubiquity_bac\tsingle_copy_bac\n\
                         alpha\tbac120\tbac\t0.99\t0.40\n\
                         beta\tbac120\tbac\t0.99\t1.00\n";
    let rated = |at: u16| {
        ContigMarkers::new(
            vec![vec![hit(at, false), hit(at, false)]],
            MarkerSet::parse(RATED),
        )
        .score(&[0])
        .contamination
    };

    assert!(
        (rated(0) - 20.0).abs() < 1e-9,
        "a marker single copy in 40 per cent of genomes carries 0.4 of a duplicate, \
         got {}",
        rated(0)
    );
    assert!(
        (rated(1) - 50.0).abs() < 1e-9,
        "a marker single copy everywhere carries a whole duplicate, got {}",
        rated(1)
    );
}

#[test]
fn counting_carriers_charges_a_second_copy_only_when_another_contig_holds_it() {
    use rosella::markers::Duplicates;

    let together = markers(vec![vec![hit(0, false), hit(0, false)]]).counting(Duplicates::Carriers);
    assert_eq!(together.score(&[0]).contamination, 0.0);
    assert_eq!(together.score(&[0]).completeness, 50.0);

    let apart =
        markers(vec![vec![hit(0, false)], vec![hit(0, false)]]).counting(Duplicates::Carriers);
    assert_eq!(apart.score(&[0, 1]).contamination, 50.0);
}
