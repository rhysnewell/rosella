use rosella::markers::{ContigMarkers, Hit, MarkerSet};
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
fn two_whole_copies_are_contamination() {
    let doubled = vec![vec![hit(0, false), hit(0, false)]];

    assert_eq!(
        markers(doubled)
            .score(&[0])
            .contamination,
        50.0
    );
}
