use rosella::markers::hmm_table::Reach;
use rosella::markers::{ContigMarkers, Hit, MarkerSet, Place};
use rosella::quality::Scorer;

const TABLE: &str = "model_name\tdomain\n\
                     alpha\tbac120\n\
                     beta\tbac120\n";

fn markers(per_contig: Vec<Vec<Hit>>) -> ContigMarkers {
    ContigMarkers::new(per_contig, MarkerSet::parse(TABLE))
}

fn hit(marker: u16, partial: bool) -> Hit {
    Hit {
        marker,
        partial,
        ..Default::default()
    }
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

    assert_eq!(markers(doubled).score(&[0]).contamination, 50.0);
}

#[test]
fn half_a_gene_beside_its_other_half_is_no_repeat_but_a_second_copy_is() {
    let over = |from, to| Hit {
        marker: 0,
        partial: true,
        place: Place {
            reach: Reach {
                model_from: from,
                model_to: to,
                model_length: 200,
                ..Default::default()
            },
            ..Default::default()
        },
    };
    let bin = markers(vec![
        vec![over(1, 100)],
        vec![over(101, 200)],
        vec![over(20, 120)],
        vec![hit(1, false)],
    ]);

    assert_eq!(bin.repeats_any(&[0, 1, 3], 1), Some(true));
    assert_eq!(bin.repeats_in_place(&[0, 1, 3], 1), Some(false));
    assert_eq!(bin.repeats_in_place(&[0, 2, 3], 2), Some(true));
}

#[test]
fn a_bin_holding_an_unsearched_contig_has_no_checkm_reading() {
    let scorer = ContigMarkers::new(vec![vec![hit(0, false)], Vec::new()], MarkerSet::embedded())
        .with_checkm(vec![Some(Vec::new()), None]);

    assert!(scorer.checkm(&[0]).is_some());
    assert!(scorer.checkm(&[0, 1]).is_none());
}

#[test]
fn a_cut_copy_fills_the_one_copy_column_and_only_whole_copies_reach_two() {
    let bin = markers(vec![
        vec![hit(0, false), hit(0, false)],
        vec![hit(1, true)],
        vec![hit(0, true)],
    ]);
    let (chosen, sets) = bin.readings(&[0, 1, 2]).unwrap();
    let gtdb = sets[chosen].gtdb;

    assert_eq!(gtdb.copies, [0, 1, 1, 0, 0, 0]);
    assert_eq!(gtdb.contamination, 50.0);
}

#[test]
fn a_duplicate_names_the_contigs_holding_its_whole_copies() {
    let bin = markers(vec![
        vec![hit(0, false)],
        vec![hit(0, false), hit(1, false)],
        vec![hit(0, true)],
    ]);
    let rows = bin
        .duplicates(&[0, 1, 2])
        .into_iter()
        .map(|row| (row.panel, row.marker, row.contig, row.copies))
        .collect::<Vec<_>>();

    assert_eq!(
        rows,
        vec![
            ("gtdb", "alpha".to_string(), 0, 1),
            ("gtdb", "alpha".to_string(), 1, 1),
        ]
    );
}
