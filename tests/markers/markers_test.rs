use rosella::external::hmmer_engine::best_hits;
use rosella::markers::orfs::translate;
use rosella::markers::{ContigMarkers, Hit, MarkerSet};

const TABLE: &str = "accession\tmodel_name\tmarker_id\tdomain\n\
PF1\tS9\tPF1\tbac120\n\
PF2\tS8\tPF2\tbac120,ar53\n\
PF3\tL3\tPF3\tar53\n";

fn hit(marker: u16, partial: bool) -> Hit {
    Hit { marker, partial }
}

fn bundle() -> ContigMarkers {
    let set = MarkerSet::parse(TABLE);
    ContigMarkers::new(
        vec![
            vec![hit(0, false), hit(1, false)],
            vec![hit(0, false)],
            vec![hit(1, true)],
            vec![],
        ],
        set,
    )
}

#[test]
fn a_second_copy_reads_as_duplication_only_when_both_are_whole() {
    let markers = bundle();
    let fused = markers.fusion(&[0, 1]).expect("two carriers");
    assert_eq!((fused.present, fused.duplicated), (2, 1));
    let whole = markers.fusion(&[0, 2]).expect("one carrier");
    assert_eq!((whole.present, whole.duplicated), (2, 0));
    assert!(markers.fusion(&[3]).is_none());
    assert_eq!(markers.distinct(&[0, 1, 2]), 2);
}

#[test]
fn the_domain_with_more_markers_is_the_one_scored() {
    let set = MarkerSet::parse(TABLE);
    let markers = ContigMarkers::new(vec![vec![hit(1, false), hit(2, false)]], set);
    let read = markers.fusion(&[0]).expect("an archaeon");
    assert_eq!(read.present, 2);
}

#[test]
fn the_best_scoring_model_owns_a_protein() {
    let table = "# comment\n\
p1 - S9 PF1 1e-30 120.0 0.1\n\
p1 - S8 PF2 1e-10 40.0 0.1\n\
p2 - L3 PF3 1e-5 20.5 0.0\n";
    let hits = best_hits(table);
    assert_eq!(hits["p1"].0, "S9");
    assert_eq!(hits["p2"], ("L3".to_string(), 20.5));
}

#[test]
fn translation_rewrites_the_start_and_drops_the_stop() {
    assert_eq!(translate(b"GTGAAATAA", true), "MK");
    assert_eq!(translate(b"GTGAAATAA", false), "VK");
    assert_eq!(translate(b"ATGNNNTAA", true), "MX");
}
