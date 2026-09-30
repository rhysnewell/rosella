use std::cmp::Ordering::{Equal, Greater, Less};
use std::collections::HashMap;

use rosella::refine::cut_report::{CutLog, owners};

use crate::scorer::MarkerScorer;

const PIECE: usize = 100_000;

// 0 to 3 hold 16 markers. 4 fills a gap, 5 repeats 0 to 3 and 6 has none.
fn markers() -> Vec<Vec<usize>> {
    vec![
        (0..4).collect(),
        (4..8).collect(),
        (8..12).collect(),
        (12..16).collect(),
        vec![16],
        (0..4).collect(),
        vec![],
    ]
}

#[test]
fn a_cut_is_priced_against_the_final_bin_its_old_bin_became() {
    let lengths = vec![PIECE; markers().len()];
    let mut log = CutLog::default();
    let pieces = owners([(0, &vec![0, 1, 2, 3]), (1, &vec![4])]);
    log.cut("split", &[0, 1, 2, 3, 4, 5, 6], &pieces, |contig| {
        lengths[contig]
    });
    log.cut("shed", &[5], &HashMap::new(), |contig| lengths[contig]);
    let finals = HashMap::from([(7, vec![0, 1, 2, 3]), (8, vec![4])]);
    let names = (0..lengths.len())
        .map(|contig| format!("c{contig}"))
        .collect::<Vec<_>>();
    let path = tempfile::NamedTempFile::new().unwrap();

    log.write(
        path.path(),
        &finals,
        &MarkerScorer::new(markers(), 20),
        2.0,
        &lengths,
        &names,
    )
    .unwrap();

    let rows = std::fs::read_to_string(path.path())
        .unwrap()
        .lines()
        .skip(1)
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            let gain = fields[12]
                .parse::<f64>()
                .ok()
                .and_then(|gain| gain.partial_cmp(&0.0));
            (
                fields[0].to_string(),
                fields[1].to_string(),
                fields[5].to_string(),
                gain,
            )
        })
        .collect::<Vec<_>>();
    let row = |stage: &str, contig: &str, target: &str, gain| {
        (
            stage.to_string(),
            contig.to_string(),
            target.to_string(),
            gain,
        )
    };
    assert_eq!(
        rows,
        vec![
            row("split", "c4", "7", Some(Greater)),
            row("split", "c5", "7", Some(Less)),
            row("split", "c6", "7", Some(Equal)),
            row("shed", "c5", "-", None),
        ]
    );
}
