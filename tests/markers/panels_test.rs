use std::fmt::Write;
use std::path::Path;

use rosella::markers::panels::{LINEAGES, Panel};
use rosella::markers::{ContigMarkers, Hit, MarkerSet};
use rosella::quality::{Bars, Scorer};

const PROKARYOTE: usize = 12;
const FUNGI: usize = 30;
const CUTOFF: f64 = 100.0;

fn model(lineage: &str, at: usize) -> String {
    format!("{at}at{lineage}")
}

fn panel(directory: &Path, chosen: &[&'static str]) -> Panel {
    for lineage in LINEAGES {
        let folder = directory.join(format!("{lineage}_odb12"));
        std::fs::create_dir_all(folder.join("hmms")).unwrap();
        let mut cutoffs = String::new();
        for at in 0..FUNGI {
            let name = model(lineage, at);
            writeln!(cutoffs, "{name}\t{CUTOFF}").unwrap();
            let text = format!("HMMER3/f\nNAME  {at}\nLENG  10\n//\n");
            std::fs::write(folder.join("hmms").join(format!("{name}.hmm")), text).unwrap();
        }
        std::fs::write(folder.join("scores_cutoff"), cutoffs).unwrap();
    }
    Panel::prepare(directory, chosen, false).unwrap()
}

fn table() -> String {
    let mut rows = String::from("model_name\tsets\n");
    for at in 0..PROKARYOTE {
        writeln!(rows, "p{at}\tbac").unwrap();
    }
    rows
}

fn decoded(set: &MarkerSet, name: &str, score: f64) -> Hit {
    stretch(set, name, score, 1, 100)
}

fn stretch(set: &MarkerSet, name: &str, score: f64, from: u32, to: u32) -> Hit {
    let field = format!("{name}:{from}:{to}:100:1:100:{score}:1:300:+:0:0");
    Hit::decode(&field, set).unwrap()
}

#[test]
fn each_chosen_lineage_is_joined_alone_under_its_busco_ids() {
    let directory = tempfile::tempdir().unwrap();
    let panel = panel(directory.path(), &["fungi", "metazoa"]);
    let names = panel.lineages.iter().map(|lineage| lineage.name);
    assert_eq!(names.collect::<Vec<_>>(), ["fungi", "metazoa"]);
    for lineage in ["fungi", "metazoa"] {
        let joined = std::fs::read_to_string(directory.path().join(format!("{lineage}.hmm")));
        let joined = joined.unwrap();
        let names = joined
            .lines()
            .filter_map(|line| line.strip_prefix("NAME"))
            .map(str::trim)
            .collect::<Vec<_>>();
        let expected = (0..FUNGI).map(|at| model(lineage, at)).collect::<Vec<_>>();
        assert_eq!(names, expected);
    }
    assert!(!directory.path().join("eukaryota.hmm").exists());
}

#[test]
fn a_fungus_reads_on_its_own_set_and_a_sparse_bin_stays_prokaryotic() {
    let directory = tempfile::tempdir().unwrap();
    let panel = panel(directory.path(), &["fungi"]);
    let set = MarkerSet::parse_with(&table(), Some(&panel));
    let fungus = (0..FUNGI)
        .map(|at| decoded(&set, &model("fungi", at), 2.0 * CUTOFF))
        .collect::<Vec<_>>();
    let sparse = vec![decoded(&set, "p0", 50.0), decoded(&set, "p1", 50.0)];
    let markers = ContigMarkers::new(vec![fungus, sparse], set).with_lengths(vec![1000; 2]);

    let held = markers.score(&[0]);
    assert_eq!(markers.set_name(held.set), "fungi");
    assert_eq!((held.completeness, held.contamination), (100.0, 0.0));
    let held = markers.score(&[1]);
    assert_eq!(
        markers.set_name(held.set),
        "bac",
        "the widest set is a prokaryote's however many models a panel brings"
    );
    let bars = Bars {
        completeness: 90.0,
        contamination: 5.0,
    };
    assert_eq!(
        markers.eukaryotic([&[0][..], &[1][..]], bars),
        [true, false]
    );
    assert!(
        !markers.carries(0),
        "a panel hit makes no prokaryote carrier"
    );
}

#[test]
fn only_a_whole_call_over_the_same_stretch_is_a_second_copy() {
    let directory = tempfile::tempdir().unwrap();
    let panel = panel(directory.path(), &["fungi"]);
    let set = MarkerSet::parse_with(&table(), Some(&panel));
    let mut fungus = (3..FUNGI)
        .map(|at| decoded(&set, &model("fungi", at), 2.0 * CUTOFF))
        .collect::<Vec<_>>();
    let short = decoded(&set, &model("fungi", 0), 0.5 * CUTOFF);
    assert!(short.partial);
    fungus.extend([
        decoded(&set, &model("fungi", 0), 2.0 * CUTOFF),
        short,
        decoded(&set, &model("fungi", 1), 2.0 * CUTOFF),
        decoded(&set, &model("fungi", 1), 2.0 * CUTOFF),
        stretch(&set, &model("fungi", 2), 2.0 * CUTOFF, 1, 50),
        stretch(&set, &model("fungi", 2), 2.0 * CUTOFF, 51, 100),
    ]);
    let markers = ContigMarkers::new(vec![fungus], set).with_lengths(vec![1000]);

    let held = markers.score(&[0]);
    assert_eq!(held.completeness, 100.0);
    assert_eq!(
        held.contamination,
        100.0 / FUNGI as f64,
        "two exons of one gene are one copy and two whole genes are two"
    );
}
