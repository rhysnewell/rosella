use rosella::markers::checkm::{Copies, Panel};
use rosella::markers::{ContigMarkers, Hit, MarkerSet};
use rosella::quality::Scorer;

const TABLE: &str = include_str!("../../data/gtdb_markers.tsv");

fn gtdb(shared: bool, model: u16, panel: &Panel, set: &MarkerSet) -> bool {
    let gtdb = panel
        .searched_name(model)
        .is_some_and(|name| set.id(name).is_some())
        || set.holds_family(panel.family(model));
    gtdb == shared
}

// A family of one lineage alone, so a contig placed in the other cannot read it.
fn model(lineage: &str, other: &str, shared: bool) -> Copies {
    let set = MarkerSet::joined(None);
    let panel = Panel::embedded(|model| set.id(model));
    let (at, theirs) = (
        panel.lineage(lineage).unwrap(),
        panel.lineage(other).unwrap(),
    );
    let model = *panel
        .models(at)
        .iter()
        .find(|model| {
            let family = panel.family(**model);
            gtdb(shared, **model, &panel, &set)
                && panel
                    .models(theirs)
                    .iter()
                    .all(|held| panel.family(*held) != family)
        })
        .unwrap();
    Copies {
        set: at as u8,
        model,
        score: 100.0,
        partial: false,
    }
}

// Cut hits place a contig without making it a GTDB carrier.
fn placing(sets: &[&str], take: usize) -> Vec<Hit> {
    let set = MarkerSet::joined(None);
    let mut hits = TABLE
        .lines()
        .skip(1)
        .filter_map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            sets.contains(&fields[10]).then(|| set.id(fields[1]))?
        })
        .map(|marker| Hit {
            marker,
            partial: true,
            ..Default::default()
        })
        .take(take)
        .collect::<Vec<_>>();
    hits.sort_unstable_by_key(|hit| hit.marker);
    hits
}

fn carriers(hits: Vec<Hit>, copy: Copies) -> (Vec<Vec<u32>>, String) {
    let held =
        ContigMarkers::new(vec![hits; 2], MarkerSet::joined(None)).with_checkm(vec![vec![copy]; 2]);
    let placed = held.set_name(held.score(&[0]).set).to_string();
    (held.carriers(&[0, 1]), placed)
}

#[test]
fn a_contig_reads_only_whole_checkm_copies_of_the_lineage_its_markers_fit() {
    let bacterial = || placing(&["bac"], usize::MAX);
    let archaeal = || placing(&["ar", "ar,dpann"], usize::MAX);
    let bac = model("bac", "ar", false);
    let ar = model("ar", "bac", false);
    let cases = [
        (bacterial(), bac, "bac", true),
        (archaeal(), ar, "ar", true),
        (bacterial(), ar, "bac", false),
        (archaeal(), bac, "ar", false),
        (placing(&["bac"], 3), bac, "bac", false),
        (
            bacterial(),
            Copies {
                partial: true,
                ..bac
            },
            "bac",
            false,
        ),
    ];
    for (at, (hits, copy, lineage, joins)) in cases.into_iter().enumerate() {
        let (found, placed) = carriers(hits, copy);

        assert_eq!(placed, lineage, "case {at}");
        assert_eq!(found == vec![vec![0, 1]], joins, "case {at}");
    }
}

#[test]
fn a_family_gtdb_holds_adds_no_second_carrier() {
    let (found, _) = carriers(placing(&["bac"], usize::MAX), model("bac", "ar", true));

    assert!(found.is_empty());
}
