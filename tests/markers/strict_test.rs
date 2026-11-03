use rosella::markers::checkm::{Copies, Panel};
use rosella::markers::replicon::Shape;
use rosella::markers::{ContigMarkers, Hit, MarkerSet, Strict};
use rosella::quality::Scorer;

struct Bac {
    gtdb: u16,
    copy: Copies,
}

fn bac() -> Bac {
    let set = MarkerSet::embedded();
    let panel = Panel::embedded(|model| set.id(model));
    let at = panel.lineage("bac").unwrap();
    Bac {
        gtdb: set.id("ATPsyn_F1gamma").unwrap(),
        copy: Copies {
            set: at as u8,
            model: panel.models(at)[0],
            copies: 1,
        },
    }
}

fn hit(marker: u16) -> Hit {
    Hit {
        marker,
        ..Default::default()
    }
}

fn markers(hits: Vec<Vec<Hit>>, checkm: Vec<Vec<Copies>>) -> ContigMarkers {
    let lengths = vec![900_000; hits.len()];
    ContigMarkers::new(hits, MarkerSet::embedded())
        .with_lengths(lengths)
        .with_checkm(checkm.into_iter().map(Some).collect())
}

#[test]
fn the_stricter_panel_sets_the_contamination() {
    let Bac { gtdb, copy } = bac();
    let checkm_doubled = markers(
        vec![vec![hit(gtdb)], Vec::new()],
        vec![vec![copy], vec![copy]],
    );
    let gtdb_doubled = markers(
        vec![vec![hit(gtdb)], vec![hit(gtdb)]],
        vec![vec![copy], vec![]],
    );

    let strict = |held: &ContigMarkers| {
        Strict {
            inner: held,
            markers: held,
            fence: None,
            markerless: false,
        }
        .score(&[0, 1])
        .contamination
    };
    let set = checkm_doubled.score(&[0, 1]).set;
    let checkm = checkm_doubled.checkm_contamination(&[0, 1], set).unwrap();

    assert_eq!(checkm_doubled.score(&[0, 1]).contamination, 0.0);
    assert!(checkm > 0.0);
    assert_eq!(strict(&checkm_doubled), checkm);
    assert_eq!(
        strict(&gtdb_doubled),
        gtdb_doubled.score(&[0, 1]).contamination
    );
}

#[test]
fn a_contig_holding_only_a_checkm_copy_the_bin_has_leaves_under_checkm() {
    let Bac { gtdb, copy } = bac();
    let held = markers(
        vec![vec![hit(gtdb)], Vec::new()],
        vec![vec![copy], vec![copy]],
    );

    assert!(held.redundant_on(&[0, 1], false).is_empty());
    assert_eq!(held.redundant_on(&[0, 1], true), vec![1]);
}

#[test]
fn a_bin_coding_like_a_eukaryote_keeps_the_gtdb_contamination() {
    let Bac { gtdb, copy } = bac();
    let coding = |coding_bases| {
        markers(
            vec![vec![hit(gtdb)], Vec::new()],
            vec![vec![copy], vec![copy]],
        )
        .with_shapes(vec![
            Some(Shape {
                coding_bases,
                genes: 1,
            });
            2
        ])
    };
    let contamination = |held: &ContigMarkers| {
        Strict {
            inner: held,
            markers: held,
            fence: Some(0.8),
            markerless: false,
        }
        .score(&[0, 1])
        .contamination
    };

    assert!(contamination(&coding(810_000)) > 0.0);
    assert_eq!(contamination(&coding(495_000)), 0.0);
}

#[test]
fn a_checkm_model_doubled_in_every_bin_counts_once() {
    let Bac { copy, .. } = bac();
    let set = MarkerSet::embedded();
    let whole = (0..set.len() as u16).map(hit).collect::<Vec<_>>();
    let panel = Panel::embedded(|model| set.id(model));
    let every = panel
        .models(usize::from(copy.set))
        .into_iter()
        .map(|model| Copies { model, ..copy })
        .collect::<Vec<_>>();
    let bins = 8;
    let mut hits = Vec::new();
    let mut checkm = Vec::new();
    for _ in 0..bins {
        hits.extend([whole.clone(), Vec::new()]);
        checkm.extend([every.clone(), vec![copy]]);
    }
    let held = markers(hits, checkm);
    let strict = || {
        Strict {
            inner: &held,
            markers: &held,
            fence: None,
            markerless: false,
        }
        .score(&[0, 1])
        .contamination
    };

    assert!(strict() > held.score(&[0, 1]).contamination);
    held.learn_multicopy(
        (0..bins)
            .map(|bin| [2 * bin, 2 * bin + 1])
            .collect::<Vec<_>>()
            .iter()
            .map(|b| b.as_slice()),
    );
    assert_eq!(strict(), held.score(&[0, 1]).contamination);
    assert!(held.redundant_on(&[0, 1], true).is_empty());
}
