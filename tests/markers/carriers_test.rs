use rosella::markers::checkm::{Copies, Panel};
use rosella::markers::{ContigMarkers, MarkerSet};

fn bac_model(shared: bool) -> Copies {
    let set = MarkerSet::embedded();
    let panel = Panel::embedded(|model| set.id(model));
    let at = panel.lineage("bac").unwrap();
    let model = *panel
        .models(at)
        .iter()
        .find(|model| {
            let gtdb = panel
                .searched_name(**model)
                .is_some_and(|name| set.id(name).is_some())
                || set.holds_family(panel.family(**model));
            gtdb == shared
        })
        .unwrap();
    Copies {
        set: at as u8,
        model,
        score: 100.0,
    }
}

fn carriers(lacking: usize, holding: usize) -> Vec<Vec<u32>> {
    let contigs = lacking.max(holding);
    let copies = (0..contigs)
        .map(|at| {
            [(at < lacking, false), (at < holding, true)]
                .into_iter()
                .filter(|(carried, _)| *carried)
                .map(|(_, shared)| bac_model(shared))
                .collect()
        })
        .collect();
    ContigMarkers::new(vec![Vec::new(); contigs], MarkerSet::embedded())
        .with_checkm(copies)
        .carriers(&(0..contigs).collect::<Vec<_>>())
}

#[test]
fn a_family_gtdb_lacks_joins_the_steps() {
    assert_eq!(carriers(2, 2), vec![vec![0, 1]]);
}

#[test]
fn a_family_gtdb_holds_adds_no_second_carrier() {
    assert!(carriers(0, 2).is_empty());
}

#[test]
fn a_family_on_more_contigs_than_any_gtdb_holds_stays_out() {
    assert!(carriers(3, 2).is_empty());
}
