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

fn carriers(copy: Copies) -> Vec<Vec<u32>> {
    ContigMarkers::new(vec![Vec::new(), Vec::new()], MarkerSet::embedded())
        .with_checkm(vec![vec![copy], vec![copy]])
        .carriers(&[0, 1])
}

#[test]
fn a_family_gtdb_lacks_joins_the_steps() {
    assert_eq!(carriers(bac_model(false)), vec![vec![0, 1]]);
}

#[test]
fn a_family_gtdb_holds_adds_no_second_carrier() {
    assert!(carriers(bac_model(true)).is_empty());
}
