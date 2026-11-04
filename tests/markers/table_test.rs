use rosella::markers::checkm::Panel;
use rosella::markers::{ContigMarkers, Hit, MarkerSet};
use rosella::quality::Scorer;

fn family(set: &MarkerSet, panel: &Panel, model: &str) -> Option<u16> {
    set.family(panel.id(model).unwrap())
}

#[test]
fn a_checkm_model_of_a_gtdb_family_is_that_gtdb_marker_and_a_cpr_twin_shares_its_family() {
    let set = MarkerSet::embedded();
    let panel = Panel::embedded(|model| set.id(model));
    for shared in ["TIGR00967", "PF00410.14", "cpr_PF00410.20"] {
        assert_eq!(family(&set, &panel, shared), None, "{shared}");
    }
    let own = family(&set, &panel, "PF00276.15").unwrap();
    assert!(set.from_checkm(own));
    assert_eq!(set.name(own), "checkm_PF00276");
    assert_eq!(family(&set, &panel, "cpr_PF00276.21"), Some(own));
}

#[test]
fn a_set_reads_a_checkm_family_only_where_checkm_holds_it() {
    let set = MarkerSet::embedded();
    let panel = Panel::embedded(|model| set.id(model));
    let archaeal = family(&set, &panel, "PF00832.15").unwrap();
    let hit = Hit {
        marker: archaeal,
        ..Default::default()
    };
    let markers = ContigMarkers::new(vec![vec![hit]], set);
    let held = |name: &str| {
        let at = (0..4u16).find(|at| markers.set_name(*at) == name).unwrap();
        markers
            .points_on(&[0], 1.0, at)
            .iter()
            .any(|(marker, _)| *marker == usize::from(archaeal))
    };
    assert!(held("ar"));
    assert!(!held("bac"));
    assert!(!held("dpann"));
}
