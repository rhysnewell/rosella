use rosella::markers::{ContigMarkers, Hit, MarkerRules, MarkerSet};
use rosella::quality::Scorer;

const WIDE: usize = 24;
const REDUCED: usize = 12;

fn table() -> String {
    let mut rows = String::from("model_name\tsets\tubiquity_wide\tubiquity_reduced\n");
    for marker in 0..WIDE {
        let in_reduced = marker < REDUCED;
        let sets = if in_reduced { "wide,reduced" } else { "wide" };
        let reduced_rate = if in_reduced { 0.98 } else { 0.50 };
        rows.push_str(&format!("m{marker}\t{sets}\t0.98\t{reduced_rate}\n"));
    }
    rows
}

fn sets_table() -> String {
    String::from(
        "set\tmedian_genome_bp\tmax_genome_bp\nwide\t3000000\t8000000\nreduced\t850000\t1800000\n",
    )
}

fn scored(present: &[usize]) -> (rosella::quality::Quality, String) {
    let hits = present
        .iter()
        .map(|marker| Hit {
            marker: *marker as u16,
            partial: false,
        })
        .collect::<Vec<_>>();
    let markers = ContigMarkers::new(
        vec![hits],
        MarkerSet::parse(&table()),
        MarkerRules::default(),
    );
    let held = markers.score(&[0]);
    let chosen = markers.set_name(held.set).to_string();
    (held, chosen)
}

#[test]
fn a_whole_reduced_genome_is_read_against_the_reduced_set() {
    let (held, chosen) = scored(&(0..REDUCED).collect::<Vec<_>>());

    assert_eq!(chosen, "reduced");
    assert_eq!(held.completeness, 100.0);
}

// A likelihood with no completeness term reads -31.9 for the reduced set against -47.2 for
// the wide one here, so it picks the reduced set. Fitting a share per set reverses that.
#[test]
fn a_half_built_wide_genome_is_not_mistaken_for_a_whole_reduced_one() {
    let every_other = (0..WIDE).step_by(2).collect::<Vec<_>>();
    let (held, chosen) = scored(&every_other);

    assert_eq!(chosen, "wide");
    assert_eq!(held.completeness, 50.0);
}

#[test]
fn too_few_markers_to_judge_falls_back_to_the_widest_set() {
    let (held, chosen) = scored(&[0, 1, 2]);

    assert_eq!(chosen, "wide");
    assert_eq!(held.completeness, 100.0 * 3.0 / WIDE as f64);
}

#[test]
fn a_bin_too_big_for_the_reduced_set_is_read_against_the_wide_one() {
    let hits = (0..REDUCED)
        .map(|marker| Hit {
            marker: marker as u16,
            partial: false,
        })
        .collect::<Vec<_>>();
    let markers = ContigMarkers::new(
        vec![hits],
        MarkerSet::parse(&table()).with_scales(&sets_table()),
        MarkerRules::default(),
    )
    .with_lengths(vec![4_000_000]);

    let held = markers.score(&[0]);

    assert_eq!(markers.set_name(held.set), "wide");
    assert_eq!(held.completeness, 50.0);
}
