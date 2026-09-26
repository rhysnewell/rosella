use rosella::markers::fragments::{Bar, Bars, accepted, complete, floor, gathering};

fn row(protein: &str, model: &str, length: u32, score: f64, from: u32, to: u32) -> String {
    let mut fields = vec!["-"; 23];
    let score = score.to_string();
    let length = length.to_string();
    let from = from.to_string();
    let to = to.to_string();
    fields[0] = protein;
    fields[3] = model;
    fields[5] = &length;
    fields[7] = &score;
    fields[13] = &score;
    fields[15] = &from;
    fields[16] = &to;
    fields.join(" ")
}

fn bars() -> Bars {
    Bars::from([(
        "Ribosomal_S9".to_string(),
        Bar {
            sequence: 100.0,
            domain: 100.0,
        },
    )])
}

#[test]
fn a_half_model_hit_needs_half_the_gathering_score() {
    let over = row("7", "Ribosomal_S9", 200, 51.0, 1, 100);
    let under = row("8", "Ribosomal_S9", 200, 49.0, 1, 100);
    let taken = accepted(&[over, under].join("\n"), &bars(), 0.05, |_| true);

    assert_eq!(taken.len(), 1);
    assert!(taken.contains_key(&7));
}

#[test]
fn a_sliver_of_a_model_is_refused_however_it_scores() {
    let sliver = row("9", "Ribosomal_S9", 200, 900.0, 1, 4);
    assert!(accepted(&sliver, &bars(), 0.05, |_| true).is_empty());
}

#[test]
fn a_model_the_table_does_not_carry_is_refused() {
    let unknown = row("9", "Nothing", 200, 900.0, 1, 200);
    assert!(accepted(&unknown, &bars(), 0.05, |_| true).is_empty());
}

#[test]
fn gathering_scores_are_read_per_model() {
    let hmm = std::env::temp_dir().join("rosella_fragments_test.hmm");
    std::fs::write(
        &hmm,
        "NAME  First\nLENG  10\nGA    22.2 22.2;\nNAME  Second\nLENG  20\nGA    800.0 640.0;\n",
    )
    .unwrap();
    let bars = gathering(&hmm).unwrap();

    assert_eq!(bars["First"].sequence, 22.2);
    assert_eq!(bars["First"].domain, 22.2);
    assert_eq!(bars["Second"].sequence, 800.0);
    assert_eq!(bars["Second"].domain, 640.0);
}

#[test]
fn the_rescue_only_looks_at_the_proteins_it_is_given() {
    let hit = row("7", "Ribosomal_S9", 200, 51.0, 1, 100);
    assert!(accepted(&hit, &bars(), 0.05, |protein| protein != 7).is_empty());
}

/// The floor has to sit under the lowest score the span rule can accept, or hmmsearch never
/// reports the hits that rule would take.
#[test]
fn the_floor_falls_with_the_span_and_the_smallest_cutoff() {
    assert_eq!(floor(&bars(), 0.3), "30.00");
    assert_eq!(floor(&Bars::new(), 0.3), "0");
}

#[test]
fn a_whole_protein_needs_both_cutoffs() {
    let split = Bars::from([(
        "secA".to_string(),
        Bar {
            sequence: 800.0,
            domain: 640.0,
        },
    )]);
    let sequence_at = |score: &str| {
        let mut fields = vec!["-"; 23];
        fields[0] = "7";
        fields[3] = "secA";
        fields[5] = "200";
        fields[7] = score;
        fields[13] = "700";
        fields[15] = "1";
        fields[16] = "200";
        fields.join(" ")
    };
    assert!(complete(&sequence_at("810"), &split).contains_key(&7));
    assert!(complete(&sequence_at("790"), &split).is_empty());
}

// A gathering cutoff is the lowest true positive in the seed and a noise cutoff is the highest
// known false positive, so the band between them is un-excluded rather than rejected.
#[test]
fn a_model_with_a_wide_noise_gap_is_searched_below_its_gathering_score() {
    let hmm = std::env::temp_dir().join("rosella_noise_gap_test.hmm");
    std::fs::write(
        &hmm,
        [
            "NAME  Wide",
            "GA    223.55 223.55;",
            "NC    101.7 101.7;",
            "//",
            "NAME  Tight",
            "GA    100.0 80.0;",
            "NC    100.0 80.0;",
            "//",
            "NAME  Noiseless",
            "GA    50.0 40.0;",
            "//",
        ]
        .join("\n"),
    )
    .unwrap();

    let bars = gathering(&hmm).unwrap();

    assert!((bars["Wide"].sequence - (223.55f64 * 101.7).sqrt()).abs() < 1e-9);
    assert_eq!(bars["Tight"].sequence, 100.0);
    assert_eq!(bars["Tight"].domain, 80.0);
    assert_eq!(bars["Noiseless"].sequence, 50.0);
    assert_eq!(bars["Noiseless"].domain, 40.0);
}
