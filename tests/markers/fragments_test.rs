use std::collections::HashMap;

use rosella::markers::fragments::{accepted, gathering};

fn row(protein: &str, model: &str, length: u32, score: f64, from: u32, to: u32) -> String {
    let mut fields = vec!["-"; 23];
    let score = score.to_string();
    let length = length.to_string();
    let from = from.to_string();
    let to = to.to_string();
    fields[0] = protein;
    fields[3] = model;
    fields[5] = &length;
    fields[13] = &score;
    fields[15] = &from;
    fields[16] = &to;
    fields.join(" ")
}

fn bars() -> HashMap<String, f64> {
    HashMap::from([("Ribosomal_S9".to_string(), 100.0)])
}

#[test]
fn a_half_model_hit_needs_half_the_gathering_score() {
    let over = row("7", "Ribosomal_S9", 200, 51.0, 1, 100);
    let under = row("8", "Ribosomal_S9", 200, 49.0, 1, 100);
    let taken = accepted(&[over, under].join("\n"), &bars(), 0.05);

    assert_eq!(taken.len(), 1);
    assert!(taken.contains_key("7"));
}

#[test]
fn a_sliver_of_a_model_is_refused_however_it_scores() {
    let sliver = row("9", "Ribosomal_S9", 200, 900.0, 1, 4);
    assert!(accepted(&sliver, &bars(), 0.05).is_empty());
}

#[test]
fn a_model_the_table_does_not_carry_is_refused() {
    let unknown = row("9", "Nothing", 200, 900.0, 1, 200);
    assert!(accepted(&unknown, &bars(), 0.05).is_empty());
}

#[test]
fn gathering_scores_are_read_per_model() {
    let hmm = std::env::temp_dir().join("rosella_fragments_test.hmm");
    std::fs::write(
        &hmm,
        "NAME  First\nLENG  10\nGA    22.2 22.2;\nNAME  Second\nLENG  20\nGA    140.0 140.0;\n",
    )
    .unwrap();
    let bars = gathering(&hmm).unwrap();

    assert_eq!(bars.get("First"), Some(&22.2));
    assert_eq!(bars.get("Second"), Some(&140.0));
}
