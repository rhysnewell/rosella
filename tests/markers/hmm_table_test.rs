use rosella::markers::hmm_table::{
    Columns, DOMAIN_HMM_FROM, DOMAIN_HMM_TO, DOMAIN_MODEL, DOMAIN_MODEL_LENGTH, DOMAIN_SCORE,
    DOMAIN_SEQUENCE_SCORE, DOMAIN_TARGET, Hits, Reach, keep_best,
};

#[test]
fn columns_reads_the_domain_layout_in_one_pass() {
    let row = "41 - 244 PF00001.2 - 310 1e-30 99.5 0.1 1 1 2e-33 4e-30 98.7 0.0 12 297 5 251";
    let mut fields = Columns::new(row);
    assert_eq!(fields.at(DOMAIN_TARGET), Some("41"));
    assert_eq!(fields.at(DOMAIN_MODEL), Some("PF00001.2"));
    assert_eq!(fields.at(DOMAIN_MODEL_LENGTH), Some("310"));
    assert_eq!(fields.at(DOMAIN_SEQUENCE_SCORE), Some("99.5"));
    assert_eq!(fields.at(DOMAIN_SCORE), Some("98.7"));
    assert_eq!(fields.at(DOMAIN_HMM_FROM), Some("12"));
    assert_eq!(fields.at(DOMAIN_HMM_TO), Some("297"));
    assert_eq!(fields.at(64), None);
}

/// A tie has to land on the same model whichever order hmmsearch listed its rows in.
#[test]
fn a_tie_falls_to_the_name_and_a_better_score_always_wins() {
    let mut forwards = Hits::new();
    keep_best(&mut forwards, 1, "bravo", 5.0, Reach::default());
    keep_best(&mut forwards, 1, "alpha", 5.0, Reach::default());

    let mut backwards = Hits::new();
    keep_best(&mut backwards, 1, "alpha", 5.0, Reach::default());
    keep_best(&mut backwards, 1, "bravo", 5.0, Reach::default());

    assert_eq!(forwards[&1].model, "alpha");
    assert_eq!(backwards[&1].model, "alpha");

    keep_best(&mut forwards, 1, "zulu", 5.5, Reach::default());
    assert_eq!(forwards[&1].model, "zulu");
}

#[test]
fn a_second_domain_of_the_kept_model_widens_its_reach() {
    let reach = |from, to| Reach {
        model_from: from,
        model_to: to,
        model_length: 300,
        protein_from: from,
        protein_to: to,
    };
    let mut best = Hits::new();
    keep_best(&mut best, 1, "alpha", 40.0, reach(10, 120));
    keep_best(&mut best, 1, "alpha", 40.0, reach(150, 290));
    keep_best(&mut best, 1, "bravo", 30.0, reach(1, 300));

    assert_eq!(best[&1].model, "alpha");
    assert_eq!(best[&1].reach, reach(10, 290));
}
