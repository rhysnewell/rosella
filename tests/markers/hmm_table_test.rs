use rosella::markers::hmm_table::{
    Columns, DOMAIN_HMM_FROM, DOMAIN_HMM_TO, DOMAIN_MODEL, DOMAIN_MODEL_LENGTH, DOMAIN_SCORE,
    DOMAIN_SEQUENCE_SCORE, DOMAIN_TARGET, Hits, keep_best,
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
    keep_best(&mut forwards, 1, "bravo", 5.0);
    keep_best(&mut forwards, 1, "alpha", 5.0);

    let mut backwards = Hits::new();
    keep_best(&mut backwards, 1, "alpha", 5.0);
    keep_best(&mut backwards, 1, "bravo", 5.0);

    assert_eq!(forwards[&1].0, "alpha");
    assert_eq!(backwards[&1].0, "alpha");

    keep_best(&mut forwards, 1, "zulu", 5.5);
    assert_eq!(forwards[&1].0, "zulu");
}
