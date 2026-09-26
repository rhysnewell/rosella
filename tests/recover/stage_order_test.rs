use rosella::recover::recover_engine::{SHIPPED_ORDER, Stage, parse_order, stage_label};

#[test]
fn an_order_keeps_repeats_and_refuses_a_name_that_is_not_a_stage() {
    assert_eq!(
        parse_order(SHIPPED_ORDER).unwrap(),
        vec![
            Stage::Dissolve,
            Stage::Join,
            Stage::Recruit,
            Stage::Audit,
            Stage::Shed
        ]
    );
    assert_eq!(
        parse_order("join, dissolve ,join").unwrap(),
        vec![Stage::Join, Stage::Dissolve, Stage::Join]
    );
    assert!(parse_order("dissolve,peel").is_err());
}

/// Dropping a stage is a legitimate arm rather than an error, which is the only reason the
/// omission warning has to carry the whole weight of catching a typo.
#[test]
fn an_order_may_leave_a_stage_out() {
    assert_eq!(
        parse_order("dissolve,join,recruit,audit").unwrap(),
        vec![Stage::Dissolve, Stage::Join, Stage::Recruit, Stage::Audit]
    );
}

#[test]
fn every_pass_of_a_stage_gets_a_census_row_of_its_own() {
    let labels = (0..3)
        .map(|pass| stage_label("shed", pass))
        .collect::<Vec<_>>();
    assert_eq!(labels, vec!["shed", "shed_2", "shed_3"]);
}
