use rosella::recover::recover_engine::{SHIPPED_ORDER, Stage, parse_order};

#[test]
fn an_order_keeps_repeats_and_refuses_a_name_that_is_not_a_stage() {
    assert_eq!(
        parse_order(SHIPPED_ORDER).unwrap(),
        vec![Stage::Dissolve, Stage::Join, Stage::Recruit]
    );
    assert_eq!(
        parse_order("join, dissolve ,join").unwrap(),
        vec![Stage::Join, Stage::Dissolve, Stage::Join]
    );
    assert!(parse_order("dissolve,shed").is_err());
}
