use rosella::seeds::consistent_sample;

#[test]
fn a_pool_that_gains_one_item_keeps_its_sample() {
    let pool = (0..3000).step_by(2).collect::<Vec<_>>();
    let grown = [pool.clone(), vec![1001]].concat();
    let before = consistent_sample(&pool, 1000, 42);
    let after = consistent_sample(&grown, 1000, 42);
    let shared = after.iter().filter(|item| before.contains(item)).count();
    assert_eq!(before.len(), 1000);
    assert!(shared >= 999, "{shared} of 1000 kept");
}
