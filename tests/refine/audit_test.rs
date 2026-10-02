use rosella::refine::audit::share;

#[test]
fn share_does_not_move_with_the_order_the_neighbours_arrived() {
    let mut arrived = [(2usize, 1e16f64), (0, 1.0), (1, 1.0)];
    let mut reversed = [(1usize, 1.0f64), (0, 1.0), (2, 1e16)];

    let naive = arrived.iter().map(|(_, weight)| *weight).sum::<f64>();
    let own = share(&mut arrived, 0).unwrap();

    assert_eq!(own, share(&mut reversed, 0).unwrap());
    assert_ne!(own, 1.0 / naive);
}

#[test]
fn nothing_binned_nearby_is_not_a_zero_share() {
    assert_eq!(share(&mut [], 0), None);
    assert_eq!(share(&mut [(0, 0.0), (1, 0.0)], 0), None);
}
