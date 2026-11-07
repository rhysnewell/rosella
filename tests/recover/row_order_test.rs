use rosella::recover::inputs::row_order;

fn laid_out(names: &[String], lengths: &[usize], defer: bool) -> Vec<String> {
    row_order(names, lengths, 1500, defer)
        .into_iter()
        .map(|row| names[row].clone())
        .collect()
}

#[test]
fn a_reordered_assembly_lays_out_the_same_rows() {
    let names = (0..300).map(|at| format!("PC{at}")).collect::<Vec<_>>();
    let lengths = (0..300)
        .map(|at| 500 + (at * 37) % 3000)
        .collect::<Vec<_>>();
    let permuted = (0..300).map(|at| (at * 7 + 11) % 300).collect::<Vec<_>>();
    let names_permuted = permuted
        .iter()
        .map(|at| names[*at].clone())
        .collect::<Vec<_>>();
    let lengths_permuted = permuted.iter().map(|at| lengths[*at]).collect::<Vec<_>>();

    for defer in [true, false] {
        let laid = laid_out(&names, &lengths, defer);
        assert_eq!(laid, laid_out(&names_permuted, &lengths_permuted, defer));
        let long = lengths.iter().filter(|length| **length >= 1500).count();
        let length_of = |name: &String| lengths[names.iter().position(|n| n == name).unwrap()];
        assert!(laid[..long].iter().all(|name| length_of(name) >= 1500));
    }
}
