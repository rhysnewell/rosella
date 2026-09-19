use ndarray::Array2;
use rosella::kmers::pca::project;

fn planted(rows: usize, cols: usize, signal: usize) -> Array2<f64> {
    Array2::from_shape_fn((rows, cols), |(row, col)| {
        let phase = (row % 7) as f64;
        if col < signal {
            phase * ((col + 1) as f64)
        } else {
            ((row * 31 + col * 17) % 5) as f64 * 1e-6
        }
    })
}

#[test]
fn a_table_whose_signal_is_a_few_directions_keeps_the_floor_not_the_whole_width() {
    let mut table = planted(400, 178, 6);
    project(&mut table, 0.75).unwrap();

    assert_eq!(table.nrows(), 400);
    assert_eq!(table.ncols(), 40);
}

#[test]
fn a_table_narrower_than_the_floor_is_not_widened() {
    let mut table = planted(80, 20, 4);
    project(&mut table, 0.75).unwrap();

    assert_eq!(table.ncols(), 20);
}

#[test]
fn a_target_the_leading_components_cannot_reach_walks_up_to_the_cap() {
    let mut table = Array2::from_shape_fn((300, 178), |(row, col)| {
        (((row * 131 + col * 79) % 997) as f64 / 997.0) - 0.5
    });
    project(&mut table, 0.99).unwrap();

    assert_eq!(table.ncols(), 60);
}
