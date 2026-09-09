//! The zero replacement in the centre log ratio transform, which has to follow contig length.

use ndarray::Array2;
use rosella::kmers::kmer_counting::KmerFrequencyTable;

const SHORT: usize = 2_000;
const LONG: usize = 20_000;
const ZERO_COLUMN: usize = 2;

fn table(n_rows: usize) -> KmerFrequencyTable {
    let rows =
        Array2::from_shape_vec((n_rows, 4), [0.5, 0.5, 0.0, 0.0].repeat(n_rows)).unwrap();
    KmerFrequencyTable::new(
        4,
        rows,
        (0..n_rows).map(|i| format!("contig_{i}")).collect(),
        String::new(),
    )
}

#[test]
fn replacement_follows_contig_length() {
    let mut table = table(2);
    table.clr(&[SHORT, LONG]).unwrap();

    let short = table.kmer_table[[0, ZERO_COLUMN]];
    let long = table.kmer_table[[1, ZERO_COLUMN]];

    assert!(
        short > long,
        "the shorter contig cannot express as small a frequency, so its replacement sits \
         higher: {short} against {long}"
    );

    let positions = |length: usize| (length - 3) as f64;
    let expected = 0.5 * (positions(LONG) / positions(SHORT)).ln();
    assert!(
        (short - long - expected).abs() < 1e-2,
        "gap {} against the {expected} the length ratio predicts",
        short - long
    );
}

#[test]
fn a_length_per_row_is_required() {
    assert!(table(2).clr(&[SHORT]).is_err());
}
