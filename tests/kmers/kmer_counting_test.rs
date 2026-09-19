use ndarray::Array2;
use rosella::kmers::kmer_counting::{KmerFrequencyTable, canonical_count};

const SHORT: usize = 2_000;
const LONG: usize = 20_000;
const ZERO_COLUMN: usize = 2;

fn table(n_rows: usize) -> KmerFrequencyTable {
    let mut row = vec![0.0; canonical_count(2)];
    row[0] = 0.5;
    row[1] = 0.5;
    let rows = Array2::from_shape_vec((n_rows, row.len()), row.repeat(n_rows)).unwrap();
    KmerFrequencyTable::new(
        2,
        rows,
        (0..n_rows).map(|i| format!("contig_{i}")).collect(),
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

    let positions = |length: usize| (length - 1) as f64;
    let expected = (positions(LONG) / positions(SHORT)).ln();
    let zeros = canonical_count(2) - 2;
    let expected = expected * (1.0 - zeros as f64 / canonical_count(2) as f64);
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

#[test]
fn a_row_is_centred_on_its_own_geometric_mean() {
    let mut table = table(1);
    table.clr(&[LONG]).unwrap();

    let row = table.kmer_table.row(0).sum();
    assert!(
        row.abs() < 1e-9,
        "a centre log ratio row sums to zero, this one sums to {row}"
    );
}

#[test]
fn a_written_table_reports_the_k_it_was_built_from() {
    let width = canonical_count(3);
    let mut written = KmerFrequencyTable::new(
        3,
        Array2::from_shape_vec((1, width), vec![1.0 / width as f64; width]).unwrap(),
        vec!["contig_0".to_string()],
    );
    let file = tempfile::NamedTempFile::new().unwrap();
    written.write(file.path()).unwrap();

    let read = KmerFrequencyTable::read(file.path()).unwrap();
    assert_eq!(read.kmer_size(), 3);
}

#[test]
fn a_width_that_no_k_reaches_is_refused() {
    let odd = canonical_count(2) + 1;
    let mut written = KmerFrequencyTable::new(
        2,
        Array2::from_shape_vec((1, odd), vec![1.0 / odd as f64; odd]).unwrap(),
        vec!["contig_0".to_string()],
    );
    let file = tempfile::NamedTempFile::new().unwrap();
    written.write(file.path()).unwrap();

    assert!(KmerFrequencyTable::read(file.path()).is_err());
}
