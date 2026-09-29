use std::io::Write;

use ndarray::Array2;
use rosella::kmers::kmer_counting::{
    KmerFrequencyTable, KmerSizes, canonical_count, count_kmers, prefixes,
};

const SHORT: usize = 2_000;
const LONG: usize = 20_000;
const ZERO_COLUMN: usize = 2;

fn sizes(list: &[usize]) -> KmerSizes {
    KmerSizes::from(list.to_vec())
}

fn table(n_rows: usize) -> KmerFrequencyTable {
    let mut row = vec![0.0; canonical_count(2)];
    row[0] = 0.5;
    row[1] = 0.5;
    let rows = Array2::from_shape_vec((n_rows, row.len()), row.repeat(n_rows)).unwrap();
    KmerFrequencyTable::new(
        sizes(&[2]),
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
fn each_block_is_centred_on_its_own() {
    let widths = [canonical_count(2), canonical_count(3)];
    let mut row = Vec::new();
    for width in widths {
        let mut block = vec![0.0; width];
        block[0] = 0.75;
        block[1] = 0.25;
        row.extend(block);
    }
    let mut table = KmerFrequencyTable::new(
        sizes(&[2, 3]),
        Array2::from_shape_vec((1, row.len()), row).unwrap(),
        vec!["contig_0".to_string()],
    );
    table.clr(&[LONG]).unwrap();

    let mut at = 0;
    for width in widths {
        let block = table
            .kmer_table
            .row(0)
            .slice(ndarray::s![at..at + width])
            .sum();
        assert!(
            block.abs() < 1e-9,
            "a centre log ratio block sums to zero, this one sums to {block}"
        );
        at += width;
    }
}

#[test]
fn a_written_table_reports_the_block_list_it_was_built_from() {
    let widths = [canonical_count(2), canonical_count(3), canonical_count(4)];
    let total = widths.iter().sum::<usize>();
    let mut written = KmerFrequencyTable::new(
        sizes(&[2, 3, 4]),
        Array2::from_shape_vec((1, total), vec![1.0 / total as f64; total]).unwrap(),
        vec!["contig_0".to_string()],
    );
    let file = tempfile::NamedTempFile::new().unwrap();
    written.write(file.path()).unwrap();

    let read = KmerFrequencyTable::read(file.path()).unwrap();
    assert_eq!(read.kmer_sizes(), sizes(&[2, 3, 4]));
}

#[test]
fn a_width_that_no_block_list_reaches_is_refused() {
    let odd = canonical_count(2) + 1;
    let mut table = KmerFrequencyTable::new(
        sizes(&[2]),
        Array2::from_shape_vec((1, odd), vec![1.0 / odd as f64; odd]).unwrap(),
        vec!["contig_0".to_string()],
    );
    assert!(table.clr(&[LONG]).is_err());
}

fn sequence(length: usize, mut state: u64) -> String {
    (0..length)
        .map(|at| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            match (at % 97 == 0, state >> 62) {
                (true, _) => 'N',
                (_, 0) => 'a',
                (_, 1) => 'C',
                (_, 2) => 'G',
                _ => 'T',
            }
        })
        .collect()
}

/// Input counts only the long contigs and attach counts a band when the walk reaches it, so
/// both must give the rows the whole-assembly count would have held, to the bit.
#[test]
fn long_rows_and_a_band_counted_late_match_the_whole_count() {
    let lengths = [2_400, 900, 1_300];
    let mut assembly = tempfile::NamedTempFile::new().unwrap();
    for (at, length) in lengths.iter().enumerate() {
        writeln!(assembly, ">c{at} extra\n{}", sequence(*length, at as u64)).unwrap();
    }
    assembly.flush().unwrap();
    let path = assembly.path().to_str().unwrap();
    let directory = tempfile::tempdir().unwrap();

    let mut whole = count_kmers(
        path,
        directory.path().to_str().unwrap(),
        lengths.len(),
        0,
        &sizes(&[4]),
        false,
    )
    .unwrap();
    whole.clr(&lengths).unwrap();
    let band = prefixes(
        path,
        &[("c2", lengths[2]), ("c1", lengths[1])],
        &sizes(&[4]),
    )
    .unwrap();

    let mut long = count_kmers(
        path,
        directory.path().to_str().unwrap(),
        2,
        1_000,
        &sizes(&[4]),
        false,
    )
    .unwrap();
    long.clr(&[lengths[0], lengths[2]]).unwrap();

    assert_eq!(band.row(0), whole.kmer_table.row(2));
    assert_eq!(band.row(1), whole.kmer_table.row(1));
    assert_eq!(long.kmer_table.nrows(), 2);
    assert_eq!(long.kmer_table.row(0), whole.kmer_table.row(0));
    assert_eq!(long.kmer_table.row(1), whole.kmer_table.row(2));
}
