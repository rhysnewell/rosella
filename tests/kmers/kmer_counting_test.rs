use std::io::Write;

use ndarray::Array2;
use rosella::kmers::kmer_counting::{
    KMER_SIZES, KmerFrequencyTable, canonical_count, canonical_index, prefixes,
};
use rosella::kmers::scan::{Floors, scan};

const SHORT: usize = 2_000;
const LONG: usize = 20_000;
const ZERO_COLUMN: usize = 2;

fn counted(path: &str, floor: usize, kmer_size: usize) -> KmerFrequencyTable {
    let floors = Floors {
        composition: Some(floor),
        ..Floors::default()
    };
    scan(path, kmer_size, floors).unwrap().composition
}

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
fn a_written_table_reports_the_k_it_was_built_from() {
    let width = canonical_count(3);
    let written = KmerFrequencyTable::new(
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
fn the_table_is_written_with_tabs_and_an_older_comma_table_still_reads() {
    let width = canonical_count(2);
    let rows = (0..width).map(|at| at as f64 / 7.0).collect::<Vec<_>>();
    let written = KmerFrequencyTable::new(
        2,
        Array2::from_shape_vec((1, width), rows).unwrap(),
        vec!["contig_0".to_string()],
    );
    let file = tempfile::NamedTempFile::new().unwrap();
    written.write(file.path()).unwrap();
    let tabbed = std::fs::read_to_string(file.path()).unwrap();
    assert_eq!(tabbed.trim_end().split('\t').count(), width + 1);

    let comma = tempfile::NamedTempFile::new().unwrap();
    std::fs::write(comma.path(), tabbed.replace('\t', ",")).unwrap();
    let read = KmerFrequencyTable::read(comma.path()).unwrap();
    assert_eq!(read.kmer_table, written.kmer_table);
}

#[test]
fn an_empty_table_is_refused_rather_than_read() {
    let file = tempfile::NamedTempFile::new().unwrap();
    assert!(KmerFrequencyTable::read(file.path()).is_err());
}

#[test]
fn a_width_that_no_k_reaches_is_refused() {
    let odd = canonical_count(2) + 1;
    let mut table = KmerFrequencyTable::new(
        2,
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

    let mut whole = counted(path, 0, 4);
    whole.clr(&lengths).unwrap();
    let band = prefixes(path, &[("c2", lengths[2]), ("c1", lengths[1])], 4).unwrap();

    let mut long = counted(path, 1_000, 4);
    long.clr(&[lengths[0], lengths[2]]).unwrap();

    assert_eq!(band.row(0), whole.kmer_table.row(2));
    assert_eq!(band.row(1), whole.kmer_table.row(1));
    assert_eq!(long.kmer_table.nrows(), 2);
    assert_eq!(long.kmer_table.row(0), whole.kmer_table.row(0));
    assert_eq!(long.kmer_table.row(1), whole.kmer_table.row(2));
}

#[test]
fn a_contig_with_no_kmer_reads_as_zeros_and_stays_finite() {
    let mut assembly = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        assembly,
        ">gap\n{}\n>real\n{}",
        "N".repeat(500),
        sequence(500, 1)
    )
    .unwrap();
    assembly.flush().unwrap();

    let path = assembly.path().to_str().unwrap();
    let mut table = counted(path, 0, 4);
    assert!(table.kmer_table.row(0).iter().all(|value| *value == 0.0));

    table.clr(&[500, 500]).unwrap();
    assert!(table.kmer_table.iter().all(|value| value.is_finite()));
}

fn noisy(length: usize, mut state: u64) -> Vec<u8> {
    const ALPHABET: &[u8] = b"ACGTACGTacgtacgtNnRy-";
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ALPHABET[(state >> 33) as usize % ALPHABET.len()]
        })
        .collect()
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'A',
    }
}

fn reference(sequence: &[u8], kmer_size: usize) -> Vec<f64> {
    let index = canonical_index(kmer_size);
    let mut counts = vec![0u32; index.len()];
    let mut total = 0u32;
    for kmer in sequence.to_ascii_uppercase().windows(kmer_size) {
        if !kmer.iter().all(|base| b"ACGT".contains(base)) {
            continue;
        }
        let reverse = kmer
            .iter()
            .rev()
            .map(|base| complement(*base))
            .collect::<Vec<_>>();
        counts[index[kmer.min(&reverse[..])]] += 1;
        total += 1;
    }
    let total = f64::from(total.max(1));
    counts
        .iter()
        .map(|count| f64::from(*count) / total)
        .collect()
}

#[test]
fn every_k_counts_what_a_window_by_window_reference_counts() {
    let lengths = [0, 1, 5, 6, 7, 64, 1_000, 5_003];
    let mut assembly = tempfile::NamedTempFile::new().unwrap();
    let contigs = lengths
        .iter()
        .enumerate()
        .map(|(at, length)| noisy(*length, at as u64 + 7))
        .collect::<Vec<_>>();
    for (at, contig) in contigs.iter().enumerate() {
        writeln!(assembly, ">c{at}\n{}", String::from_utf8_lossy(contig)).unwrap();
    }
    assembly.flush().unwrap();
    let path = assembly.path().to_str().unwrap();

    for kmer_size in *KMER_SIZES.start() as usize..=*KMER_SIZES.end() as usize {
        let table = counted(path, 0, kmer_size);
        for (at, contig) in contigs.iter().enumerate() {
            assert_eq!(
                table.kmer_table.row(at).to_vec(),
                reference(contig, kmer_size),
                "contig {at} at k={kmer_size}"
            );
        }
    }
}

// The weight reads each half from counts taken in the one scan. A half has to give the row
// that counting it as a contig of its own gives.
#[test]
fn a_half_counted_in_the_scan_matches_the_half_counted_alone() {
    let contigs = [5_003, 6_000, 4_001]
        .iter()
        .enumerate()
        .map(|(at, length)| noisy(*length, at as u64 + 3))
        .collect::<Vec<_>>();
    let mut whole = tempfile::NamedTempFile::new().unwrap();
    let mut split = tempfile::NamedTempFile::new().unwrap();
    let mut lengths = Vec::new();
    for (at, contig) in contigs.iter().enumerate() {
        writeln!(whole, ">c{at}\n{}", String::from_utf8_lossy(contig)).unwrap();
        let (first, second) = contig.split_at(contig.len() / 2);
        for (side, piece) in [first, second].iter().enumerate() {
            writeln!(split, ">c{at}_{side}\n{}", String::from_utf8_lossy(piece)).unwrap();
            lengths.push(piece.len());
        }
    }
    whole.flush().unwrap();
    split.flush().unwrap();

    let floors = Floors {
        halves: Some(4_500),
        ..Floors::default()
    };
    let halves = scan(whole.path().to_str().unwrap(), 4, floors)
        .unwrap()
        .halves;
    let [first, second] = halves.composition(&["c1", "c0"], 4).unwrap();
    assert!(halves.composition(&["c2"], 4).is_err());

    let mut alone = counted(split.path().to_str().unwrap(), 0, 4);
    alone.clr(&lengths).unwrap();
    assert_eq!(first.row(0), alone.kmer_table.row(2));
    assert_eq!(second.row(0), alone.kmer_table.row(3));
    assert_eq!(first.row(1), alone.kmer_table.row(0));
    assert_eq!(second.row(1), alone.kmer_table.row(1));
}
