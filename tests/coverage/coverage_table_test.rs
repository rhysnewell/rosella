//! Layouts captured from CoverM 0.8.0. Two bugs live here, and both are invisible at one
//! sample. The long read table repeats the contig length once per sample, which the record
//! parser read as coverage. And both parsers grouped the row as every mean then every
//! variance, while the coverage distance reads a row as interleaved pairs.

use std::io::Write;

use rosella::coverage::coverage_table::CoverageTable;

const LONG_TWO_SAMPLES: &str = "\
Contig\tref.fa/s1.fq Length\tref.fa/s1.fq Trimmed Mean\tref.fa/s1.fq Variance\tref.fa/s2.fq Length\tref.fa/s2.fq Trimmed Mean\tref.fa/s2.fq Variance
c1\t1500\t10.380756\t12.270764\t1500\t11.032072\t13.403855
c2\t1500\t10.449013\t9.653534\t1500\t10.011513\t11.240647
";

const SHORT_TWO_SAMPLES: &str = "\
contigName\tcontigLen\ttotalAvgDepth\tref.fa/s1.fq.bam\tref.fa/s1.fq.bam-var\tref.fa/s2.fq.bam\tref.fa/s2.fq.bam-var
c1\t1500\t10.8163\t10.5133\t12.2708\t11.1193\t13.4039
c2\t1500\t10.1537\t10.3393\t9.6535\t9.9681\t11.2406
";

fn write(contents: &str) -> tempfile::NamedTempFile {
    let mut file = tempfile::NamedTempFile::new().unwrap();
    file.write_all(contents.as_bytes()).unwrap();
    file.flush().unwrap();
    file
}

#[test]
fn long_read_triples_survive_more_than_one_sample() {
    let file = write(LONG_TWO_SAMPLES);
    let table = CoverageTable::from_any_file(file.path()).unwrap();

    assert_eq!(table.sample_names, ["s1.fq", "s2.fq"]);
    assert_eq!(table.contig_names, ["c1", "c2"]);
    assert_eq!(table.contig_lengths, [1500, 1500]);
    assert_eq!(table.table.dim(), (2, 4));
    assert_eq!(
        table.table.row(0).to_vec(),
        [10.380756, 12.270764, 11.032072, 13.403855]
    );
    assert!((table.average_depths[0] - (10.380756 + 11.032072) / 2.0).abs() < 1e-9);
}

#[test]
fn depth_table_columns_alternate_coverage_and_variance() {
    let file = write(SHORT_TWO_SAMPLES);
    let table = CoverageTable::from_any_file(file.path()).unwrap();

    assert_eq!(table.sample_names, ["s1.fq", "s2.fq"]);
    assert_eq!(table.table.dim(), (2, 4));
    assert_eq!(
        table.table.row(1).to_vec(),
        [10.3393, 9.6535, 9.9681, 11.2406]
    );
    assert_eq!(table.average_depths, [10.8163, 10.1537]);
}

#[test]
fn an_unrecognised_header_names_both_layouts() {
    let file = write("name\tvalue\nc1\t1\n");
    let error = match CoverageTable::from_any_file(file.path()) {
        Ok(_) => panic!("a two column table parsed as a coverage table"),
        Err(error) => error.to_string(),
    };
    assert!(error.contains("contigName"), "{error}");
    assert!(error.contains("Contig"), "{error}");
}

/// `.replace(".bam", "")` mangled any sample whose name contained the substring.
#[test]
fn only_a_trailing_bam_suffix_is_stripped() {
    let file = write(
        "contigName\tcontigLen\ttotalAvgDepth\trun.bam2.bam\trun.bam2.bam-var\nc1\t1500\t1.0\t1.0\t0.5\n",
    );
    let table = CoverageTable::from_any_file(file.path()).unwrap();
    assert_eq!(table.sample_names, ["run.bam2"]);
}

/// The coverage distance reads a row as interleaved per-sample mean and variance. The
/// parser used to group every mean ahead of every variance, so from two samples up the
/// distance was computed with one sample's variance standing in for the next sample's mean.
#[test]
fn a_parsed_row_is_interleaved_mean_and_variance() {
    let file = write(SHORT_TWO_SAMPLES);
    let table = CoverageTable::from_any_file(file.path()).unwrap();

    let row = table.table.row(0).to_vec();
    assert_eq!(row, [10.5133, 12.2708, 11.1193, 13.4039]);
    assert_eq!(
        row.iter().step_by(2).copied().collect::<Vec<_>>(),
        [10.5133, 11.1193]
    );
}

/// Rows under the cutoff are read back from the written table while longer rows stay in
/// memory, so the row has to come back in its header's order and to the bit. A share the tie
/// split left decaying toward zero has to write as 0, not as dozens of zeros.
#[test]
fn a_written_table_reads_back_unchanged() {
    let mut table = CoverageTable::new(
        vec!["c1".into(), "c2".into()],
        vec![1500, 820],
        vec!["s1".into()],
        ndarray::array![[2.215080346106304, 0.7769790108498187], [2.95e-32, -1e-17]],
    );
    let out = tempfile::NamedTempFile::new().unwrap();
    table.write(out.path()).unwrap();

    let text = std::fs::read_to_string(out.path()).unwrap();
    let rows = text.lines().skip(1).collect::<Vec<_>>();
    assert_eq!(
        rows,
        ["c1\t1500\t2.21508\t2.21508\t0.776979", "c2\t820\t0\t0\t0"]
    );
    let back = CoverageTable::from_any_file(out.path()).unwrap();
    assert_eq!(back.sample_names, table.sample_names);
    assert_eq!(back.table, table.table);
    assert_eq!(back.average_depths, table.average_depths);
}

#[test]
fn a_gzipped_table_reads_the_same_rows() {
    let packed = tempfile::NamedTempFile::new().unwrap();
    let mut encoder = flate2::write::GzEncoder::new(packed.as_file(), flate2::Compression::fast());
    encoder.write_all(SHORT_TWO_SAMPLES.as_bytes()).unwrap();
    encoder.finish().unwrap();
    let plain = write(SHORT_TWO_SAMPLES);

    let read = |path| CoverageTable::from_any_file(path).unwrap();
    let (zipped, flat) = (read(packed.path()), read(plain.path()));
    assert_eq!(zipped.contig_names, flat.contig_names);
    assert_eq!(zipped.table, flat.table);
}

#[test]
fn a_band_is_read_in_the_order_asked_and_a_stray_name_is_refused() {
    let file = write(SHORT_TWO_SAMPLES);
    let whole = CoverageTable::from_any_file(file.path()).unwrap();
    let band = CoverageTable::rows_named(file.path(), &["c2", "c1"]).unwrap();

    assert_eq!(band.row(0), whole.table.row(1));
    assert_eq!(band.row(1), whole.table.row(0));
    assert!(CoverageTable::rows_named(file.path(), &["c3"]).is_err());
}

#[test]
fn a_sample_already_in_the_table_is_not_merged_in_twice() {
    let read = || CoverageTable::from_any_file(write(SHORT_TWO_SAMPLES).path());
    let mut held = read().unwrap();

    assert!(held.merge(read().unwrap()).is_err());
}
