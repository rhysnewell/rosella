//! Layouts captured from CoverM 0.8.0. Two bugs live here, and both are invisible at one
//! sample. The long read table repeats the contig length once per sample, which the record
//! parser read as coverage. And both parsers grouped the row as every mean then every
//! variance, while `metabat` reads a row as interleaved pairs.

use std::io::Write;

use rosella::coverage::coverage_table::CoverageTable;
use rosella::external::coverm_engine::MappingMode;

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
    let table = CoverageTable::from_file(file.path(), MappingMode::LongBam).unwrap();

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
fn metabat_columns_alternate_coverage_and_variance() {
    let file = write(SHORT_TWO_SAMPLES);
    let table = CoverageTable::from_file(file.path(), MappingMode::ShortBam).unwrap();

    assert_eq!(table.sample_names, ["s1.fq", "s2.fq"]);
    assert_eq!(table.table.dim(), (2, 4));
    assert_eq!(
        table.table.row(1).to_vec(),
        [10.3393, 9.6535, 9.9681, 11.2406]
    );
    assert_eq!(table.average_depths, [10.8163, 10.1537]);
}

/// A `--coverage-file` rosella did not write can be either layout, and passing the long
/// read one through the metabat parser used to read a length column as coverage.
#[test]
fn layout_comes_from_the_header_when_the_mode_is_unknown() {
    let long = write(LONG_TWO_SAMPLES);
    let short = write(SHORT_TWO_SAMPLES);

    let from_long = CoverageTable::from_any_file(long.path()).unwrap();
    let from_short = CoverageTable::from_any_file(short.path()).unwrap();

    assert_eq!(from_long.sample_names, ["s1.fq", "s2.fq"]);
    assert_eq!(from_long.table.row(0)[0], 10.380756);
    assert_eq!(from_short.sample_names, ["s1.fq", "s2.fq"]);
    assert_eq!(from_short.table.row(0)[0], 10.5133);
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
    let table = CoverageTable::from_file(file.path(), MappingMode::ShortBam).unwrap();
    assert_eq!(table.sample_names, ["run.bam2"]);
}

/// `metabat` reads a row as interleaved per-sample mean and variance. The parser used to
/// group every mean ahead of every variance, so from two samples up the distance was
/// computed with one sample's variance standing in for the next sample's mean.
#[test]
fn a_parsed_row_is_interleaved_mean_and_variance() {
    let file = write(SHORT_TWO_SAMPLES);
    let table = CoverageTable::from_file(file.path(), MappingMode::ShortBam).unwrap();

    let row = table.table.row(0).to_vec();
    assert_eq!(row, [10.5133, 12.2708, 11.1193, 13.4039]);
    assert_eq!(
        row.iter().step_by(2).copied().collect::<Vec<_>>(),
        [10.5133, 11.1193]
    );
}

/// The header rosella writes names the columns in interleaved pairs, so the row it writes
/// has to be interleaved too or a table cannot survive its own round trip.
#[test]
fn a_written_table_reads_back_unchanged() {
    let file = write(SHORT_TWO_SAMPLES);
    let mut table = CoverageTable::from_file(file.path(), MappingMode::ShortBam).unwrap();
    let out = tempfile::NamedTempFile::new().unwrap();
    table.write(out.path()).unwrap();

    let back = CoverageTable::from_any_file(out.path()).unwrap();
    assert_eq!(back.sample_names, table.sample_names);
    // `write` rounds to three decimals, so only the ordering survives exactly.
    for (written, read) in table.table.iter().zip(back.table.iter()) {
        assert!((written - read).abs() < 1e-3, "{written} became {read}");
    }
}
