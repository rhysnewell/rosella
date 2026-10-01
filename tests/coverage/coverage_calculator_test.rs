use std::io::Write;

use rosella::cli::CoverageSource;
use rosella::coverage::coverage_calculator::ReadCollection;
use rosella::coverage::coverage_table::CoverageTable;

/// A rerun skips the samples the table already holds, so a BAM has to be named the way CoverM
/// heads its column or every rerun appends that sample again.
#[test]
fn a_bam_is_named_as_the_column_coverm_writes_for_it() {
    let source = CoverageSource {
        coverage_file: None,
        read1: Vec::new(),
        read2: Vec::new(),
        coupled: Vec::new(),
        interleaved: Vec::new(),
        single: Vec::new(),
        longreads: Vec::new(),
        bam_files: vec!["mapped/s1.bam".to_string()],
        longread_bam_files: Vec::new(),
    };
    let mut table = tempfile::NamedTempFile::new().unwrap();
    write!(
        table,
        "contigName\tcontigLen\ttotalAvgDepth\ts1\ts1-var\nc1\t1500\t1.0\t1.0\t0.5\n"
    )
    .unwrap();
    table.flush().unwrap();

    let reads = ReadCollection::new(&source).unwrap();
    assert_eq!(
        reads.sample_names(),
        CoverageTable::sample_names_in(table.path()).unwrap()
    );
}
