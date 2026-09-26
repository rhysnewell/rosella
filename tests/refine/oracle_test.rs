//! Reading a contig to genome map out of a CAMI binning file.

use std::io::Write;

use rosella::refine::oracle::read_groups;

#[test]
fn the_header_lines_are_skipped_and_the_contigs_group_by_bin() {
    let mut file = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        file,
        "@Version:0.9.1\n@SampleID:probe\n\n@@SEQUENCEID\tBINID\tTAXID\t_LENGTH\n\
         C0\tOTU_1\t1485\t498\nC2\tOTU_1\t1485\t550\nC1\tOTU_2\t1485\t251\n\
         C9\tOTU_3\t1485\t300"
    )
    .unwrap();

    let names = ["C0", "C1", "C2"].map(String::from);
    let groups = read_groups(file.path().to_str().unwrap(), &names).unwrap();

    assert_eq!(groups, vec![vec![0, 2], vec![1]]);
}
