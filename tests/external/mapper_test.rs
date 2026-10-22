use std::fs::File;
use std::io::Write;

use flate2::Compression;
use flate2::write::GzEncoder;
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write as _};
use rosella::coverage::bam_depth::complete;
use rosella::external::mapper::{index_parts, intact, write_bam};

fn member(reads: usize) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    for read in 0..reads {
        writeln!(
            encoder,
            "@r{read}\n{}\n+\n{}",
            "ACGT".repeat(40 + read % 7),
            "I".repeat(160 + 4 * (read % 7))
        )
        .unwrap();
    }
    encoder.finish().unwrap()
}

#[test]
fn a_read_file_damaged_past_its_first_gzip_member_is_refused() {
    let directory = tempfile::tempdir().unwrap();
    let (first, second) = (member(500), member(500));
    let whole = [first.as_slice(), second.as_slice()].concat();
    let mut damaged = whole.clone();
    let middle = first.len() + second.len() / 2;
    damaged[middle..middle + 16]
        .iter_mut()
        .for_each(|byte| *byte ^= 0xff);
    let plain = b"@r0\nACGT\n+\nIIII\n".to_vec();
    for (name, bytes, ok) in [
        ("whole.fq.gz", whole, true),
        ("damaged.fq.gz", damaged, false),
        ("plain.fq", plain, true),
    ] {
        let path = directory.path().join(name);
        std::fs::write(&path, bytes).unwrap();
        assert_eq!(intact(&path).is_ok(), ok, "{name}");
    }
}

#[test]
fn only_finished_indexes_of_this_assembly_and_mapper_are_reused() {
    let prefix = "assembly.fasta.gz.minimap2-sr.";
    for (name, parts) in [
        ("assembly.fasta.gz.minimap2-sr.mmi", Some(1)),
        ("assembly.fasta.gz.minimap2-sr.7.mmi", Some(7)),
        ("assembly.fasta.gz.minimap2-sr.7.mmi.partial", None),
        ("._assembly.fasta.gz.minimap2-sr.7.mmi", None),
        ("assembly.fasta.gz.minimap2-ont.7.mmi", None),
        ("assembly.fasta.gz.minimap2-sr.x.mmi", None),
    ] {
        assert_eq!(index_parts(name, prefix), parts, "{name}");
    }
}

#[test]
fn mapped_sam_is_kept_as_a_whole_bam_with_every_hit_in_mapper_order() {
    let records = [
        "r1\t0\tc2\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\tAS:i:10",
        "r1\t256\tc1\t40\t0\t4M2I4M\t*\t0\t0\t*\t*\tAS:i:10",
        "r2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII",
    ];
    let body = format!("{}\n", records.join("\n"));
    let sam = format!("@SQ\tSN:c1\tLN:100\n@SQ\tSN:c2\tLN:50\n{body}");
    let home = tempfile::tempdir().unwrap();
    let path = home.path().join("sample.bam");
    write_bam(sam.as_bytes(), &path, 3).unwrap();

    assert!(complete(&path).unwrap());
    let mut reader = bam::io::Reader::new(File::open(&path).unwrap());
    let header = reader.read_header().unwrap();
    let names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|name| name.to_string())
        .collect();
    assert_eq!(names, ["c1", "c2"]);
    let mut text = sam::io::Writer::new(Vec::new());
    for record in reader.records() {
        text.write_alignment_record(&header, &record.unwrap())
            .unwrap();
    }
    assert_eq!(String::from_utf8(text.into_inner()).unwrap(), body);
}
