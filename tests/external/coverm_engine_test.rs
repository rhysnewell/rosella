use std::io::Write;

use flate2::Compression;
use flate2::write::GzEncoder;
use rosella::external::coverm_engine::intact;

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
