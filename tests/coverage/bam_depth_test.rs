use std::fs::File;
use std::num::NonZero;
use std::path::Path;
use std::time::{Duration, SystemTime};

use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        io::Write as _,
        record::{
            Flags,
            cigar::{Op, op::Kind},
        },
        record_buf::RecordBuf,
    },
    header::record::value::{Map, map::ReferenceSequence},
};

use rosella::coverage::bam_depth::{Sample, read};
use rosella::learned::windows::Windows;

const LENGTHS: [usize; 3] = [250, 120, 2000];

struct Hit {
    contig: usize,
    start: usize,
    ops: Vec<(Kind, usize)>,
    flags: Flags,
}

fn hit(contig: usize, start: usize, ops: &[(Kind, usize)]) -> Hit {
    Hit {
        contig,
        start,
        ops: ops.to_vec(),
        flags: Flags::empty(),
    }
}

fn header() -> sam::Header {
    LENGTHS
        .iter()
        .enumerate()
        .fold(sam::Header::builder(), |builder, (at, length)| {
            builder.add_reference_sequence(
                format!("c{at}"),
                Map::<ReferenceSequence>::new(NonZero::new(*length).unwrap()),
            )
        })
        .build()
}

fn write_bam(path: &Path, hits: &[Hit]) {
    let header = header();
    let mut writer = bam::io::Writer::new(File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    for hit in hits {
        let record = RecordBuf::builder()
            .set_flags(hit.flags)
            .set_reference_sequence_id(hit.contig)
            .set_alignment_start(Position::try_from(hit.start + 1).unwrap())
            .set_cigar(
                hit.ops
                    .iter()
                    .map(|(kind, length)| Op::new(*kind, *length))
                    .collect(),
            )
            .build();
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn hand_depth(hits: &[Hit]) -> Vec<Vec<f64>> {
    let mut depth = LENGTHS.map(|length| vec![0.0; length]).to_vec();
    for hit in hits.iter().filter(|hit| hit.flags.is_empty()) {
        let span = hit
            .ops
            .iter()
            .filter(|(kind, _)| matches!(kind, Kind::Match | Kind::Deletion | Kind::Skip))
            .map(|(_, length)| length)
            .sum::<usize>();
        let end = (hit.start + span).min(LENGTHS[hit.contig]);
        for base in &mut depth[hit.contig][hit.start..end] {
            *base += 1.0;
        }
    }
    depth
}

fn f16_value(bits: u16) -> f64 {
    let exponent = i32::from((bits >> 10) & 0x1f);
    let mantissa = f64::from(bits & 0x3ff);
    match exponent {
        0 => mantissa * 2f64.powi(-24),
        _ => (1.0 + mantissa / 1024.0) * 2f64.powi(exponent - 15),
    }
}

fn hits() -> Vec<Hit> {
    let mut hits = vec![
        hit(
            0,
            10,
            &[
                (Kind::SoftClip, 5),
                (Kind::Match, 30),
                (Kind::Deletion, 4),
                (Kind::Match, 20),
            ],
        ),
        hit(
            0,
            30,
            &[(Kind::Match, 40), (Kind::Insertion, 3), (Kind::Match, 37)],
        ),
        hit(0, 220, &[(Kind::Match, 100)]),
        hit(2, 0, &[(Kind::Match, 150)]),
        hit(
            2,
            90,
            &[(Kind::Match, 60), (Kind::Skip, 500), (Kind::Match, 40)],
        ),
        hit(2, 1950, &[(Kind::Match, 50)]),
    ];
    for flags in [Flags::SECONDARY, Flags::SUPPLEMENTARY, Flags::UNMAPPED] {
        hits.push(Hit {
            flags,
            ..hit(2, 1000, &[(Kind::Match, 150)])
        });
    }
    hits.sort_by_key(|hit| (hit.contig, hit.start));
    hits
}

#[test]
fn a_bam_reads_as_plain_depth_with_secondary_supplementary_and_unmapped_left_out() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("s1.bam");
    let hits = hits();
    write_bam(&bam, &hits);
    let depth = directory.path().join("depth");
    let samples = [Sample {
        name: "s1".to_string(),
        path: bam,
    }];

    let table = read(&samples, &depth, 2).unwrap();

    let windows = std::fs::read(depth.join("s1.win.f16")).unwrap();
    let mut at = 0;
    for (contig, bases) in hand_depth(&hits).iter().enumerate() {
        let mean = bases.iter().sum::<f64>() / bases.len() as f64;
        let variance =
            bases.iter().map(|base| base * base).sum::<f64>() / bases.len() as f64 - mean * mean;
        let row = table.table.row(contig);
        assert!((row[0] - mean).abs() < 1e-12, "c{contig} mean {}", row[0]);
        assert!(
            (row[1] - variance).abs() < 1e-9,
            "c{contig} variance {}",
            row[1]
        );
        for window in bases.chunks(100) {
            let expected = window.iter().sum::<f64>() / window.len() as f64;
            let given = f16_value(u16::from_le_bytes([windows[at], windows[at + 1]]));
            assert!(
                (given - expected).abs() <= expected / 1024.0,
                "c{contig} window {} reads {given} against {expected}",
                at / 2
            );
            at += 2;
        }
    }
    assert_eq!(at, windows.len());
    assert!(
        Windows::read(&depth, &["s1".to_string()]).is_ok(),
        "the learned views read what the BAM pass wrote"
    );
}

#[test]
fn kept_depth_is_reused_until_its_bam_is_newer() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("s1.bam");
    let depth = directory.path().join("depth");
    let samples = [Sample {
        name: "s1".to_string(),
        path: bam.clone(),
    }];
    write_bam(&bam, &hits());
    let first = read(&samples, &depth, 1).unwrap();

    write_bam(&bam, &[hit(0, 0, &[(Kind::Match, 250)])]);
    let past = SystemTime::now() - Duration::from_secs(3600);
    File::options()
        .write(true)
        .open(&bam)
        .unwrap()
        .set_modified(past)
        .unwrap();
    assert_eq!(read(&samples, &depth, 1).unwrap().table, first.table);

    File::options()
        .write(true)
        .open(&bam)
        .unwrap()
        .set_modified(SystemTime::now() + Duration::from_secs(3600))
        .unwrap();
    assert_eq!(read(&samples, &depth, 1).unwrap().table.row(0)[0], 1.0);
}

#[test]
fn a_bam_out_of_coordinate_order_is_refused() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("s1.bam");
    write_bam(
        &bam,
        &[
            hit(2, 0, &[(Kind::Match, 50)]),
            hit(0, 0, &[(Kind::Match, 50)]),
        ],
    );
    let samples = [Sample {
        name: "s1".to_string(),
        path: bam,
    }];
    let Err(error) = read(&samples, &directory.path().join("depth"), 1) else {
        panic!("a BAM out of coordinate order was read");
    };
    assert!(format!("{error:#}").contains("not sorted"), "{error:#}");
}

#[test]
fn a_depth_directory_made_one_sample_at_a_time_stands_in_for_its_bams() {
    use clap::Parser;
    use rosella::cli::{Cli, Command};
    use rosella::coverage::coverage_calculator::{CoverageInputs, calculate_coverage};
    use rosella::coverage::coverage_table::CoverageTable;
    use rosella::coverage::make::run_make;

    let directory = tempfile::tempdir().unwrap();
    let made = directory.path().join("made");
    let parse = |args: &[&str]| Cli::try_parse_from(args).unwrap().command;
    for (name, hits) in [
        ("s2", vec![hit(0, 0, &[(Kind::Match, 250)])]),
        ("s1", hits()),
    ] {
        let bam = directory.path().join(format!("{name}.bam"));
        write_bam(&bam, &hits);
        let Command::Make(make) = parse(&[
            "rosella",
            "make",
            "-b",
            bam.to_str().unwrap(),
            "-o",
            made.to_str().unwrap(),
        ]) else {
            panic!("make parsed as another command");
        };
        run_make(&make).unwrap();
        std::fs::remove_file(&bam).unwrap();
    }

    let out = directory.path().join("out");
    let Command::Recover(recover) = parse(&[
        "rosella",
        "recover",
        "-r",
        "unread.fasta",
        "--depth-directory",
        made.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]) else {
        panic!("recover parsed as another command");
    };
    let coverage = calculate_coverage(&CoverageInputs {
        assembly: &recover.assembly,
        output_directory: &recover.common.output_directory,
        depth: &made,
        threads: 1,
        reads: &recover.reads,
    })
    .unwrap();
    assert_eq!(coverage.windows.as_deref(), Some(made.as_path()));
    assert_eq!(coverage.table.sample_names, ["s1", "s2"]);
    let first = CoverageTable::from_any_file(made.join("s1.depth.tsv")).unwrap();
    assert_eq!(coverage.table.table.column(0), first.table.column(0));
    assert_eq!(coverage.table.table.row(0)[2], 1.0);
}
