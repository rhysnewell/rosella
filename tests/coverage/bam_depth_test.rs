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
            data::field::Tag,
        },
        record_buf::{RecordBuf, data::field::Value},
    },
    header::record::value::{
        Map,
        map::{
            Header, ReferenceSequence,
            header::tag::{GROUP_ORDER, SORT_ORDER},
        },
    },
};

use rosella::coverage::bam_depth::{Sample, read};
use rosella::coverage::depth_file::{self, Contigs};
use rosella::learned::windows::Windows;

const LENGTHS: [usize; 3] = [250, 120, 2000];

struct Hit {
    contig: usize,
    start: usize,
    ops: Vec<(Kind, usize)>,
    flags: Flags,
    name: String,
    score: i32,
    edits: i32,
    others: Option<&'static str>,
}

fn hit(contig: usize, start: usize, ops: &[(Kind, usize)]) -> Hit {
    Hit {
        contig,
        start,
        ops: ops.to_vec(),
        flags: Flags::empty(),
        name: "r".to_string(),
        score: 0,
        edits: 0,
        others: None,
    }
}

fn named(name: &str, flags: Flags, contig: usize, start: usize, score: i32) -> Hit {
    Hit {
        flags,
        name: name.to_string(),
        score,
        ..hit(contig, start, &[(Kind::Match, 50)])
    }
}

fn header(in_read_order: bool) -> sam::Header {
    let builder = match in_read_order {
        true => sam::Header::builder().set_header(
            Map::<Header>::builder()
                .insert(SORT_ORDER, "unsorted")
                .insert(GROUP_ORDER, "query")
                .build()
                .unwrap(),
        ),
        false => sam::Header::builder(),
    };
    LENGTHS
        .iter()
        .enumerate()
        .fold(builder, |builder, (at, length)| {
            builder.add_reference_sequence(
                format!("c{at}"),
                Map::<ReferenceSequence>::new(NonZero::new(*length).unwrap()),
            )
        })
        .build()
}

fn write_bam(path: &Path, hits: &[Hit]) {
    write_with(path, &header(false), hits);
}

fn write_with(path: &Path, header: &sam::Header, hits: &[Hit]) {
    let mut writer = bam::io::Writer::new(File::create(path).unwrap());
    writer.write_header(header).unwrap();
    for hit in hits {
        let record = RecordBuf::builder()
            .set_name(hit.name.as_str())
            .set_flags(hit.flags)
            .set_data(
                [
                    Some((Tag::ALIGNMENT_SCORE, Value::from(hit.score))),
                    Some((Tag::EDIT_DISTANCE, Value::from(hit.edits))),
                    hit.others
                        .map(|others| (Tag::OTHER_ALIGNMENTS, Value::from(others))),
                ]
                .into_iter()
                .flatten()
                .collect(),
            )
            .set_reference_sequence_id(hit.contig)
            .set_alignment_start(Position::try_from(hit.start + 1).unwrap())
            .set_cigar(
                hit.ops
                    .iter()
                    .map(|(kind, length)| Op::new(*kind, *length))
                    .collect(),
            )
            .build();
        writer.write_alignment_record(header, &record).unwrap();
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

fn unpacked(path: &Path) -> Vec<u8> {
    zstd::decode_all(File::open(path).unwrap()).unwrap()
}

fn words(path: &Path) -> Vec<u32> {
    unpacked(path)
        .as_chunks::<4>()
        .0
        .iter()
        .map(|word| u32::from_le_bytes(*word))
        .collect()
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
    for (at, hit) in hits.iter_mut().enumerate() {
        hit.name = format!("r{at}");
    }
    hits[6].name = "r3".to_string();
    hits[6].score = -10;
    hits[0].edits = 3;
    hits[1].edits = 1;
    hits[3].others = Some("c1,11,+,50M,60,0;c2,700,-,50M,60,0;");
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

    let file = depth_file::path(&depth, "s1");
    let windows = depth_file::windows(&file, &Contigs::read(&depth).unwrap()).unwrap();
    assert_eq!(
        depth_file::read(&file, LENGTHS.len()).unwrap().reads,
        [3, 0, 3],
        "each primary alignment counts once, on its contig"
    );
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
    let summary = std::fs::read_to_string(depth.join("s1.summary.tsv")).unwrap();
    assert!(summary.contains("\ns1\t7\t6\t0\t"), "{summary}");
    assert_eq!(
        words(&depth.join("s1.links.zst")),
        [1, 2, 1],
        "a supplementary on another contig links it, one on its own does not"
    );
}

#[test]
fn kept_depth_is_reused_until_its_bam_changes() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("s1.bam");
    let depth = directory.path().join("depth");
    let samples = [Sample {
        name: "s1".to_string(),
        path: bam.clone(),
    }];
    let touch = |when: SystemTime| {
        File::options()
            .write(true)
            .open(&bam)
            .unwrap()
            .set_modified(when)
            .unwrap()
    };
    let made = |depth: &Path| {
        std::fs::metadata(depth_file::path(depth, "s1"))
            .unwrap()
            .modified()
            .unwrap()
    };
    write_bam(&bam, &hits());
    let first = read(&samples, &depth, 1).unwrap();
    let written = made(&depth);

    touch(SystemTime::now() + Duration::from_secs(3600));
    assert_eq!(read(&samples, &depth, 1).unwrap().table, first.table);
    assert_eq!(made(&depth), written, "a touched BAM is not read again");

    write_bam(&bam, &[hit(0, 0, &[(Kind::Match, 250)])]);
    touch(SystemTime::now() - Duration::from_secs(3600));
    assert_eq!(
        read(&samples, &depth, 1).unwrap().table.row(0)[0],
        1.0,
        "a BAM copied in with an older mtime is read again"
    );
}

fn tied_hits() -> Vec<Hit> {
    let (primary, secondary) = (Flags::empty(), Flags::SECONDARY);
    let mut hits = (0..14)
        .map(|at| named(&format!("u0_{at}"), primary, 0, at * 13 % 200, 100))
        .collect::<Vec<_>>();
    hits.push(named("better", primary, 2, 0, 90));
    hits.push(named("better", secondary, 0, 100, 100));
    for at in 0..40 {
        hits.push(named(&format!("u2_{at}"), primary, 2, at * 45, 100));
    }
    hits.insert(16, named("u2_0", secondary, 1, 10, 80));
    for at in 0..10 {
        let name = format!("tie_{at}");
        hits.push(named(&name, primary, 0, at * 20, 100));
        hits.push(named(&name, secondary, 2, 1000 + at * 50, 100));
    }
    hits.push(named("tie_9", Flags::SUPPLEMENTARY, 1, 0, 100));
    let (first, last) = (
        Flags::SEGMENTED | Flags::FIRST_SEGMENT,
        Flags::SEGMENTED | Flags::LAST_SEGMENT,
    );
    hits.push(named("pair", first, 2, 1500, 100));
    hits.push(named("pair", last, 0, 150, 100));
    hits.push(named("pair", last | secondary, 2, 1600, 100));
    hits.push(named("repeat", primary, 2, 1900, 100));
    hits.push(named("repeat", secondary, 2, 300, 100));
    hits
}

#[test]
fn a_bam_in_mapper_order_splits_each_tied_read_by_its_contigs_depth() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("s1.bam");
    write_with(&bam, &header(true), &tied_hits());
    let depth = directory.path().join("depth");
    let samples = [Sample {
        name: "s1".to_string(),
        path: bam,
    }];

    let table = read(&samples, &depth, 2).unwrap();

    let file = depth_file::path(&depth, "s1");
    assert_eq!(
        depth_file::read(&file, LENGTHS.len()).unwrap().reads,
        [25, 0, 43]
    );
    let means = (0..3)
        .map(|row| table.table.row(row)[0])
        .collect::<Vec<_>>();
    assert_eq!(means[1], 0.0, "a worse or supplementary hit lays no depth");
    let windows = depth_file::windows(&file, &Contigs::read(&depth).unwrap()).unwrap();
    let last = windows.len() - 2;
    assert_eq!(
        f16_value(u16::from_le_bytes([windows[last], windows[last + 1]])),
        0.25,
        "a read with two equal copies on one contig lays half its depth on each"
    );
    let (unique, tied) = (750.0, 500.0);
    let share = (means[0] * LENGTHS[0] as f64 - unique) / tied;
    let settled = means[0] / (means[0] + means[2]);
    assert!(
        (share - settled).abs() < 1e-4 && settled > 0.8,
        "c0 took {share} of the tied bases at a depth share of {settled}"
    );
    assert_eq!(words(&depth.join("s1.ties.zst")), [0, 2, 10]);
    let summary = std::fs::read_to_string(depth.join("s1.summary.tsv")).unwrap();
    assert!(summary.contains("\ns1\t68\t68\t10\t"), "{summary}");
}

#[test]
fn a_bam_in_any_other_order_splits_its_tied_reads_as_mapper_order_does() {
    let directory = tempfile::tempdir().unwrap();
    let depth_of = |name: &str, header: &sam::Header, hits: &[Hit]| {
        let bam = directory.path().join(format!("{name}.bam"));
        write_with(&bam, header, hits);
        let depth = directory.path().join(name);
        let samples = [Sample {
            name: "s1".to_string(),
            path: bam,
        }];
        let table = read(&samples, &depth, 2).unwrap().table;
        let file = depth_file::path(&depth, "s1");
        let reads = depth_file::read(&file, LENGTHS.len()).unwrap().reads;
        (table, reads, words(&depth.join("s1.ties.zst")))
    };
    let (table, reads, ties) = depth_of("mapper", &header(true), &tied_hits());
    let mut sorted = tied_hits();
    sorted.sort_by_key(|hit| (hit.contig, hit.start));
    let mut reversed = tied_hits();
    reversed.reverse();
    for (name, hits) in [("sorted", sorted), ("reversed", reversed)] {
        let (other, other_reads, other_ties) = depth_of(name, &header(false), &hits);
        assert_eq!((&other_reads, &other_ties), (&reads, &ties), "{name}");
        assert!(
            table
                .iter()
                .zip(&other)
                .all(|(want, got)| (want - got).abs() < 1e-9),
            "{name} reads {other} where mapper order reads {table}"
        );
    }
}

#[test]
fn a_bam_whose_read_is_split_apart_is_refused() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("s1.bam");
    write_with(
        &bam,
        &header(true),
        &[
            named("a", Flags::empty(), 0, 0, 100),
            named("b", Flags::empty(), 2, 0, 100),
            named("a", Flags::SECONDARY, 2, 100, 100),
        ],
    );
    let samples = [Sample {
        name: "s1".to_string(),
        path: bam,
    }];
    let Err(error) = read(&samples, &directory.path().join("depth"), 1) else {
        panic!("a read whose hits are apart was read");
    };
    assert!(
        format!("{error:#}").contains("not next to each other"),
        "{error:#}"
    );
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
    let run = &recover.run;
    let source = &run.reads.coverage;
    let coverage = calculate_coverage(
        &CoverageInputs {
            assembly: &run.assembly,
            output_directory: &run.common.output_directory,
            depth: &made,
            threads: 1,
            samples: &source.samples,
            mapping: &run.reads.mapping,
        },
        source.coverage_file.as_ref(),
        source.depth_directory.is_some(),
    )
    .unwrap();
    assert_eq!(coverage.windows.as_deref(), Some(made.as_path()));
    assert_eq!(coverage.table.sample_names, ["s1", "s2"]);
    let exported = CoverageTable::from_any_file(made.join(depth_file::TABLE)).unwrap();
    assert_eq!(exported.sample_names, ["s1", "s2"]);
    assert_eq!(exported.table, coverage.table.table);
    assert_eq!(coverage.table.table.row(0)[2], 1.0);
}
