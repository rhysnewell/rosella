use std::collections::BTreeMap;
use std::fs::File;
use std::path::Path;

use noodles_bam as bam;
use noodles_sam::alignment::record::data::field::Tag;
use rosella::coverage::coverage_calculator::Reads;
use rosella::external::provenance::seeded_alike;
use rosella::external::strobealign::{index_bytes, map, parts, settings};
use strobealign::seeding::SeedingParameters;

const DATA: &str = "tests/data/minimap2";

type Hits = BTreeMap<(String, bool), (usize, i64, Vec<(usize, usize)>)>;

fn hits(path: &Path) -> Hits {
    let mut reader = bam::io::Reader::new(File::open(path).unwrap());
    reader.read_header().unwrap();
    let mut found = Hits::new();
    for record in reader.records() {
        let record = record.unwrap();
        let flags = record.flags();
        let name = String::from_utf8(record.name().unwrap().to_vec()).unwrap();
        let mate =
            found
                .entry((name, flags.is_last_segment()))
                .or_insert((0, i64::MIN, Vec::new()));
        if flags.is_unmapped() {
            continue;
        }
        mate.0 += usize::from(!flags.is_secondary());
        let score = record
            .data()
            .get(&Tag::ALIGNMENT_SCORE)
            .unwrap()
            .unwrap()
            .as_int()
            .unwrap();
        let at = (
            record.reference_sequence_id().unwrap().unwrap(),
            usize::from(record.alignment_start().unwrap().unwrap()),
        );
        if score > mate.1 {
            (mate.1, mate.2) = (score, Vec::new());
        }
        if score == mate.1 {
            mate.2.push(at);
            mate.2.sort_unstable();
            mate.2.dedup();
        }
    }
    found
}

#[test]
fn an_assembly_mapped_in_parts_keeps_each_reads_best_hits_and_one_primary() {
    let home = tempfile::tempdir().unwrap();
    let set = Reads::Paired(
        format!("{DATA}/reads_1.fq.gz"),
        format!("{DATA}/reads_2.fq.gz"),
    );
    let assembly = format!("{DATA}/assembly.fa");
    let (whole, split) = (home.path().join("whole.bam"), home.path().join("split.bam"));
    map(&assembly, &set, 2, 5, &whole, u64::MAX).unwrap();
    map(&assembly, &set, 2, 5, &split, 1).unwrap();
    let header = bam::io::Reader::new(File::open(&split).unwrap())
        .read_header()
        .unwrap();
    assert!(seeded_alike(&header, &settings(&set, 5).unwrap()));
    assert!(!seeded_alike(&header, &settings(&set, 20).unwrap()));

    let (whole, split) = (hits(&whole), hits(&split));
    assert_eq!(whole.len(), 600);
    assert!(whole.values().filter(|mate| mate.0 == 1).count() > 500);
    for (mate, (primaries, _, _)) in &split {
        assert!(*primaries <= 1, "{mate:?} has {primaries} primaries");
    }
    let kept = |(score, hits): (i64, &[(usize, usize)]), split: Option<&(usize, i64, Vec<_>)>| {
        split.is_some_and(|(_, best, found)| {
            *best == score && hits.iter().all(|hit| found.contains(hit))
        })
    };
    let differ = whole
        .iter()
        .filter(|(mate, (_, score, hits))| !kept((*score, hits), split.get(*mate)))
        .map(|(mate, hits)| format!("{mate:?} whole {hits:?} split {:?}", split.get(mate)))
        .collect::<Vec<_>>();
    assert!(differ.is_empty(), "{}", differ.join("\n"));
}

#[test]
fn a_contig_longer_than_a_part_gets_a_part_of_its_own() {
    let bytes = index_bytes(&SeedingParameters::new(150));
    let lengths = [100, 5_000_000, 200, 300, 400, 1_000_000];
    let budget = bytes(lengths.iter().sum()).div_ceil(2);
    assert_eq!(
        parts(&lengths, &bytes, budget),
        [vec![1], vec![0, 2, 3, 4, 5]]
    );
}
