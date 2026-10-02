use rosella::markers::hmm_table::Best;
use rosella::markers::{Hit, heterogeneity, identity, rebuild};
use rosella::quality::orfs::call_over;

const ASSEMBLY: &str = "tests/data/ben/1.fna.gz";

#[test]
fn identity_skips_gapped_ends_and_counts_a_gap_inside_as_a_mismatch() {
    assert_eq!(identity(b"--ABCD-", b"-XABCDE"), 1.0);
    assert_eq!(identity(b"AB-D", b"ABCD"), 0.75);
    assert!((identity(b"A--D", b"A-CD") - 2.0 / 3.0).abs() < 1e-12);
}

#[test]
fn heterogeneity_is_the_share_of_pairs_strictly_above_ninety_per_cent() {
    assert_eq!(heterogeneity(&[0.95, 0.9, 0.5, 1.0]), 50.0);
    assert_eq!(heterogeneity(&[]), 0.0);
}

#[test]
fn a_protein_rebuilt_from_its_hit_is_the_protein_that_was_called() {
    rosella::pool::init(2).ok();
    let mut sequences = Vec::new();
    let mut reader = needletail::parse_fastx_file(ASSEMBLY).unwrap();
    while let Some(record) = reader.next() {
        sequences.push(record.unwrap().seq().to_vec());
    }
    let best = Best {
        model: String::new(),
        score: 0.0,
        reach: Default::default(),
    };
    let (mut reverse, mut cut) = (0, 0);
    call_over(ASSEMBLY, 0..usize::MAX, |batch| {
        for orf in batch {
            let place = Hit::called(0, &orf, &best).place;
            reverse += usize::from(orf.reverse);
            cut += usize::from(orf.partial());
            assert_eq!(rebuild(&sequences[orf.contig], &place), Some(orf.protein));
        }
        Ok(())
    })
    .unwrap();

    assert!(reverse > 0 && cut > 0);
}
