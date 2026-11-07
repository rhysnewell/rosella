use rosella::coverage::bridges::{Piece, bridges};

const LENGTH: usize = 10_000;

#[test]
fn a_bridge_needs_the_read_to_leave_one_contig_end_and_enter_another() {
    let cases = [
        (
            "leaves c0 at its end, enters c1 at its start",
            Piece::parsed(7_000, "3000M1000S", false),
            Piece::parsed(0, "3000S1000M", false),
            true,
        ),
        (
            "the second contig is read reversed",
            Piece::parsed(7_000, "3000M1000S", false),
            Piece::parsed(LENGTH - 1_000, "1000M3000S", true),
            true,
        ),
        (
            "the supplementary lands mid-contig",
            Piece::parsed(7_000, "3000M1000S", false),
            Piece::parsed(4_000, "3000S1000M", false),
            false,
        ),
        (
            "both cover one stretch of the read",
            Piece::parsed(7_000, "3000M1000S", false),
            Piece::parsed(0, "1000M3000S", false),
            false,
        ),
        (
            "the clip is an aligner trim",
            Piece::parsed(7_000, "3000M50S", false),
            Piece::parsed(0, "3000S50M", false),
            false,
        ),
    ];
    for (case, one, other, expected) in cases {
        let (one, other) = (one.unwrap(), other.unwrap());
        assert_eq!(bridges(&one, LENGTH, &other, LENGTH), expected, "{case}");
        assert_eq!(
            bridges(&other, LENGTH, &one, LENGTH),
            expected,
            "{case}, read the other way"
        );
    }
}

#[test]
fn a_cigar_with_an_unknown_op_is_refused() {
    assert!(Piece::parsed(0, "10M5Q", false).is_none());
}
