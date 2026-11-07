use rosella::recover::circular::{TERMINAL_REPEAT, closes, marked};

#[test]
fn a_header_flag_names_a_circle_only_as_a_true_value_of_its_own_field() {
    for (header, circle) in [
        ("1 length=4200 depth=1.00x circular=true", true),
        ("ctg7 length=88123 coverage=12 circular=yes", true),
        ("ctg7 Circular:Y", true),
        ("ctg7 circular=no", false),
        ("circular=true", false),
        ("ctg7 notcircular=true", false),
    ] {
        assert_eq!(marked(header.as_bytes()), circle, "{header}");
    }
}

#[test]
fn a_terminal_repeat_closes_a_contig_only_from_twenty_clean_bases() {
    let body = b"GATTACAGGCCTTAAGCTAGCTTGACCATGCAAGTCGATCGGATCCATG".to_vec();
    let closing = |repeat: usize| [&body[..], &body[..repeat]].concat();
    let mut masked = closing(TERMINAL_REPEAT);
    masked[..TERMINAL_REPEAT].make_ascii_lowercase();
    let mut gapped = closing(TERMINAL_REPEAT + 5);
    gapped[3] = b'N';
    let last = gapped.len() - (TERMINAL_REPEAT + 5) + 3;
    gapped[last] = b'N';
    for (sequence, circle, case) in [
        (closing(TERMINAL_REPEAT), true, "the shortest repeat"),
        (closing(TERMINAL_REPEAT + 9), true, "a longer repeat"),
        (masked, true, "a soft-masked start"),
        (closing(TERMINAL_REPEAT - 1), false, "one base short"),
        (body.clone(), false, "no repeat"),
        (gapped, false, "a gap in the start"),
    ] {
        assert_eq!(closes(&sequence), circle, "{case}");
    }
}
