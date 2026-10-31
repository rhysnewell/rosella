// Direct terminal repeats of 20 bp or more mark a closed genome (Nayfach et al. 2021).
pub const TERMINAL_REPEAT: usize = 20;

const TRUE: [&str; 3] = ["true", "yes", "y"];

pub fn marked(header: &[u8]) -> bool {
    let Ok(header) = std::str::from_utf8(header) else {
        return false;
    };
    header.split_whitespace().skip(1).any(|field| {
        field.split_once(['=', ':']).is_some_and(|(key, value)| {
            key.eq_ignore_ascii_case("circular")
                && TRUE.iter().any(|truth| value.eq_ignore_ascii_case(truth))
        })
    })
}

// An assembler that closes a circle keeps the overlap it closed on at both ends.
pub fn closes(sequence: &[u8]) -> bool {
    let length = sequence.len();
    if length < 2 * TERMINAL_REPEAT
        || !sequence[..TERMINAL_REPEAT]
            .iter()
            .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    {
        return false;
    }
    (length / 2..=length - TERMINAL_REPEAT)
        .any(|start| sequence[start..].eq_ignore_ascii_case(&sequence[..length - start]))
}
