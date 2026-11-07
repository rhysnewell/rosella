use noodles_sam::alignment::record::cigar::op::Kind;

// Aligners trim a few bases off a read's end, so an exit is a read leaving within this of the
// contig's end. 50 and 200 gave the same bridges on CAMI III gold assemblies.
const END_SLOP: usize = 200;

pub struct Piece {
    start: usize,
    end: usize,
    clips: [usize; 2],
    read: [usize; 2],
    reverse: bool,
}

impl Piece {
    pub fn of(start: usize, ops: impl IntoIterator<Item = (Kind, usize)>, reverse: bool) -> Self {
        let ops = ops.into_iter().collect::<Vec<_>>();
        let clip = |op: Option<&(Kind, usize)>| match op {
            Some((Kind::SoftClip | Kind::HardClip, length)) => *length,
            _ => 0,
        };
        let clips = [clip(ops.first()), clip(ops.last())];
        let span = ops
            .iter()
            .filter(|(kind, _)| kind.consumes_reference())
            .map(|(_, length)| length)
            .sum::<usize>();
        let query = ops
            .iter()
            .filter(|(kind, _)| kind.consumes_read() && *kind != Kind::SoftClip)
            .map(|(_, length)| length)
            .sum::<usize>();
        let first = if reverse { clips[1] } else { clips[0] };
        Self {
            start,
            end: start + span,
            clips,
            read: [first, first + query],
            reverse,
        }
    }

    pub fn parsed(start: usize, cigar: &str, reverse: bool) -> Option<Self> {
        let mut ops = Vec::new();
        let mut length = 0usize;
        for byte in cigar.bytes() {
            match byte {
                b'0'..=b'9' => length = length * 10 + usize::from(byte - b'0'),
                _ => {
                    ops.push((kind(byte)?, length));
                    length = 0;
                }
            }
        }
        Some(Self::of(start, ops, reverse))
    }

    fn exits(&self, length: usize) -> [bool; 2] {
        let left = self.start <= END_SLOP && self.clips[0] > END_SLOP;
        let right = length.saturating_sub(self.end) <= END_SLOP && self.clips[1] > END_SLOP;
        match self.reverse {
            false => [left, right],
            true => [right, left],
        }
    }
}

// Repeats and sequence two strains share land mid-contig or align one stretch of a read twice.
pub fn bridges(one: &Piece, one_length: usize, other: &Piece, other_length: usize) -> bool {
    let (first, second, first_length, second_length) = match one.read[0] <= other.read[0] {
        true => (one, other, one_length, other_length),
        false => (other, one, other_length, one_length),
    };
    first.read[1] <= second.read[0] + END_SLOP
        && first.exits(first_length)[1]
        && second.exits(second_length)[0]
}

fn kind(byte: u8) -> Option<Kind> {
    Some(match byte {
        b'M' => Kind::Match,
        b'I' => Kind::Insertion,
        b'D' => Kind::Deletion,
        b'N' => Kind::Skip,
        b'S' => Kind::SoftClip,
        b'H' => Kind::HardClip,
        b'P' => Kind::Pad,
        b'=' => Kind::SequenceMatch,
        b'X' => Kind::SequenceMismatch,
        _ => return None,
    })
}
