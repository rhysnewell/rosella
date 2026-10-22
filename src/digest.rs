use std::hash::Hasher;

const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const PRIME: u64 = 0x0000_0100_0000_01b3;

// The hasher in std is seeded per process and free to change between compiler versions. Two
// runs of the same binary have to agree. A rebuild must not rename every cache entry.
pub struct Fnv1a(u64);

impl Default for Fnv1a {
    fn default() -> Self {
        Self(OFFSET)
    }
}

impl Hasher for Fnv1a {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for byte in bytes {
            self.0 = (self.0 ^ u64::from(*byte)).wrapping_mul(PRIME);
        }
    }
}

pub fn fold(bytes: &[u8]) -> u64 {
    let mut hash = Fnv1a::default();
    hash.write(bytes);
    hash.finish()
}

pub fn mix(value: u64) -> u64 {
    let mut mixed = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    mixed ^ (mixed >> 31)
}

// Summed over the contigs, so one set of contigs in any order reads alike. Lines are fed apart
// so a wrapped copy of a contig reads as the unwrapped one.
pub fn contig(name: &str, raw: &[u8], bases: usize) -> u64 {
    let mut sequence = flate2::Crc::new();
    for line in raw.split(|byte| *byte == b'\n') {
        sequence.update(line.strip_suffix(b"\r").unwrap_or(line));
    }
    let mut hash = Fnv1a::default();
    hash.write(name.as_bytes());
    hash.write(&sequence.sum().to_le_bytes());
    hash.write(&(bases as u64).to_le_bytes());
    mix(hash.finish())
}
