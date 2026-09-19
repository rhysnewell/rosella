use std::hash::Hasher;

const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const PRIME: u64 = 0x0000_0100_0000_01b3;

/// The hasher in std is seeded per process and free to change between compiler versions. Two
/// runs of the same binary have to agree. A rebuild must not rename every cache entry.
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
