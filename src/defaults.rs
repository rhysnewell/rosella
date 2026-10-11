pub const MIN_CONTIG_SIZE: usize = 1500;
pub const THREADS: usize = 10;
pub const FASTA_EXTENSION: &str = "fna";
pub const QUALITY_FILE: &str = "quality.tsv";

// Forty epochs read 338/413/492/516/543 on mouse gut 16 against 332/394/478/507/537 at twenty.
pub const EPOCHS: u32 = 40;

// Two flat epochs stopped mouse gut 16 at 27 of 40 and cost 14 t1, so the stop waits for five.
pub const PATIENCE: u32 = 5;

// Five and seven stretches a step read level with three on mouse gut, for 40 and 100 per cent
// more training.
pub const STRETCHES: u32 = 3;

// The assembly is read a chunk at a time so it is never held whole beside what is built from it,
// which was the memory peak on a multi-sample assembly.
pub const CHUNK_BASES: usize = 64 << 20;

// Far enough from the `+ 1` the visit order shuffles on that a sampler seeded from it is not
// replaying the same stream one node later.
pub const SEED_STRIDE: u64 = 0x9E37_79B9_7F4A_7C15;

// The partition ensemble and the neighbour search draw from these rather than from --seed, so
// two runs of one assembly agree and nobody has to pick a seed or run three.
pub const PARTITION_SEED: u64 = 42;
pub const KNN_SEED: u64 = 42;
// minimap2 -x sr mapped every set rosella was measured on, and map-ont is the long read preset
// that also takes PacBio reads.
pub const SHORT_READ_MAPPER: &str = "minimap2-sr";
pub const LONG_READ_MAPPER: &str = "minimap2-ont";
