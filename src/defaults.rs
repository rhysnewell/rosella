pub const MIN_CONTIG_SIZE: usize = 1500;
pub const THREADS: usize = 10;
pub const FASTA_EXTENSION: &str = "fna";
pub const QUALITY_FILE: &str = "quality.tsv";

/// Far enough from the `+ 1` the visit order shuffles on that a sampler seeded from it is not
/// replaying the same stream one node later.
pub const SEED_STRIDE: u64 = 0x9E37_79B9_7F4A_7C15;

/// The partition ensemble and the neighbour search draw from these rather than from --seed, so
/// two runs of one assembly agree and nobody has to pick a seed or run three.
pub const PARTITION_SEED: u64 = 42;
pub const KNN_SEED: u64 = 42;
