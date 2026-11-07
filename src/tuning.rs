// Thresholds no sweep has ever moved, in one file so the unmeasured set is one edit. One with a
// derivation, a citation or a standard behind it stays beside the code that reads it.

pub const MISPLACED_LENGTH: usize = 1_000_000;

// All pairs up to here. Past it every contig is scored against one shared sample, since sampling
// pairs would leave most contigs with no estimate at all.
pub const EXACT_LIMIT: usize = 2_000;

pub const REFERENCE_SAMPLE: usize = 1_000;

// Density validity says a labelling separates well, not that the bin was chimeric, and a pure
// genome separates happily. Without this margin the split takes good bins apart.
pub const REQUIRED_IMPROVEMENT: f64 = 0.9;

// Noise above this fraction of the original bin means the split threw away more than it
// explained.
pub const MAX_NOISE_FRACTION: f64 = 0.6;

pub const FAMILY_ALPHA: f64 = 0.05;

pub const BISECT_SIZE_MULTIPLE: usize = 2;

// Two candidates a hundredth of a bin apart are the same proposal to the bar, so a lineage
// only re-enters the heap once it has grown enough to be a different answer.
pub const LINEAGE_GROWTH: f64 = 1.01;

pub const SWEEP_WIDTH: usize = 10;

pub const CONVERGENCE_FRACTION: f64 = 0.01;

// The short contigs' own kNN sets only the share's denominator, and half the main k costs no tier.
pub const ATTACH_NEIGHBOURS: usize = 50;

// Floor on a contig's per-sample coverage variance, so a contig reported with none does not
// read as infinitely certain about its depth.
pub const MIN_VAR: f64 = 1.0;

pub const PEEL_SIGMA: f64 = 1.0;

pub const LADDER_COARSEST: f64 = 2.0;
pub const LADDER_FINEST: f64 = 512.0;

// A sample size, not an all-pairs threshold.
pub const DIP_SAMPLE: usize = 2_000;

// Depth below this share of a pair's deepest sample counts as absent, and that sample leaves
// the pair's coverage distance. At one sample mutual absence otherwise reads as agreement.
pub const PRESENCE_FRACTION: f64 = 0.01;

// Below this many observed markers a bin carries no evidence about which lineage it is, so
// the widest set is used rather than letting a handful of absences pick a reduced one.
pub const MARKERS_TO_CHOOSE_A_SET: usize = 10;
