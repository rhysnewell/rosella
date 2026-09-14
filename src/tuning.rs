//! Thresholds that decide a bin and that no sweep has ever moved. They live together so the
//! unmeasured set is one file rather than twenty numbers spread across fourteen modules, and
//! so a sweep over any of them is one edit. A threshold with a derivation, a citation or a
//! standard behind it stays beside the code that reads it.

/// Bases of out of place contigs that force a split on their own.
pub const MISPLACED_LENGTH: usize = 1_000_000;

/// All pairs up to here. Past it every contig is scored against one shared sample instead,
/// which keeps the per-contig figures usable where sampling pairs would leave most contigs
/// with no estimate at all.
pub const EXACT_LIMIT: usize = 2_000;

/// Size of that shared sample.
pub const REFERENCE_SAMPLE: usize = 1_000;

/// The pieces have to be this much tighter than the bin they came out of. Density validity
/// says a labelling separates well, not that the bin was chimeric, and a pure genome
/// separates perfectly happily. Without this the split takes good bins apart.
pub const REQUIRED_IMPROVEMENT: f64 = 0.9;

/// Noise above this fraction of the original bin means the split threw away more than it
/// explained.
pub const MAX_NOISE_FRACTION: f64 = 0.6;

/// Aggregate distance a leftover piece must hold within before it is kept as its own bin.
pub const LEFTOVER_AGGREGATE: f64 = 0.5;

/// Family-wise level the bisect dip test is corrected to.
pub const FAMILY_ALPHA: f64 = 0.05;

/// Multiple of the output floor a bin must reach before bisecting it is worth the work.
pub const BISECT_SIZE_MULTIPLE: usize = 2;

/// Two candidates a hundredth of a bin apart are the same proposal to the bar, so a lineage
/// only re-enters the heap once it has grown enough to be a different answer.
pub const LINEAGE_GROWTH: f64 = 1.01;

/// A pair joins, then the pair it made can take a third piece, but the chain is short and
/// every pass costs a full sweep of the boosters.
pub const JOIN_PASSES: usize = 4;

/// Rungs in the resolution ladder every partition is drawn from.
pub const SWEEP_WIDTH: usize = 10;

/// Share of the rows still moving below which the neighbour descent has converged.
pub const CONVERGENCE_FRACTION: f64 = 0.001;

/// Floor on a contig's per-sample coverage variance, so a contig reported with none does not
/// read as infinitely certain about its depth.
pub const MIN_VAR: f64 = 1.0;

/// Bit score floor hmmsearch is given when rescuing a marker cut by a contig end, where the
/// model's own gathering cutoff is out of reach by construction.
pub const DOMAIN_FLOOR: &str = "10";

/// Standard deviations past its bin's mean a contig has to sit before the peel takes it.
pub const PEEL_SIGMA: f64 = 1.0;

/// Coarsest and finest community the resolution ladder spans, as a divisor of the graph's
/// total node mass.
pub const LADDER_COARSEST: f64 = 2.0;
pub const LADDER_FINEST: f64 = 512.0;

/// The rescue ladder's size floor falls by this much a rung and stops here, so the later
/// rungs relax the bar without also letting a smaller bin through.
pub const RUNG_FLOOR_STEP: f64 = 0.25;
pub const RUNG_FLOOR_FLOOR: f64 = 0.5;
