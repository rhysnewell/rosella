use clap::Args;

use crate::cli::runtime::{percentage, unit_interval};

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Rescue pool")]
pub struct RescueParams {
    /// Whether every bin short of the bars goes back in the pot with the unbinned and is
    /// embedded again as one pool
    #[arg(long = "dissolve", value_parser = crate::recover::settings::DISSOLVE_NAMES,
          default_value = "on")]
    pub dissolve: String,

    /// What the pool keeps out of the pot: bins over every bar, bins at genome scale, or those
    /// also under the tier's contamination
    #[arg(long = "dissolve-hold", value_parser = crate::recover::settings::HOLD_NAMES,
          default_value = "bars", hide_short_help = true)]
    pub dissolve_hold: String,

    /// Completeness a candidate needs before the pool adopts it
    #[arg(long = "min-completeness", default_value_t = crate::refine::rung::DEFAULT_COMPLETENESS,
          value_parser = percentage, hide_short_help = true)]
    pub min_completeness: f64,

    /// Contamination a candidate may carry before the pool refuses it
    #[arg(long = "max-contamination", default_value_t = crate::refine::rung::DEFAULT_CONTAMINATION,
          value_parser = percentage, hide_short_help = true)]
    pub max_contamination: f64,

    /// Neighbour widths the pool is searched at, halving each round so a genome the dense
    /// graph buries can still form its own community. One graph is built per pass and each
    /// round narrows it, so a round costs a partition rather than a search
    #[arg(long = "dissolve-rounds", default_value_t = 6,
          value_parser = clap::value_parser!(u16).range(1..=8), hide_short_help = true)]
    pub dissolve_rounds: u16,

    /// Cap on the passes over the pool. The passes stop on their own once one finds bins the
    /// model scores worse than the last
    #[arg(long = "dissolve-passes", default_value_t = 3,
          value_parser = clap::value_parser!(u16).range(1..=32), hide_short_help = true)]
    pub dissolve_passes: u16,

    /// Search the pool graph again on every pass rather than inducing it from the first build.
    /// An induced graph cannot reach a neighbour outside that build's top 100 and takes its
    /// width from its narrowest surviving row. Measured within one per cent on every tier for
    /// about three times the pool search, so it is off
    #[arg(long = "dissolve-reembed", action = clap::ArgAction::SetTrue, hide_short_help = true)]
    pub dissolve_reembed: bool,

    /// How far a pass walks its rungs: every rung, or stopping at the first that adopts
    #[arg(long = "dissolve-rung-walk",
          value_parser = crate::recover::settings::RUNG_WALK_NAMES,
          default_value = "walk", hide_short_help = true)]
    pub dissolve_rung_walk: String,

    /// Partition seeds the ladder is built at. Each seed and arm contributes its best rung
    /// to the per-bin combination, not every labelling it made
    #[arg(long = "partition-seeds", default_value_t = 3,
          value_parser = clap::value_parser!(u16).range(1..=16), hide_short_help = true)]
    pub partition_seeds: u16,

    /// Rank the ensemble's candidates on marker F1 and drain them in contamination tiers,
    /// so a pure candidate claims its contigs before a dirtier one is offered them
    #[arg(long = "peel", action = clap::ArgAction::SetTrue, hide_short_help = true)]
    pub peel: bool,

    /// Weight on contamination when ranking rescue candidates by worth
    #[arg(long = "worth-contamination",
          default_value_t = crate::refine::rung::DEFAULT_WORTH_CONTAMINATION,
          hide_short_help = true)]
    pub worth_contamination: f64,

    /// Let a bin short of the bars take single contigs back off its neighbours, while the
    /// pair of bins is worth more after the move than before it
    #[arg(long = "recruit", action = clap::ArgAction::SetTrue, help_heading = "Refinement")]
    pub recruit: bool,

    /// Order the refine cycle runs its stages in. A stage may be named twice to run it twice,
    /// or left out to skip it
    #[arg(long = "stage-order", default_value = crate::recover::recover_engine::SHIPPED_ORDER,
          help_heading = "Refinement", hide_short_help = true)]
    pub stage_order: String,

    /// Split a bin the shed judges fused into two rather than writing the redundant copy
    /// unbinned. Falls back to the eviction when the split finds one cloud
    #[arg(long = "shed-split", action = clap::ArgAction::SetTrue, help_heading = "Refinement",
          hide_short_help = true)]
    pub shed_split: bool,

    /// Refuse to shed a contig longer than this many marker spacings of its own set, where a
    /// spacing is the set's median genome over its marker count. 0 shreds at any length
    #[arg(long = "shed-length-multiple", default_value_t = 0.0, hide_short_help = true,
          help_heading = "Refinement")]
    pub shed_length_multiple: f64,

    /// Let the shed and the rung walk run even where most bins arrived over the bars. Only
    /// for attributing the gate, which is on by default
    #[arg(long = "no-finished-gate", action = clap::ArgAction::SetTrue,
          help_heading = "Refinement", hide_short_help = true)]
    pub no_finished_gate: bool,

    /// Rungs the rescue ladder walks. Scaffolding for pricing a shape nothing ever measured
    #[arg(long = "rungs", default_value_t = crate::refine::rung::RUNGS as u16,
          value_parser = clap::value_parser!(u16).range(1..=16), hide_short_help = true,
          help_heading = "Refinement")]
    pub rungs: u16,

    /// Cap on the multiple of --max-contamination a loose rung may reach. 2 is MIMAG's ceiling
    /// for a medium quality bin, the ladder's own is 3
    #[arg(long = "rung-contamination-cap", default_value_t = f64::INFINITY,
          hide_short_help = true, help_heading = "Refinement")]
    pub rung_contamination_cap: f64,

    /// How far the ladder's size floor falls a rung and where it stops, as shares of genome scale
    #[arg(long = "rung-floor-step", default_value_t = crate::tuning::RUNG_FLOOR_STEP,
          value_parser = unit_interval, hide_short_help = true, help_heading = "Refinement")]
    pub rung_floor_step: f64,

    #[arg(long = "rung-floor-floor", default_value_t = crate::tuning::RUNG_FLOOR_FLOOR,
          value_parser = unit_interval, hide_short_help = true, help_heading = "Refinement")]
    pub rung_floor_floor: f64,

    /// Completeness a bin needs before it may recruit, as a share of --min-completeness
    #[arg(long = "recruit-floor", default_value_t = crate::refine::recruit::DEFAULT_FLOOR,
          value_parser = unit_interval, hide_short_help = true)]
    pub recruit_floor: f64,

    /// How sure the coverage overlap and the markers together have to be that the receiving
    /// bin owns a contig before it is taken off the bin holding it
    #[arg(long = "recruit-confidence",
          default_value_t = crate::refine::recruit::DEFAULT_CONFIDENCE,
          value_parser = unit_interval, hide_short_help = true)]
    pub recruit_confidence: f64,
}
