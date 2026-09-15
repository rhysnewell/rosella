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

    /// What the pool keeps out of the pot: bins over every bar, bins at genome scale, those
    /// also under the tier's contamination, or those also over the completeness bar
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

    /// Completeness bar of the pool's last rung, as a share of the full bar
    #[arg(long = "rung-floor", default_value_t = crate::refine::rung::DEFAULT_RUNG_FLOOR,
          value_parser = unit_interval, hide_short_help = true)]
    pub rung_floor: f64,

    /// Searches of the pool, each one over the whole of it, with the neighbour count halving
    /// each round so a genome the dense graph buries can still form its own community
    #[arg(long = "dissolve-rounds", default_value_t = 6,
          value_parser = clap::value_parser!(u16).range(1..=8), hide_short_help = true)]
    pub dissolve_rounds: u16,

    /// Cap on the passes over the pool, each one re-embedding what the pass before it left
    /// unclaimed. The passes stop on their own once one finds bins the model scores worse
    /// than the last
    #[arg(long = "dissolve-passes", default_value_t = 3,
          value_parser = clap::value_parser!(u16).range(1..=32), hide_short_help = true)]
    pub dissolve_passes: u16,

    /// Partition seeds the ladder is built at. Every labelling from every seed reaches the
    /// per-bin combination
    #[arg(long = "partition-seeds", default_value_t = 3,
          value_parser = clap::value_parser!(u16).range(1..=16), hide_short_help = true)]
    pub partition_seeds: u16,

    /// Weight on contamination when ranking rescue candidates by worth
    #[arg(long = "worth-contamination",
          default_value_t = crate::refine::rung::DEFAULT_WORTH_CONTAMINATION,
          hide_short_help = true)]
    pub worth_contamination: f64,

    /// Contamination a bin may carry before worth charges it any
    #[arg(long = "worth-allowance", default_value_t = 0.0, value_parser = percentage,
          hide_short_help = true)]
    pub worth_allowance: f64,

    /// Completeness points below the bar a bin may sit and still draw contigs from the bins
    /// under --min-bin-size, which are discarded anyway
    #[arg(long = "recruit-near-bar", value_parser = percentage, hide_short_help = true)]
    pub recruit_near_bar: Option<f64>,

    /// Put a dissolved bin back whole when the pool broke it into pieces that all miss the bar
    #[arg(long = "dissolve-restore", hide_short_help = true)]
    pub dissolve_restore: bool,

    /// Rungs the pool walks past its last one, each loosening contamination further
    #[arg(long = "dissolve-extra-rungs", default_value_t = 0,
          value_parser = clap::value_parser!(u16).range(0..=6), hide_short_help = true)]
    pub dissolve_extra_rungs: u16,
}
