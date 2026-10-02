use clap::{Args, ValueEnum};

use crate::cli::runtime::percentage;

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Rescue pool")]
pub struct RescueParams {
    /// Whether every bin short of the bars goes back in the pot with the unbinned and is
    /// embedded again as one pool
    #[arg(long = "dissolve", value_enum, default_value_t = Switch::On)]
    pub dissolve: Switch,

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

    /// Partition seeds the ladder is built at. Each seed and arm contributes its best rung
    /// to the per-bin combination, not every labelling it made
    #[arg(long = "partition-seeds", default_value_t = 3,
          value_parser = clap::value_parser!(u16).range(1..=16), hide_short_help = true)]
    pub partition_seeds: u16,

    /// Weight on contamination when ranking rescue candidates by worth
    #[arg(long = "worth-contamination",
          default_value_t = crate::refine::rung::DEFAULT_WORTH_CONTAMINATION,
          hide_short_help = true)]
    pub worth_contamination: f64,

    /// Order the refine cycle runs its stages in. A stage may be named twice to run it twice,
    /// or left out to skip it
    #[arg(long = "stage-order", default_value = crate::recover::recover_engine::SHIPPED_ORDER,
          help_heading = "Refinement", hide_short_help = true)]
    pub stage_order: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Switch {
    On,
    Off,
}
