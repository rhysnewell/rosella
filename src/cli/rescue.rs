use clap::{ArgAction, Args, ValueEnum};

use crate::cli::runtime::percentage;

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Rescue pool")]
pub struct RescueParams {
    /// Pool every bin short of the bars with the unbinned contigs and bin them again
    #[arg(long = "dissolve", value_enum, default_value_t = Switch::On)]
    pub dissolve: Switch,

    /// Least completeness for the pool to keep a candidate
    #[arg(long = "min-completeness", default_value_t = crate::refine::rung::DEFAULT_COMPLETENESS,
          value_parser = percentage, hide_short_help = true)]
    pub min_completeness: f64,

    /// Most contamination the pool lets a candidate carry
    #[arg(long = "max-contamination", default_value_t = crate::refine::rung::DEFAULT_CONTAMINATION,
          value_parser = percentage, hide_short_help = true)]
    pub max_contamination: f64,

    /// Neighbour widths the pool searches at. Each round halves the width, so a genome the
    /// dense graph buries can form its own community. One graph is built per pass
    #[arg(long = "dissolve-rounds", default_value_t = 6,
          value_parser = clap::value_parser!(u16).range(1..=8), hide_short_help = true)]
    pub dissolve_rounds: u16,

    /// Most passes over the pool. Passes stop early once one scores worse than the last
    #[arg(long = "dissolve-passes", default_value_t = 3,
          value_parser = clap::value_parser!(u16).range(1..=32), hide_short_help = true)]
    pub dissolve_passes: u16,

    /// Weight on contamination when ranking rescue candidates by worth
    #[arg(long = "worth-contamination",
          default_value_t = crate::refine::rung::DEFAULT_WORTH_CONTAMINATION,
          hide_short_help = true)]
    pub worth_contamination: f64,

    /// Skip the splitter that cuts and peels chimeric bins
    #[arg(long = "no-refine", action = ArgAction::SetTrue, help_heading = "Refinement")]
    pub no_refine: bool,

    /// Let the scorer fuse whole bin pairs. Off by default: the markers cannot tell two halves
    /// of one genome from two genomes
    #[arg(long = "join", action = ArgAction::SetTrue, help_heading = "Refinement")]
    pub join: bool,

    /// Order of the refine cycle's stages. Name a stage twice to run it twice, or leave it out
    /// to skip it
    #[arg(long = "stage-order", default_value = crate::recover::recover_engine::SHIPPED_ORDER,
          help_heading = "Refinement", hide_short_help = true)]
    pub stage_order: String,

    /// Let refinement cut a bin into the communities inside it
    #[arg(long = "split", value_enum, default_value_t = Switch::On,
          help_heading = "Refinement", hide_short_help = true)]
    pub split: Switch,

    /// Let refinement peel lone contigs off a bin
    #[arg(long = "peel", value_enum, default_value_t = Switch::On,
          help_heading = "Refinement", hide_short_help = true)]
    pub peel: Switch,

    /// Take marker-free passengers out of their bins before writing
    #[arg(long = "passengers", value_enum, default_value_t = Switch::On,
          help_heading = "Refinement", hide_short_help = true)]
    pub passengers: Switch,

    /// Write small replicons as bins of their own
    #[arg(long = "replicons", value_enum, default_value_t = Switch::On,
          help_heading = "Refinement", hide_short_help = true)]
    pub replicons: Switch,

    #[arg(long = "strays", value_enum, default_value_t = Switch::Off, hide = true)]
    pub strays: Switch,

    #[arg(long = "adrift", value_enum, default_value_t = Switch::Off, hide = true)]
    pub adrift: Switch,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Switch {
    On,
    Off,
}
