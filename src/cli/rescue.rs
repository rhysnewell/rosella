use clap::{ArgAction, Args};

use crate::cli::runtime::percentage;

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Rescue pool")]
pub struct RescueParams {
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

    #[arg(long = "strain-split", action = ArgAction::SetTrue, hide = true)]
    pub strain_split: bool,

    /// Skip the splitter that cuts and peels chimeric bins
    #[arg(long = "no-refine", action = ArgAction::SetTrue, help_heading = "Refinement")]
    pub no_refine: bool,
}
