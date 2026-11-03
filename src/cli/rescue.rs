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

    /// Let a replicon leave its bin on composition only when its depth disagrees too
    #[arg(long = "replicon-depth", value_enum, default_value_t = Switch::On,
          help_heading = "Refinement", hide = true)]
    pub replicon_depth: Switch,

    /// Count a marker on a sparse coding contig as present, never as a second copy
    #[arg(long = "sparse-presence", value_enum, default_value_t = Switch::On,
          help_heading = "Refinement", hide = true)]
    pub sparse_presence: Switch,

    /// Stages that count a marker gene cut by a contig end as a second copy
    #[arg(long = "cut-copies", value_enum, default_value_t = CutCopies::Both,
          help_heading = "Refinement", hide = true)]
    pub cut_copies: CutCopies,

    /// Judge contamination on the higher of the GTDB and CheckM sets. Searches CheckM's models
    /// before binning
    #[arg(long = "checkm-contamination", action = ArgAction::SetTrue,
          help_heading = "Refinement", hide = true)]
    pub checkm_contamination: bool,

    /// Stages that read CheckM under --checkm-contamination
    #[arg(long = "checkm-stages", value_enum, default_value_t = CheckmContamination::Both,
          help_heading = "Refinement", hide = true)]
    pub checkm_stages: CheckmContamination,

    /// Leave CheckM's contamination off a bin that holds no GTDB marker
    #[arg(long = "checkm-markerless", value_enum, default_value_t = Switch::Off,
          help_heading = "Refinement", hide = true)]
    pub checkm_markerless: Switch,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum CheckmContamination {
    Off,
    Combine,
    Shed,
    Both,
}

impl CheckmContamination {
    pub fn searched(self) -> bool {
        self != Self::Off
    }

    pub fn combine(self) -> bool {
        matches!(self, Self::Combine | Self::Both)
    }

    pub fn shed(self) -> bool {
        matches!(self, Self::Shed | Self::Both)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum CutCopies {
    None,
    Combine,
    Pool,
    Claims,
    Both,
    BothClaims,
    All,
}

impl CutCopies {
    pub fn combine(self) -> bool {
        matches!(self, Self::Combine | Self::Both | Self::BothClaims)
    }

    pub fn pool(self) -> bool {
        matches!(
            self,
            Self::Pool | Self::Claims | Self::Both | Self::BothClaims
        )
    }

    pub fn holds(self) -> bool {
        matches!(self, Self::Pool | Self::Both)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Switch {
    On,
    Off,
}
