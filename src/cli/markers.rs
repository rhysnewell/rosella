use clap::Args;

use crate::cli::runtime::{percentage, unit_interval};

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Single copy markers")]
pub struct MarkerParams {
    /// Pieces the protein file is cut into, each searched by its own hmmsearch. The default is
    /// a quarter of the thread count, since a shard costs a master thread plus its workers
    #[arg(long = "hmm-shards", value_parser = clap::value_parser!(u16).range(1..=64),
          hide_short_help = true)]
    pub hmm_shards: Option<u16>,

    /// Least share of a model a cut gene must span before its hit is trusted
    #[arg(long = "marker-fragment-span", default_value_t = crate::markers::fragments::DEFAULT_SPAN,
          value_parser = unit_interval, hide_short_help = true)]
    pub marker_fragment_span: f64,

    /// Reuse the single copy marker annotation across runs over the same assembly, keyed on
    /// the build and every setting that changes it
    #[arg(long = "marker-cache", hide_short_help = true)]
    pub marker_cache: Option<String>,

    /// How much lower the marker bar sits than the requested completeness
    #[arg(long = "marker-bar-offset", default_value_t = crate::markers::DEFAULT_BAR_OFFSET,
          value_parser = percentage, hide_short_help = true)]
    pub marker_bar_offset: f64,
}
