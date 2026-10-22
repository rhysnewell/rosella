use clap::Args;

use crate::cli::runtime::unit_interval;

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Single copy markers")]
pub struct MarkerParams {
    /// Least share of a marker model that a gene cut by a contig end must cover to count
    #[arg(long = "marker-fragment-span", default_value_t = crate::markers::fragments::DEFAULT_SPAN,
          value_parser = unit_interval, hide_short_help = true)]
    pub marker_fragment_span: f64,

    /// Directory to keep the marker annotation in. A later run over the same assembly skips
    /// gene calling
    #[arg(long = "marker-cache", hide_short_help = true)]
    pub marker_cache: Option<String>,

    /// Also score each bin on CheckM1's lineage marker sets in quality.tsv
    #[arg(long = "checkm")]
    pub checkm: bool,
}
