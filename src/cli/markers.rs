use clap::Args;

use crate::cli::runtime::unit_interval;

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Single copy markers")]
pub struct MarkerParams {
    /// Least share of a model a cut gene must span before its hit is trusted
    #[arg(long = "marker-fragment-span", default_value_t = crate::markers::fragments::DEFAULT_SPAN,
          value_parser = unit_interval, hide_short_help = true)]
    pub marker_fragment_span: f64,

    /// Reuse the single copy marker annotation across runs over the same assembly, keyed on
    /// the build and every setting that changes it
    #[arg(long = "marker-cache", hide_short_help = true)]
    pub marker_cache: Option<String>,

    /// Also score each bin on CheckM1's lineage marker sets in quality.tsv, searched over the
    /// binned contigs once the bins are known
    #[arg(long = "checkm")]
    pub checkm: bool,
}
