use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result};
use clap::{Args, ValueEnum};

use crate::cli::runtime::unit_interval;
use crate::markers::panels::{LINEAGES, Panel};

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

    #[command(flatten)]
    pub auxiliary: AuxiliaryMarkers,
}

#[derive(Args, Debug, Clone)]
pub struct AuxiliaryMarkers {
    #[arg(
        long = "auxiliary-markers",
        value_enum,
        value_delimiter = ',',
        help_heading = "Single copy markers",
        help = "Also score bins on BUSCO's lineage sets, given as a comma separated list"
    )]
    pub auxiliary_markers: Vec<Auxiliary>,

    #[arg(
        long = "auxiliary-markers-dir",
        help_heading = "Single copy markers",
        hide_short_help = true,
        help = "Where the BUSCO sets are kept and fetched on first use [default: ~/.cache/rosella/auxiliary_markers]"
    )]
    pub auxiliary_markers_dir: Option<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Auxiliary {
    Fungi,
    Euk,
    Plant,
    Animal,
    All,
}

impl Auxiliary {
    fn holds(self, lineage: &str) -> bool {
        match self {
            Self::Fungi => lineage == "fungi",
            Self::Euk => lineage == "eukaryota",
            Self::Plant => lineage == "viridiplantae",
            Self::Animal => lineage == "metazoa",
            Self::All => true,
        }
    }
}

impl AuxiliaryMarkers {
    pub fn prepare(&self) -> Result<Option<Arc<Panel>>> {
        if self.auxiliary_markers.is_empty() {
            return Ok(None);
        }
        let chosen = LINEAGES
            .into_iter()
            .filter(|lineage| {
                self.auxiliary_markers
                    .iter()
                    .any(|choice| choice.holds(lineage))
            })
            .collect::<Vec<_>>();
        let directory = match &self.auxiliary_markers_dir {
            Some(directory) => PathBuf::from(directory),
            None => default_directory()?,
        };
        Ok(Some(Arc::new(Panel::prepare(&directory, &chosen)?)))
    }
}

fn default_directory() -> Result<PathBuf> {
    let cache = std::env::var_os("XDG_CACHE_HOME")
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|home| PathBuf::from(home).join(".cache")))
        .context("neither XDG_CACHE_HOME nor HOME is set, so pass --auxiliary-markers-dir")?;
    Ok(cache.join("rosella").join("auxiliary_markers"))
}
