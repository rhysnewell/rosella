use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Result, bail};
use clap::{Args, ValueEnum};
use log::warn;

use crate::cli::runtime::unit_interval;
use crate::markers::genome_tree::Data;
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

    #[command(flatten)]
    pub data: Databases,
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
    pub fn prepare(&self, data: Option<&Located>) -> Result<Option<Arc<Panel>>> {
        if self.auxiliary_markers.is_empty() {
            return Ok(None);
        }
        let Some(data) = data else {
            bail!("the BUSCO sets are not downloaded, so run once with --download-data <DIR>");
        };
        let chosen = LINEAGES
            .into_iter()
            .filter(|lineage| {
                self.auxiliary_markers
                    .iter()
                    .any(|choice| choice.holds(lineage))
            })
            .collect::<Vec<_>>();
        let panel = Panel::prepare(
            &data.directory.join("auxiliary_markers"),
            &chosen,
            data.fetch,
        )?;
        Ok(Some(Arc::new(panel)))
    }
}

#[derive(Args, Debug, Clone)]
pub struct Databases {
    #[arg(
        long = "download-data",
        value_name = "DIR",
        help_heading = "Single copy markers",
        help = "Download CheckM's lineage sets and any BUSCO sets asked for into DIR, and remember DIR for later runs"
    )]
    pub download_data: Option<PathBuf>,
}

pub struct Located {
    pub directory: PathBuf,
    pub fetch: bool,
}

impl Databases {
    pub fn locate(&self) -> Result<Option<Located>> {
        let Some(directory) = &self.download_data else {
            return Ok(crate::data_dir::remembered().map(|directory| Located {
                directory,
                fetch: false,
            }));
        };
        std::fs::create_dir_all(directory)?;
        crate::data_dir::remember(directory)?;
        Ok(Some(Located {
            directory: directory.clone(),
            fetch: true,
        }))
    }
}

// Without the downloaded sets CheckM still reports, on the domain sets rosella ships.
pub fn checkm_lineage(data: Option<&Located>) -> Option<Arc<Data>> {
    let data = data?;
    match Data::prepare(&data.directory.join("checkm_lineage"), data.fetch) {
        Ok(found) => found.map(Arc::new),
        Err(error) => {
            warn!("CheckM lineage sets: {error:#}, so CheckM reads the domain sets");
            None
        }
    }
}
