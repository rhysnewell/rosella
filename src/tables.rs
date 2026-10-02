use std::collections::HashMap;

use anyhow::{Result, bail};
use log::{debug, info};

use crate::assembly_graph::Link;
use crate::cli::binning::DistanceParams;
use crate::cli::coverage::ReadArgs;
use crate::cli::runtime::Common;
use crate::coverage::coverage_calculator::{CoverageInputs, calculate_coverage};
use crate::coverage::coverage_table::CoverageTable;
use crate::embedding::features::{ContigFeatures, FeatureSettings};
use crate::kmers::halves::Halves;
use crate::kmers::kmer_counting::{KmerFrequencyTable, kept_table, table_path};
use crate::kmers::scan::{Floors, scan};
use crate::kmers::sketch::ContigSketches;

pub struct Sources<'a> {
    pub assembly: &'a str,
    pub common: &'a Common,
    pub reads: &'a ReadArgs,
    pub distance: &'a DistanceParams,
    pub lengths: Lengths,
    pub threads: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct Lengths {
    pub min_contig: usize,
    pub composition: usize,
    pub sketch: Option<usize>,
    pub halves: Option<usize>,
}

impl Lengths {
    pub fn from(min_contig: usize, composition: usize) -> Self {
        Self {
            min_contig,
            composition,
            sketch: None,
            halves: None,
        }
    }
}

pub struct Tables {
    pub assembly: String,
    pub coverage: CoverageTable,
    pub coverage_file: String,
    pub tnf: KmerFrequencyTable,
    pub sketches: Option<ContigSketches>,
    pub halves: Halves,
    pub links: Option<Vec<Link>>,
}

impl Tables {
    // Both subcommands need the same row-aligned pair, so the guard, the stage timers and
    // the alignment checks live here rather than once each and differently.
    pub fn build(sources: &Sources<'_>) -> Result<Self> {
        if !std::path::Path::new(sources.assembly).is_file() {
            bail!("no assembly file at {}", sources.assembly);
        }
        let output_directory = &sources.common.output_directory;
        crate::bins::refuse_used(output_directory)?;
        std::fs::create_dir_all(output_directory)?;

        debug!("Calculating contig coverages.");
        let (mut coverage, coverage_file) = {
            let _timer = crate::timing::scope("coverage");
            calculate_coverage(&CoverageInputs {
                assembly: sources.assembly,
                output_directory,
                threads: sources.threads,
                reads: sources.reads,
            })?
        };
        let n_contigs = coverage.table.nrows();

        let filtered = {
            let _timer = crate::timing::scope("length_filter");
            coverage.filter_by_length(sources.lengths.min_contig)
        };
        if coverage.table.nrows() == 0 {
            bail!(
                "none of the {n_contigs} contigs in {} reach --min-contig-size {}, so there is \
                 nothing to bin",
                sources.assembly,
                sources.lengths.min_contig
            );
        }

        let counted = coverage
            .contig_lengths
            .iter()
            .zip(&coverage.contig_names)
            .filter(|(length, _)| **length >= sources.lengths.composition)
            .collect::<Vec<_>>();
        let scanned = {
            let _timer = crate::timing::scope("kmers");
            let held = match &sources.common.kmer_frequency_file {
                Some(path) => {
                    debug!("Reading TNF table.");
                    Some(KmerFrequencyTable::read(path)?)
                }
                None => kept_table(output_directory, counted.len(), sources.distance.kmer_size)?,
            };
            let kmer_size = held
                .as_ref()
                .map_or(sources.distance.kmer_size, KmerFrequencyTable::kmer_size);
            let floors = Floors {
                composition: held.is_none().then_some(sources.lengths.composition),
                sketch: sources.lengths.sketch,
                halves: sources.lengths.halves,
            };
            let mut scanned = scan(sources.assembly, kmer_size, floors)?;
            match held {
                Some(held) => scanned.composition = held,
                None if sources.distance.write_kmer_table => scanned
                    .composition
                    .write(table_path(output_directory, kmer_size))?,
                None => {}
            }
            scanned
        };
        let mut tnf = scanned.composition;

        // A coverage file need not list contigs in assembly order, so rows are matched by name.
        let row_of = tnf
            .contig_names
            .iter()
            .enumerate()
            .map(|(row, name)| (name.as_str(), row))
            .collect::<HashMap<_, _>>();
        let stray = counted
            .iter()
            .find(|(_, name)| !row_of.contains_key(name.as_str()));
        if stray.is_some() || row_of.len() < tnf.contig_names.len() {
            bail!(
                "the composition table does not hold the contigs of the coverage table, {}, so \
                 they were not built from {}",
                match stray {
                    Some((_, name)) => format!("starting with {name}"),
                    None => "and holds a contig twice".to_string(),
                },
                sources.assembly
            );
        }
        let order = counted
            .iter()
            .map(|(_, name)| row_of[name.as_str()])
            .collect::<Vec<_>>();
        if order.len() < tnf.contig_names.len()
            || order.iter().enumerate().any(|(at, row)| at != *row)
        {
            tnf.take_rows(&order);
        }
        let lengths = counted
            .iter()
            .map(|(length, _)| **length)
            .collect::<Vec<_>>();
        tnf.clr(&lengths)?;

        info!(
            "{} valid contigs, {} filtered contigs.",
            coverage.table.nrows(),
            filtered
        );
        Ok(Self {
            assembly: sources.assembly.to_string(),
            coverage,
            coverage_file,
            tnf,
            sketches: scanned.sketches,
            halves: scanned.halves,
            links: None,
        })
    }

    // Read after any reordering of the rows, since a link names its contigs by row.
    pub fn read_links(&mut self, graph: Option<&str>) -> Result<()> {
        self.links = graph
            .map(|path| crate::assembly_graph::read_links(path, &self.coverage.contig_names))
            .transpose()?;
        Ok(())
    }

    pub fn features(&self, settings: FeatureSettings) -> ContigFeatures<'_> {
        ContigFeatures::new(
            &self.coverage.table,
            &self.tnf.kmer_table,
            &self.coverage.contig_lengths,
        )
        .with_distance(settings.distance)
        .with_links(self.links.as_deref(), settings.link_weight)
        .with_sketches(self.sketches.as_ref())
    }
}
