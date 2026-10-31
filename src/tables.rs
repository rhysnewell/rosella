use std::collections::{HashMap, HashSet};
use std::ops::Range;
use std::sync::{Mutex, PoisonError};
use std::thread::JoinHandle;

use anyhow::{Result, bail};
use log::{debug, info};
use ndarray::Array2;

use crate::assembly_graph::Link;
use crate::cli::binning::DistanceParams;
use crate::cli::coverage::ReadArgs;
use crate::cli::runtime::Common;
use crate::coverage::coverage_calculator::{Coverage, CoverageInputs, calculate_coverage};
use crate::coverage::coverage_table::CoverageTable;
use crate::embedding::features::{ContigFeatures, FeatureSettings};
use crate::kmers::halves::Halves;
use crate::kmers::kmer_counting::{KmerFrequencyTable, kept_table, prefixes, table_path};
use crate::kmers::scan::{Floors, scan};
use crate::kmers::sketch::ContigSketches;
use crate::learned::Space;

pub struct Sources<'a> {
    pub assembly: &'a str,
    pub common: &'a Common,
    pub reads: &'a ReadArgs,
    pub distance: &'a DistanceParams,
    pub lengths: Lengths,
    pub threads: usize,
    pub only: Option<&'a HashSet<String>>,
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
    pub circles: HashSet<String>,
    pub learned: Option<Space>,
    pub windows: Option<std::path::PathBuf>,
    pub content: Option<u64>,
    ahead: Mutex<Option<(Range<usize>, JoinHandle<Rows>)>>,
}

type Rows = Result<(Array2<f64>, Array2<f64>)>;

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

        let Coverage {
            table: mut coverage,
            file: coverage_file,
            windows,
        } = {
            let _timer = crate::timing::scope("coverage");
            let source = &sources.reads.coverage;
            calculate_coverage(
                &CoverageInputs {
                    assembly: sources.assembly,
                    output_directory,
                    depth: &source.depth_directory.as_ref().map_or_else(
                        || std::path::Path::new(output_directory).join("depth"),
                        std::path::PathBuf::from,
                    ),
                    threads: sources.threads,
                    samples: &source.samples,
                    mapping: &sources.reads.mapping,
                },
                source.coverage_file.as_ref(),
                source.depth_directory.is_some(),
            )?
        };
        let n_contigs = coverage.table.nrows();

        let filtered = {
            let _timer = crate::timing::scope("length_filter");
            coverage.filter_by_length(sources.lengths.min_contig)
        };
        if let Some(only) = sources.only {
            coverage.retain_named(|name| only.contains(name));
        }
        if coverage.table.nrows() == 0 {
            bail!(
                "none of the {n_contigs} contigs in {} reach --min-contig-size {}",
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
                    debug!("k-mer table: {path}");
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
                "the composition table was not built from {}. {}",
                sources.assembly,
                match stray {
                    Some((_, name)) => format!("It lacks {name}"),
                    None => "It holds a contig twice".to_string(),
                }
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
        let lengths = drop_blank(&mut coverage, &mut tnf, lengths);
        tnf.clr(&lengths)?;

        match filtered {
            0 => info!("Contigs: {}", coverage.table.nrows()),
            _ => info!(
                "Contigs: {} ({filtered} under {} bp left out)",
                coverage.table.nrows(),
                sources.lengths.min_contig
            ),
        }
        Ok(Self {
            assembly: sources.assembly.to_string(),
            coverage,
            coverage_file,
            tnf,
            sketches: scanned.sketches,
            halves: scanned.halves,
            links: None,
            circles: HashSet::new(),
            learned: None,
            windows,
            content: scanned.content,
            ahead: Mutex::new(None),
        })
    }

    // Read after any reordering of the rows, since a link names its contigs by row.
    pub fn read_links(&mut self, graph: Option<&str>) -> Result<()> {
        let Some(path) = graph else {
            return Ok(());
        };
        let graph = crate::assembly_graph::read_graph(path, &self.coverage.contig_names)?;
        self.links = Some(graph.links);
        self.circles = graph.circles;
        Ok(())
    }

    pub fn read_rows(&self, rows: Range<usize>) -> Rows {
        let waiting = self
            .ahead
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .take();
        if let Some((ahead, reader)) = waiting
            && ahead == rows
        {
            return reader
                .join()
                .unwrap_or_else(|panic| std::panic::resume_unwind(panic));
        }
        read_rows(
            &self.assembly,
            &self.coverage_file,
            &self.coverage.contig_names[rows.clone()],
            &self.coverage.contig_lengths[rows],
            self.tnf.kmer_size(),
        )
    }

    // Inflating the assembly is one core's work, so the next band is read while this one is
    // measured. A read the caller never asks for only spends that core.
    pub fn read_ahead(&self, rows: Range<usize>) {
        let (assembly, file) = (self.assembly.clone(), self.coverage_file.clone());
        let names = self.coverage.contig_names[rows.clone()].to_vec();
        let lengths = self.coverage.contig_lengths[rows.clone()].to_vec();
        let kmer_size = self.tnf.kmer_size();
        let reader =
            std::thread::spawn(move || read_rows(&assembly, &file, &names, &lengths, kmer_size));
        *self.ahead.lock().unwrap_or_else(PoisonError::into_inner) = Some((rows, reader));
    }

    pub fn features(&self, settings: &FeatureSettings) -> ContigFeatures<'_> {
        ContigFeatures::new(
            &self.coverage.table,
            &self.tnf.kmer_table,
            &self.coverage.contig_lengths,
        )
        .with_distance(settings.distance.clone())
        .with_links(self.links.as_deref(), settings.link_weight)
        .with_sketches(self.sketches.as_ref())
    }
}

// An assembler writes a gap it could not fill as a contig of Ns. With no clean k-mer its zero row
// sits equally far from everything, so it is left out of the bins rather than offered to them.
fn drop_blank(
    coverage: &mut CoverageTable,
    tnf: &mut KmerFrequencyTable,
    lengths: Vec<usize>,
) -> Vec<usize> {
    let blank = tnf.blank_rows();
    if blank.is_empty() {
        return lengths;
    }
    debug!(
        "Blank contigs: {} (no run of {} clean bases)",
        blank.len(),
        tnf.kmer_size()
    );
    let names = blank
        .iter()
        .map(|row| tnf.contig_names[*row].clone())
        .collect::<HashSet<_>>();
    let kept = (0..lengths.len())
        .filter(|row| blank.binary_search(row).is_err())
        .collect::<Vec<_>>();
    tnf.take_rows(&kept);
    coverage.retain_named(|name| !names.contains(name));
    kept.iter().map(|row| lengths[*row]).collect()
}

fn read_rows(
    assembly: &str,
    coverage_file: &str,
    names: &[String],
    lengths: &[usize],
    kmer_size: usize,
) -> Rows {
    let names = names.iter().map(String::as_str).collect::<Vec<_>>();
    let coverage = {
        let _timer = crate::timing::scope("coverage");
        CoverageTable::rows_named(coverage_file, &names)?
    };
    let composition = {
        let _timer = crate::timing::scope("kmers");
        let pieces = names
            .iter()
            .copied()
            .zip(lengths.iter().copied())
            .collect::<Vec<_>>();
        prefixes(assembly, &pieces, kmer_size)?
    };
    Ok((composition, coverage))
}
