use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use log::{info, warn};
use needletail::parse_fastx_file;

use crate::cli::ScoreArgs;
use crate::markers::{MarkerAnnotation, MarkerRules};


struct Layout {
    names: Vec<String>,
    lengths: Vec<usize>,
    bins: BTreeMap<String, Vec<usize>>,
}

fn read_bins(paths: &[PathBuf], min_contig_size: usize) -> Result<Layout> {
    let mut held = Layout {
        names: Vec::new(),
        lengths: Vec::new(),
        bins: BTreeMap::new(),
    };
    let mut short = 0;
    for path in paths {
        let mut contigs = Vec::new();
        let mut reader = parse_fastx_file(path)?;
        while let Some(record) = reader.next() {
            let record = record?;
            let length = record.seq().len();
            if length < min_contig_size {
                short += 1;
                continue;
            }
            contigs.push(held.names.len());
            held.names.push(crate::contig_id(record.id())?.to_string());
            held.lengths.push(length);
        }
        if !contigs.is_empty() {
            held.bins.insert(crate::bins::stem(path), contigs);
        }
    }
    if short > 0 {
        info!("{short} contigs under the minimum size were left out of the scored bins.");
    }
    if held.names.is_empty() {
        bail!("every contig in every bin is under {min_contig_size} bp");
    }
    Ok(held)
}

pub fn run_score(args: ScoreArgs) -> Result<()> {
    let paths = crate::bins::discover(
        &args.genome_fasta_files,
        args.genome_fasta_directory.as_ref(),
        &args.genome_fasta_extension,
    )?;
    let held = read_bins(&paths, args.min_contig_size)?;
    info!(
        "Scoring {} bins over {} contigs.",
        held.bins.len(),
        held.names.len()
    );

    let annotation = MarkerAnnotation::build(
        &args.assembly,
        args.min_contig_size,
        crate::quality::orfs::GeneRules {
            min_length: args.markers.gene_min_length,
            model_depth: args.markers.gene_model_depth,
        },
        args.runtime.threads,
        args.markers.hmm_shards.map(usize::from),
        MarkerRules {
            fragment_span: args.markers.marker_fragment_span,
            bar_offset: args.markers.marker_bar_offset,
        },
        args.markers.marker_cache.as_deref().map(Path::new),
    )?;

    let annotated = annotation
        .names()
        .iter()
        .cloned()
        .collect::<std::collections::HashSet<_>>();
    let absent = held
        .names
        .iter()
        .filter(|name| !annotated.contains(*name))
        .count();
    if absent > 0 {
        warn!("{absent} binned contigs are not in the assembly, so they score as featureless.");
    }
    let scorer = annotation.select_present(&held.names)?;

    crate::quality::write_report(
        &scorer,
        held.bins
            .iter()
            .map(|(name, contigs)| (name.clone(), contigs.as_slice())),
        &held.lengths,
        Path::new(&args.output_file),
    )?;
    info!("Wrote the quality table to {}.", args.output_file);

    let output = Path::new(&args.output_file);
    crate::timing::report(
        output
            .parent()
            .unwrap_or(Path::new("."))
            .join(crate::timing::TIMINGS_FILE),
    )?;
    Ok(())
}
