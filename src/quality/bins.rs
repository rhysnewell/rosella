use std::collections::{BTreeMap, HashMap};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use log::info;
use needletail::parse_fastx_file;

use crate::cli::ScoreArgs;
use crate::markers::{Annotator, MarkerRules};

struct Layout {
    names: Vec<String>,
    lengths: Vec<usize>,
    bins: BTreeMap<String, Vec<usize>>,
}

// A contig named in two bins is called once and counted in both.
fn read_bins(paths: &[PathBuf], contigs: &Path) -> Result<Layout> {
    let mut held = Layout {
        names: Vec::new(),
        lengths: Vec::new(),
        bins: BTreeMap::new(),
    };
    let mut index = HashMap::new();
    let mut sink = BufWriter::new(std::fs::File::create(contigs)?);
    for path in paths {
        let mut members = Vec::new();
        let mut reader = parse_fastx_file(path)?;
        while let Some(record) = reader.next() {
            let record = record?;
            let name = crate::contig_id(record.id())?.to_string();
            let at = *index.entry(name.clone()).or_insert(held.names.len());
            if at == held.names.len() {
                let sequence = record.seq();
                writeln!(sink, ">{name}")?;
                sink.write_all(&sequence)?;
                writeln!(sink)?;
                held.names.push(name);
                held.lengths.push(sequence.len());
            }
            members.push(at);
        }
        if !members.is_empty() {
            held.bins.insert(crate::bins::stem(path), members);
        }
    }
    sink.flush()?;
    if held.names.is_empty() {
        bail!("the bins hold no contigs");
    }
    Ok(held)
}

pub fn run_score(args: &ScoreArgs) -> Result<()> {
    let paths = crate::bins::discover(
        &args.genome_fasta_files,
        args.genome_fasta_directory.as_ref(),
        &args.genome_fasta_extension,
    )?;
    let directory = tempfile::tempdir()?;
    let contigs = directory.path().join("binned.fna");
    let held = read_bins(&paths, &contigs)?;
    info!(
        "Scoring {} bins over {} contigs.",
        held.bins.len(),
        held.names.len()
    );

    let annotator = Annotator {
        assembly: contigs.to_string_lossy().into_owned(),
        threads: args.runtime.threads,
        shards: None,
        rules: MarkerRules::default(),
        cache: None,
        checkm: true,
    };
    let mut scorer = annotator
        .annotate(0..usize::MAX)?
        .select(&held.names)?
        .with_lengths(held.lengths.clone());
    let every = (0..held.names.len()).collect::<Vec<_>>();
    annotator.complete_checkm(&mut scorer, &every, &held.names)?;

    if let Some(path) = &args.marker_report {
        scorer.report(&held.names, Path::new(path))?;
        info!("Wrote every marker hit to {path}.");
    }
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
