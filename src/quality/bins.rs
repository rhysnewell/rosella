use std::collections::btree_map::Entry;
use std::collections::{BTreeMap, HashMap};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use log::info;
use needletail::parse_fastx_file;

use crate::cli::ScoreArgs;
use crate::markers::{Annotator, MarkerRules};
use crate::quality::bases::Bases;
use crate::quality::report::{Bin, Scored};

struct Layout {
    names: Vec<String>,
    lengths: Vec<usize>,
    bases: HashMap<usize, Bases>,
    bins: BTreeMap<String, Vec<usize>>,
}

// A contig named in two bins is called once and counted in both.
fn read_bins(paths: &[PathBuf], contigs: &Path) -> Result<Layout> {
    let mut held = Layout {
        names: Vec::new(),
        lengths: Vec::new(),
        bases: HashMap::new(),
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
                held.bases.insert(at, Bases::count(&sequence));
            }
            members.push(at);
        }
        if members.is_empty() {
            continue;
        }
        match held.bins.entry(crate::bins::stem(path)) {
            Entry::Occupied(taken) => {
                bail!(
                    "{} and another bin are both named {}",
                    path.display(),
                    taken.key()
                )
            }
            Entry::Vacant(slot) => {
                slot.insert(members);
            }
        }
    }
    sink.flush()?;
    if held.names.is_empty() {
        bail!("the bins hold no contigs");
    }
    Ok(held)
}

pub fn run_score(args: &ScoreArgs) -> Result<()> {
    let paths = args.genomes.discover()?;
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
    let members = held.bins.values().map(Vec::as_slice).collect::<Vec<_>>();
    let strain = annotator.strain_heterogeneity(&scorer, &members, &held.names)?;

    if let Some(path) = &args.marker_report {
        let placed = held
            .bins
            .iter()
            .flat_map(|(name, contigs)| contigs.iter().map(move |contig| (name.as_str(), *contig)));
        scorer.report(placed, &held.names, path)?;
        info!("Wrote every marker hit to {}.", path.display());
    }
    let bins = held
        .bins
        .iter()
        .zip(strain)
        .map(|((name, contigs), strain)| Bin {
            name: name.clone(),
            contigs,
            strain,
        })
        .collect::<Vec<_>>();
    Scored {
        markers: &scorer,
        names: &held.names,
        lengths: &held.lengths,
        bases: &held.bases,
    }
    .write(&bins, Path::new(&args.output_file))?;
    info!("Wrote the quality tables beside {}.", args.output_file);

    let output = Path::new(&args.output_file);
    crate::timing::report(
        output
            .parent()
            .unwrap_or(Path::new("."))
            .join(crate::timing::TIMINGS_FILE),
    )?;
    Ok(())
}
