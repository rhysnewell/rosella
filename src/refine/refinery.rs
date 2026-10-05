use std::{
    collections::{BTreeMap, HashMap, HashSet},
    path,
};

use anyhow::Result;
use log::{debug, info, warn};
use needletail::parse_fastx_file;

use crate::{
    bin_files::BinFiles,
    cli::RefineArgs,
    embedding::features::FeatureSettings,
    recover::recover_engine::UNBINNED,
    refine::{
        quality_table::read_contamination,
        splitter::{RefineSettings, Refiner},
    },
    tables::{Lengths, Sources, Tables},
};

pub fn run_refine(args: &RefineArgs) -> Result<()> {
    RefineEngine::new(args)?.run()
}

struct RefineEngine {
    tables: Tables,
    genomes: Vec<path::PathBuf>,
    bin_quality: Option<String>,
    bin_tag: String,
    settings: RefineSettings,
    features: FeatureSettings,
    output_directory: String,
}

impl RefineEngine {
    fn new(args: &RefineArgs) -> Result<Self> {
        let mut tables = Tables::build(&Sources {
            assembly: &args.assembly,
            common: &args.common,
            reads: &args.reads,
            distance: &args.distance,
            lengths: Lengths::from(args.binning.cutoff(), args.binning.cutoff()),
            threads: args.runtime.threads,
        })?;
        tables.read_links(args.graph.assembly_graph.as_deref())?;
        Ok(Self {
            tables,
            genomes: args.genomes.discover()?,
            bin_quality: args.bin_quality.clone(),
            bin_tag: args.bin_tag.clone(),
            features: crate::recover::settings::feature_settings(&args.graph),
            output_directory: args.common.output_directory.clone(),
            settings: RefineSettings {
                min_bin_size: args.binning.min_bin_size,
                max_bin_size: args.binning.max_bin_size,
                knn: crate::recover::settings::knn_settings(&args.graph),
                max_retries: args.refine.max_retries,
                seeds: crate::recover::settings::seeds(&args.seeds),
                max_contamination: Some(args.split_contamination),
                partition: args.binning.partition,
            },
        })
    }

    fn run(self) -> Result<()> {
        let indices = self
            .tables
            .coverage
            .contig_names
            .iter()
            .enumerate()
            .map(|(index, name)| (name.as_str(), index))
            .collect::<HashMap<_, _>>();

        // Every genome goes to the refiner, as in recover, so a closed genome in one contig
        // still sets the genome floor the splits are held to.
        let mut bins = BTreeMap::new();
        let mut names = HashMap::new();
        let mut too_short = HashSet::new();
        for (position, genome) in self.genomes.iter().enumerate() {
            let (contigs, skipped) = self.contigs_in(genome, &indices)?;
            too_short.extend(skipped);
            names.insert(position, crate::bins::stem(genome));
            bins.insert(position, contigs);
        }

        let contamination = self.contamination(&names)?;
        let features = self.tables.features(&self.features);
        let refined = Refiner::new(features, self.settings, bins, Vec::new())
            .with_contamination(contamination)
            .run();

        let labelled = refined.bins.into_values().collect::<Vec<_>>();
        self.write(&labelled, &refined.unbinned, &too_short)?;

        crate::timing::report(
            path::Path::new(&self.output_directory).join(crate::timing::TIMINGS_FILE),
        )
    }

    // The second list keeps length-filtered contigs out of the refined bins without also
    // losing them from the output.
    fn contigs_in(
        &self,
        genome: &path::Path,
        indices: &HashMap<&str, usize>,
    ) -> Result<(Vec<usize>, Vec<String>)> {
        let mut reader = parse_fastx_file(genome)?;
        let mut contigs = Vec::new();
        let mut skipped = Vec::new();
        while let Some(record) = reader.next() {
            let seqrec = record?;
            let name = crate::contig_id(seqrec.id())?;
            match indices.get(name) {
                Some(index) => contigs.push(*index),
                None => skipped.push(name.to_string()),
            }
        }
        contigs.sort_unstable();
        Ok((contigs, skipped))
    }

    fn contamination(&self, names: &HashMap<usize, String>) -> Result<HashMap<usize, f64>> {
        let Some(path) = &self.bin_quality else {
            return Ok(HashMap::new());
        };

        let stats = read_contamination(path)?;
        let mut contamination = HashMap::new();
        for (bin_id, name) in names.iter() {
            match stats.get(name) {
                Some(contaminated) => {
                    contamination.insert(*bin_id, *contaminated);
                }
                None => debug!("{} has no row in {}", name, path),
            }
        }
        info!("Read quality estimates for {} bins.", contamination.len());
        Ok(contamination)
    }

    // Written by contig name off the assembly, so the bin files and the coverage table
    // never have to agree on an ordering.
    fn write(
        &self,
        bins: &[Vec<usize>],
        unbinned: &[usize],
        too_short: &HashSet<String>,
    ) -> Result<()> {
        let mut labels = HashMap::new();
        for (label, contigs) in bins.iter().enumerate() {
            for index in contigs {
                labels.insert(
                    self.tables.coverage.contig_names[*index].as_str(),
                    Some(label),
                );
            }
        }
        for index in unbinned {
            labels.insert(self.tables.coverage.contig_names[*index].as_str(), None);
        }

        let directory = path::Path::new(&self.output_directory);
        let mut files = BinFiles::new(|label: &Option<usize>| {
            directory.join(format!(
                "rosella_{}_{}.{}",
                self.bin_tag,
                label.map_or_else(|| UNBINNED.to_string(), |label| label.to_string()),
                crate::defaults::FASTA_EXTENSION
            ))
        });
        let mut reader = parse_fastx_file(path::Path::new(&self.tables.assembly))?;
        let mut written = 0;
        let mut short = 0;
        let mut skipped = 0;
        while let Some(record) = reader.next() {
            let seqrec = record?;
            let name = crate::contig_id(seqrec.id())?;
            let label = match labels.get(name) {
                Some(label) => *label,
                None if too_short.contains(name) => {
                    short += 1;
                    None
                }
                None => {
                    skipped += 1;
                    continue;
                }
            };
            files.write(&label, seqrec.id(), &seqrec.seq())?;
            written += 1;
        }
        files.finish()?;

        info!(
            "Wrote {} contigs into {} bins from {} input genomes.",
            written,
            bins.len(),
            self.genomes.len()
        );
        if short > 0 {
            warn!(
                "{} contigs of the input genomes are under --min-contig-size, so they could \
                 not be refined and were written to {}",
                short, UNBINNED
            );
        }
        if skipped > 0 {
            warn!(
                "{} assembly contigs belong to no input genome, so nothing was written for them",
                skipped
            );
        }
        Ok(())
    }
}
