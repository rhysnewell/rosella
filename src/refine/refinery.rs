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
    coverage::coverage_table::CoverageTable,
    embedding::features::ContigFeatures,
    kmers::kmer_counting::KmerFrequencyTable,
    recover::recover_engine::UNBINNED,
    refine::{
        quality_table::read_contamination,
        splitter::{RefineSettings, Refiner},
    },
};

pub fn run_refine(args: &RefineArgs) -> Result<()> {
    RefineEngine::new(args)?.run()
}

struct RefineEngine {
    output_directory: String,
    assembly: String,
    coverage_table: CoverageTable,
    tnf_table: KmerFrequencyTable,
    genomes: Vec<std::path::PathBuf>,
    bin_quality: Option<String>,
    bin_tag: String,
    settings: RefineSettings,
    distance: crate::embedding::metrics::DistanceSettings,
    links: Option<Vec<crate::assembly_graph::Link>>,
    link_weight: f32,
}

impl RefineEngine {
    fn new(args: &RefineArgs) -> Result<Self> {
        let output_directory = args.common.output_directory.clone();
        let tables = crate::tables::Tables::build(&crate::tables::Sources {
            assembly: &args.assembly,
            common: &args.common,
            min_contig_size: args.binning.cutoff(),
            composition_from: args.binning.cutoff(),
            sketch_from: None,
            halves_from: None,
            coverage: &args.coverage,
            mapping: &args.mapping,
            filtering: &args.filtering,
            alignment: &args.alignment,
            trimming: &args.trimming,
            distance: &args.distance,
            threads: args.runtime.threads,
        })?;
        let (coverage_table, tnf_table, distance) = (tables.coverage, tables.tnf, tables.distance);

        let partition = args.binning.partition;

        let links = args
            .graph
            .assembly_graph
            .as_ref()
            .map(|path| crate::assembly_graph::read_links(path, &coverage_table.contig_names))
            .transpose()?;

        let genomes = args.genomes.discover()?;
        Ok(Self {
            assembly: args.assembly.clone(),
            output_directory,
            coverage_table,
            tnf_table,
            genomes,
            bin_quality: args.bin_quality.clone(),
            bin_tag: args.bin_tag.clone(),
            distance,
            links,
            link_weight: args.graph.assembly_graph_weight as f32,
            settings: RefineSettings {
                min_bin_size: args.binning.min_bin_size,
                max_bin_size: args.binning.max_bin_size,
                n_neighbours: args.graph.n_neighbours,
                knn_candidates: args.graph.candidates(),
                max_retries: args.refine.max_retries,
                seeds: crate::recover::settings::seeds(&args.seeds),
                max_contamination: Some(args.split_contamination),
                partition,
            },
        })
    }

    fn run(self) -> Result<()> {
        let indices = self
            .coverage_table
            .contig_names
            .iter()
            .enumerate()
            .map(|(index, name)| (name.as_str(), index))
            .collect::<HashMap<_, _>>();

        let mut bins = BTreeMap::new();
        let mut unchanged = Vec::new();
        let mut names = HashMap::new();
        let mut too_short = HashSet::new();
        for (position, genome) in self.genomes.iter().enumerate() {
            let (contigs, skipped) = self.contigs_in(genome, &indices)?;
            too_short.extend(skipped);
            if contigs.len() < crate::refine::bar::MIN_SPLIT_CONTIGS {
                debug!("{} has too few contigs to refine", genome.display());
                unchanged.push(contigs);
                continue;
            }
            names.insert(position, crate::bins::stem(genome));
            bins.insert(position, contigs);
        }

        if bins.is_empty() {
            info!("No genomes large enough to refine.");
        }

        let contamination = self.contamination(&names)?;
        let features = ContigFeatures::new(
            &self.coverage_table.table,
            &self.tnf_table.kmer_table,
            &self.coverage_table.contig_lengths,
        )
        .with_distance(self.distance)
        .with_links(self.links.as_deref(), self.link_weight);
        let mut refiner = Refiner::new(features, self.settings, bins, Vec::new())
            .with_contamination(contamination);
        refiner.run();

        let mut labelled = refiner.bins.into_values().collect::<Vec<_>>();
        labelled.extend(unchanged);
        self.write(&labelled, &refiner.unbinned, &too_short)?;

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
                    self.coverage_table.contig_names[*index].as_str(),
                    Some(label),
                );
            }
        }
        for index in unbinned {
            labels.insert(self.coverage_table.contig_names[*index].as_str(), None);
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
        let mut reader = parse_fastx_file(path::Path::new(&self.assembly))?;
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
