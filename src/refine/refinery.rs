use std::{
    collections::{BTreeMap, HashMap, hash_map::Entry},
    fs::{File, OpenOptions},
    io::BufWriter,
    path,
};

use anyhow::Result;
use log::{debug, info};
use needletail::{
    parse_fastx_file,
    parser::{LineEnding, write_fasta},
};

use crate::{
    cli::RefineArgs,
    clustering::objective::ObjectiveChoice,
    coverage::{
        coverage_calculator::{CoverageInputs, calculate_coverage},
        coverage_table::CoverageTable,
    },
    embedding::features::ContigFeatures,
    kmers::kmer_counting::{KmerFrequencyTable, count_kmers},
    recover::recover_engine::{RECOVER_FASTA_EXTENSION, UNBINNED},
    refine::{
        checkm::read_checkm,
        splitter::{RefineSettings, Refiner},
    },
};

pub fn run_refine(args: RefineArgs) -> Result<()> {
    RefineEngine::new(&args)?.run()
}

/// Replaces bird_tool_utils' version, which read a `genome-fasta-list` argument rosella
/// never defined. Sorted, because directory order is not stable and the bin names are.
fn genomes_to_refine(args: &RefineArgs) -> Result<Vec<String>> {
    if !args.genome_fasta_files.is_empty() {
        return Ok(args.genome_fasta_files.clone());
    }

    let Some(directory) = &args.genome_fasta_directory else {
        bail!("Pass the bins to refine with --genome-fasta-files or --genome-fasta-directory");
    };
    let wanted = args.genome_fasta_extension.trim_start_matches('.');
    let mut genomes = std::fs::read_dir(directory)?
        .filter_map(|entry| {
            let path = entry.ok()?.path();
            let found = path.extension()?.to_str()?;
            (found == wanted).then(|| path.to_string_lossy().into_owned())
        })
        .collect::<Vec<_>>();
    if genomes.is_empty() {
        bail!("{} holds no .{} files to refine", directory, wanted);
    }
    genomes.sort_unstable();
    Ok(genomes)
}

struct RefineEngine {
    output_directory: String,
    assembly: String,
    coverage_table: CoverageTable,
    tnf_table: KmerFrequencyTable,
    genomes: Vec<String>,
    checkm_results: Option<String>,
    min_contig_count: usize,
    bin_tag: String,
    settings: RefineSettings,
    distance: crate::embedding::metrics::DistanceSettings,
    objective: ObjectiveChoice,
}

impl RefineEngine {
    fn new(args: &RefineArgs) -> Result<Self> {
        let output_directory = args.common.output_directory.clone();
        std::fs::create_dir_all(&output_directory)?;

        let min_contig_size = args.binning.min_contig_size;
        let mut coverage_table = calculate_coverage(&CoverageInputs {
            assembly: args.assembly.as_deref(),
            output_directory: &output_directory,
            threads: args.common.threads,
            coverage: &args.coverage,
            mapping: &args.mapping,
            filtering: &args.filtering,
            alignment: &args.alignment,
            trimming: &args.trimming,
        })?;
        let n_contigs = coverage_table.table.nrows();
        let filtered_contigs = coverage_table.filter_by_length(min_contig_size)?;

        let mut tnf_table = if let Some(path) = &args.common.kmer_frequency_file {
            info!("Reading TNF table.");
            KmerFrequencyTable::read(path)?
        } else {
            info!("Calculating TNF table.");
            let assembly = args
                .assembly
                .as_deref()
                .ok_or_else(|| anyhow!("Counting tetranucleotides needs --assembly"))?;
            count_kmers(assembly, &output_directory, Some(n_contigs))?
        };
        tnf_table.filter_by_name(&filtered_contigs)?;
        if coverage_table.table.nrows() != tnf_table.kmer_table.nrows() {
            bail!(
                "Coverage table has {} contigs and the TNF table {}. Refinement indexes both \
                 by the same contig, so the mismatch surfaces as a panic inside the splitter.",
                coverage_table.table.nrows(),
                tnf_table.kmer_table.nrows()
            );
        }

        let genomes = genomes_to_refine(args)?;

        Ok(Self {
            assembly: args
                .assembly
                .clone()
                .ok_or_else(|| anyhow!("Writing the refined bins needs --assembly"))?,
            output_directory,
            coverage_table,
            tnf_table,
            genomes,
            checkm_results: args.checkm_results.clone(),
            min_contig_count: args.min_contig_count,
            bin_tag: args.bin_tag.clone(),
            distance: crate::recover::recover_engine::distance_settings(&args.distance),
            settings: RefineSettings {
                min_bin_size: args.binning.min_bin_size,
                max_bin_size: args.binning.max_bin_size,
                n_neighbours: args.binning.n_neighbours,
                max_retries: args.binning.max_retries,
                seeds: crate::recover::recover_engine::seeds(args.common.seed, &args.seeds),
                overrides: crate::recover::recover_engine::embed_overrides(&args.overrides),
                max_contamination: Some(args.max_contamination),
                largest_cluster: args.binning.max_cluster_size,
            },
            objective: ObjectiveChoice::parse(&args.binning.objective)
                .ok_or_else(|| anyhow!("unknown objective {}", args.binning.objective))?,
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
        for (position, genome) in self.genomes.iter().enumerate() {
            let contigs = self.contigs_in(genome, &indices)?;
            if contigs.len() < self.min_contig_count {
                debug!("{} has too few contigs to refine", genome);
                unchanged.push(contigs);
                continue;
            }
            names.insert(position, stem(genome));
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
        .with_distance(self.distance);
        let scorer = self.objective.build(
            &self.coverage_table.contig_lengths,
            self.settings.min_bin_size,
        );
        let mut refiner = Refiner::new(features, None, &scorer, self.settings, bins, Vec::new())
            .with_contamination(contamination);
        refiner.run();

        let mut labelled = refiner.bins.into_values().collect::<Vec<_>>();
        labelled.extend(unchanged);
        self.write(labelled, refiner.unbinned)?;

        crate::timing::report(
            path::Path::new(&self.output_directory).join(crate::timing::TIMINGS_FILE),
        )
    }

    /// Contigs of a genome as indices into the coverage table. Anything the length filter
    /// dropped is not there to be refined, so it is left out.
    fn contigs_in(&self, genome: &str, indices: &HashMap<&str, usize>) -> Result<Vec<usize>> {
        let mut reader = parse_fastx_file(path::Path::new(genome))?;
        let mut contigs = Vec::new();
        while let Some(record) = reader.next() {
            let seqrec = record?;
            let name = std::str::from_utf8(seqrec.id())?;
            if let Some(index) = indices.get(name) {
                contigs.push(*index);
            }
        }
        contigs.sort_unstable();
        Ok(contigs)
    }

    fn contamination(&self, names: &HashMap<usize, String>) -> Result<HashMap<usize, f64>> {
        let Some(path) = &self.checkm_results else {
            return Ok(HashMap::new());
        };

        let stats = read_checkm(path)?;
        let mut contamination = HashMap::new();
        for (bin_id, name) in names.iter() {
            match stats.get(name) {
                Some((_, contaminated)) => {
                    contamination.insert(*bin_id, *contaminated);
                }
                None => debug!("{} has no row in {}", name, path),
            }
        }
        info!("Read quality estimates for {} bins.", contamination.len());
        Ok(contamination)
    }

    /// Written by contig name off the assembly, so the bin files and the coverage table
    /// never have to agree on an ordering.
    fn write(&self, bins: Vec<Vec<usize>>, unbinned: Vec<usize>) -> Result<()> {
        let mut labels = HashMap::new();
        for (label, contigs) in bins.iter().enumerate() {
            for index in contigs.iter() {
                labels.insert(
                    self.coverage_table.contig_names[*index].as_str(),
                    format!("{label}"),
                );
            }
        }
        for index in unbinned.iter() {
            labels.insert(
                self.coverage_table.contig_names[*index].as_str(),
                UNBINNED.to_string(),
            );
        }

        let mut reader = parse_fastx_file(path::Path::new(&self.assembly))?;
        let mut writers: HashMap<String, BufWriter<File>> = HashMap::new();
        let mut written = 0;
        while let Some(record) = reader.next() {
            let seqrec = record?;
            let name = std::str::from_utf8(seqrec.id())?;
            let Some(label) = labels.get(name) else {
                continue;
            };

            let writer = match writers.entry(label.clone()) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => {
                    let bin_path = path::Path::new(&self.output_directory).join(format!(
                        "rosella_{}_{}{}",
                        self.bin_tag,
                        entry.key(),
                        RECOVER_FASTA_EXTENSION
                    ));
                    let file = OpenOptions::new()
                        .append(true)
                        .create(true)
                        .open(bin_path)?;
                    entry.insert(BufWriter::new(file))
                }
            };
            write_fasta(seqrec.id(), &seqrec.seq(), writer, LineEnding::Unix)?;
            written += 1;
        }

        info!(
            "Wrote {} contigs into {} bins from {} input genomes.",
            written,
            bins.len(),
            self.genomes.len()
        );
        Ok(())
    }
}

fn stem(genome: &str) -> String {
    path::Path::new(genome)
        .file_stem()
        .map(|stem| stem.to_string_lossy().into_owned())
        .unwrap_or_else(|| genome.to_string())
}
