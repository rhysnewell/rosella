use std::{
    collections::{BTreeMap, HashMap, hash_map::Entry},
    fs::{File, OpenOptions},
    io::BufWriter,
    path,
};

use anyhow::Result;
use bird_tool_utils::clap_utils::parse_list_of_genome_fasta_files;
use log::{debug, info};
use needletail::{
    parse_fastx_file,
    parser::{LineEnding, write_fasta},
};

use crate::{
    clustering::objective::Dbcv,
    coverage::{coverage_calculator::calculate_coverage, coverage_table::CoverageTable},
    embedding::features::ContigFeatures,
    kmers::kmer_counting::{KmerFrequencyTable, count_kmers},
    recover::recover_engine::{RECOVER_FASTA_EXTENSION, UNBINNED},
    refine::{
        checkm::read_checkm,
        splitter::{RefineSettings, Refiner},
    },
};

pub fn run_refine(m: &clap::ArgMatches) -> Result<()> {
    RefineEngine::new(m)?.run()
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
}

impl RefineEngine {
    fn new(m: &clap::ArgMatches) -> Result<Self> {
        let output_directory = m.get_one::<String>("output-directory").unwrap().clone();
        std::fs::create_dir_all(&output_directory)?;

        let assembly = m.get_one::<String>("assembly").unwrap().clone();
        let min_contig_size = *m.get_one::<usize>("min-contig-size").unwrap();
        let mut coverage_table = calculate_coverage(m)?;
        let n_contigs = coverage_table.table.nrows();
        let filtered_contigs = coverage_table.filter_by_length(min_contig_size)?;

        let mut tnf_table = if let Some(path) = m.get_one::<String>("kmer-frequency-file") {
            info!("Reading TNF table.");
            KmerFrequencyTable::read(path)?
        } else {
            info!("Calculating TNF table.");
            count_kmers(m, Some(n_contigs))?
        };
        tnf_table.filter_by_name(&filtered_contigs)?;

        let genomes = match parse_list_of_genome_fasta_files(m, true) {
            Ok(genomes) => genomes,
            Err(e) => bail!("Failed to parse the genomes to refine: {}", e),
        };

        Ok(Self {
            output_directory,
            assembly,
            coverage_table,
            tnf_table,
            genomes,
            checkm_results: m.get_one::<String>("checkm-results").cloned(),
            min_contig_count: *m.get_one::<usize>("min-contig-count").unwrap(),
            bin_tag: m.get_one::<String>("bin-tag").unwrap().clone(),
            settings: RefineSettings {
                min_bin_size: *m.get_one::<usize>("min-bin-size").unwrap(),
                max_bin_size: *m.get_one::<usize>("max-bin-size").unwrap(),
                n_neighbours: *m.get_one::<usize>("n-neighbours").unwrap(),
                max_retries: *m.get_one::<usize>("max-retries").unwrap(),
                seed: *m.get_one::<u64>("seed").unwrap(),
                overrides: crate::recover::recover_engine::embed_overrides(m),
                max_contamination: m.get_one::<f64>("max-contamination").copied(),
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
        );
        let mut refiner = Refiner::new(features, None, &Dbcv, self.settings, bins, Vec::new())
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
