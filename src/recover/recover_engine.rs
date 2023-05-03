use std::process::Command;

use anyhow::Result;
use log::info;

use crate::{coverage::{coverage_calculator::calculate_coverage, coverage_table::CoverageTable}, kmers::kmer_counting::{count_kmers, KmerFrequencyTable}, sketch::contig_sketcher::{ContigSketchResult, sketch_contigs}};


pub fn run_recover(m: &clap::ArgMatches) -> Result<()> {
    let mut recover_engine = RecoverEngine::new(m)?;
    recover_engine.run()?;
    Ok(())
}

struct RecoverEngine {
    output_directory: String,
    coverage_table: CoverageTable,
    contig_sketches: ContigSketchResult,
}

impl RecoverEngine {
    pub fn new(m: &clap::ArgMatches) -> Result<Self> {
        // create output directory
        let output_directory = m.get_one::<String>("output-directory").unwrap().clone();
        // create the output directory but do not fail if it already exists
        std::fs::create_dir_all(&output_directory)?;
        info!("Calculating contig coverages.");
        let coverage_table = calculate_coverage(m)?;
        info!("Counting kmers.");
        let contig_sketches = sketch_contigs(m)?;

        Ok(
            Self {
                output_directory,
                coverage_table,
                contig_sketches,
            }
        )
    }

    /// Runs through the rosella bin recovery pipeline
    pub fn run(&mut self) -> Result<()> {
        
        Ok(())
    }

    fn run_flight(&self) -> Result<()> {
        let mut flight_command = Command::new("flight");
        // flight_command.arg("bin")
        Ok(())
    }
}