use std::path::Path;

use anyhow::Result;
use log::info;

use super::bam_depth;
use super::coverage_calculator::{CoverageInputs, samples};
use super::depth_file::{self, TABLE};
use crate::cli::MakeArgs;

// The merged table is rebuilt from every sample in the directory, so a last run with no reads
// gathers samples that runs in parallel wrote.
pub fn run_make(args: &MakeArgs) -> Result<()> {
    let depth = Path::new(&args.output_directory);
    let inputs = CoverageInputs {
        assembly: args.assembly.as_deref().unwrap_or_default(),
        output_directory: &args.output_directory,
        depth,
        threads: args.runtime.threads,
        samples: &args.samples,
        mapping: &args.mapping,
    };
    let samples = samples(&inputs)?;
    if !samples.is_empty() {
        bam_depth::read(&samples, depth, inputs.threads)?;
    }
    let mut table = depth_file::kept(depth)?;
    table.write(depth.join(TABLE))?;
    info!(
        "Depth: {} samples by {} contigs in {}",
        table.sample_names.len(),
        table.contig_names.len(),
        args.output_directory
    );
    Ok(())
}
