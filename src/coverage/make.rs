use std::path::Path;

use anyhow::{Result, bail};
use log::info;

use super::bam_depth;
use super::coverage_calculator::{CoverageInputs, samples};
use crate::cli::MakeArgs;

pub fn run_make(args: &MakeArgs) -> Result<()> {
    let source = &args.reads.coverage;
    if source.coverage_file.is_some() || source.depth_directory.is_some() {
        bail!(
            "make writes depth from reads or BAMs into --output-directory, so -C and --depth-directory are not read"
        );
    }
    let depth = Path::new(&args.output_directory);
    let inputs = CoverageInputs {
        assembly: args.assembly.as_deref().unwrap_or_default(),
        output_directory: &args.output_directory,
        depth,
        threads: args.runtime.threads,
        reads: &args.reads,
    };
    let samples = samples(&inputs)?;
    if samples.is_empty() {
        bail!("make needs reads or BAMs");
    }
    let table = bam_depth::read(&samples, depth, inputs.threads)?;
    info!(
        "Wrote the depth of {} samples over {} contigs to {}.",
        table.sample_names.len(),
        table.contig_names.len(),
        args.output_directory
    );
    Ok(())
}
