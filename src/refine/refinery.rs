use std::collections::HashMap;

use anyhow::Result;

use crate::cli::RefineArgs;
use crate::recover::origin::{Given, Origin};
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::quality_table::read_contamination;

pub fn run_refine(args: &RefineArgs) -> Result<()> {
    let contamination = match &args.bin_quality {
        Some(path) => read_contamination(path)?,
        None => HashMap::new(),
    };
    let given = Given::read(
        &args.genomes.discover()?,
        format!("rosella_{}_", args.bin_tag),
        contamination,
        args.bin_quality.as_ref().map(|_| args.split_contamination),
    )?;
    RecoverEngine::new(&args.run, Origin::Given(given))?.run()
}
