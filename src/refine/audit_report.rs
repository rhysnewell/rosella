use std::collections::BTreeMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::refine::report_context::{Context, Inputs};

pub fn write(path: &Path, bins: &BTreeMap<usize, Vec<usize>>, inputs: &Inputs<'_>) -> Result<()> {
    let context = Context::of(bins, inputs.features, inputs.quality);
    let mut out = BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        out,
        "contig\tlength\town_bin\town_share\trival_bin\trival_share\tclaim"
    )?;
    for (contig, label) in context.owned() {
        let Some((own, rivals)) = context.rivals(contig, label, inputs) else {
            continue;
        };
        for rival in rivals {
            writeln!(
                out,
                "{}\t{}\t{}\t{:.6}\t{}\t{:.6}\t{:.6}",
                inputs.names[contig],
                inputs.lengths[contig],
                label,
                own,
                rival.bin,
                rival.share,
                rival.claim
            )?;
        }
    }
    out.flush()?;
    Ok(())
}
