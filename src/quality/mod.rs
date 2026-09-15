pub mod bins;
pub mod orfs;

use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;


#[derive(Debug, Clone, Copy, Default)]
pub struct Quality {
    pub completeness: f64,
    pub contamination: f64,
}

impl Quality {
    pub fn score(&self, weight: f64) -> f64 {
        self.completeness - weight * self.contamination
    }
}

pub trait Scorer: Sync {
    fn score(&self, contigs: &[usize]) -> Quality;

    /// A partner that brings no feature the bin lacks cannot raise its completeness, which is
    /// most pairs, so this keeps the join off the scorer.
    fn features(&self, contigs: &[usize]) -> std::collections::HashSet<u32>;

    fn completeness_bar(&self, requested: f64) -> f64 {
        requested
    }
}

pub fn write_report<'a>(
    scorer: &dyn Scorer,
    bins: impl IntoIterator<Item = (String, &'a [usize])>,
    lengths: &[usize],
    path: &Path,
) -> Result<()> {
    let mut sink = BufWriter::new(std::fs::File::create(path)?);
    writeln!(sink, "bin\tcontigs\tbp\tcompleteness\tcontamination")?;
    for (bin, contigs) in bins {
        let held = scorer.score(contigs);
        let bp = contigs.iter().map(|contig| lengths[*contig]).sum::<usize>();
        writeln!(
            sink,
            "{bin}\t{}\t{bp}\t{:.2}\t{:.2}",
            contigs.len(),
            held.completeness,
            held.contamination
        )?;
    }
    sink.flush()?;
    Ok(())
}
