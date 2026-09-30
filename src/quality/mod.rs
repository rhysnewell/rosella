pub mod bins;
pub mod orfs;

use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

#[derive(Debug, Clone, Copy, Default)]
pub struct Quality {
    pub completeness: f64,
    pub contamination: f64,
    pub set: u16,
}

impl Quality {
    pub fn score(&self, weight: f64) -> f64 {
        self.completeness - weight * self.contamination
    }

    pub fn clears(&self, bars: Bars) -> bool {
        self.completeness >= bars.completeness && self.contamination <= bars.contamination
    }
}

// The standard error of a squared worth sum under resampling the marker catalogue, by the
// infinitesimal jackknife (Jaeckel 1972), so it needs no replicates and no seed.
pub fn edge_spread(sides: impl IntoIterator<Item = (f64, Vec<(usize, f64)>)>) -> f64 {
    let mut gradient = std::collections::BTreeMap::<usize, f64>::new();
    for (sign, points) in sides {
        let count = points.len() as f64;
        let worth = points.iter().map(|(_, point)| point).sum::<f64>() / count;
        if points.is_empty() || worth <= 0.0 {
            continue;
        }
        for (marker, point) in points {
            *gradient.entry(marker).or_default() += sign * 2.0 * worth * (point - worth) / count;
        }
    }
    gradient
        .values()
        .map(|slope| slope * slope)
        .sum::<f64>()
        .sqrt()
}

#[derive(Debug, Clone, Copy)]
pub struct Bars {
    pub completeness: f64,
    pub contamination: f64,
}

pub trait Scorer: Sync {
    fn score(&self, contigs: &[usize]) -> Quality;

    /// A partner that brings no feature the bin lacks cannot raise its completeness, which is
    /// most pairs, so this keeps the join off the scorer.
    fn features(&self, contigs: &[usize]) -> std::collections::HashSet<u32>;

    fn set_name(&self, _set: u16) -> &str {
        ""
    }

    // A scorer without a marker catalogue has nothing to resample, so its edges carry no spread.
    fn points(&self, _contigs: &[usize], _weight: f64) -> Vec<(usize, f64)> {
        Vec::new()
    }
}

pub fn write_report<'a>(
    scorer: &dyn Scorer,
    bins: impl IntoIterator<Item = (String, &'a [usize])>,
    lengths: &[usize],
    path: &Path,
) -> Result<()> {
    let mut sink = BufWriter::new(std::fs::File::create(path)?);
    writeln!(sink, "bin\tcontigs\tbp\tcompleteness\tcontamination\tset")?;
    for (bin, contigs) in bins {
        let held = scorer.score(contigs);
        let bp = contigs.iter().map(|contig| lengths[*contig]).sum::<usize>();
        writeln!(
            sink,
            "{bin}\t{}\t{bp}\t{:.2}\t{:.2}\t{}",
            contigs.len(),
            held.completeness,
            held.contamination,
            scorer.set_name(held.set)
        )?;
    }
    sink.flush()?;
    Ok(())
}
