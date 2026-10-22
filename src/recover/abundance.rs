use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use ndarray::s;

use crate::coverage::coverage_table::CoverageTable;

pub const ABUNDANCE_FILE: &str = "abundance.tsv";

pub fn write(bins: &[(String, Vec<usize>)], coverage: &CoverageTable, path: &Path) -> Result<()> {
    let depths = coverage.table.slice(s![.., ..;2]);
    let samples = depths.ncols();
    let mut totals = vec![0.0; samples];
    let held = bins
        .iter()
        .map(|(_, contigs)| {
            let mut bp = 0usize;
            let mut weighted = vec![0.0; samples];
            for contig in contigs {
                let length = coverage.contig_lengths[*contig];
                bp += length;
                if *contig >= depths.nrows() {
                    continue;
                }
                for (sample, depth) in depths.row(*contig).iter().enumerate() {
                    weighted[sample] += depth * length as f64;
                }
            }
            for (total, sum) in totals.iter_mut().zip(&weighted) {
                *total += sum;
            }
            (bp, weighted)
        })
        .collect::<Vec<_>>();

    let mut sink = BufWriter::new(std::fs::File::create(path)?);
    write!(sink, "bin\tbp")?;
    for sample in 0..samples {
        let name = coverage
            .sample_names
            .get(sample)
            .cloned()
            .unwrap_or_else(|| format!("sample_{}", sample + 1));
        write!(sink, "\t{name}_depth\t{name}_share")?;
    }
    writeln!(sink)?;
    for ((name, _), (bp, weighted)) in bins.iter().zip(held) {
        write!(sink, "{name}\t{bp}")?;
        for (sum, total) in weighted.iter().zip(&totals) {
            let depth = if bp == 0 { 0.0 } else { sum / bp as f64 };
            let share = if *total > 0.0 {
                100.0 * sum / total
            } else {
                0.0
            };
            write!(sink, "\t{depth:.4}\t{share:.4}")?;
        }
        writeln!(sink)?;
    }
    sink.flush()?;
    Ok(())
}
