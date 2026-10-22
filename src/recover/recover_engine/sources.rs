use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use super::RecoverEngine;
use super::offered::{PARTITION, SetRow};
use super::tables::Scoring;

pub const SOURCES_FILE: &str = "bin_sources.tsv";
pub const ENSEMBLE_FILE: &str = "ensemble_info.tsv";

struct Match {
    bin: usize,
    jaccard: f64,
}

impl RecoverEngine {
    pub(super) fn write_sources(&self, scoring: &Scoring, offered: &[SetRow]) -> Result<()> {
        let lengths = self.lengths();
        let owners = self
            .candidates
            .iter()
            .map(|set| {
                let mut owner = vec![u32::MAX; lengths.len()];
                for (bin, (_, contigs)) in set.bins.iter().enumerate() {
                    for contig in contigs {
                        owner[*contig] = bin as u32;
                    }
                }
                let bases = set
                    .bins
                    .iter()
                    .map(|(_, contigs)| contigs.iter().map(|contig| lengths[*contig]).sum())
                    .collect::<Vec<usize>>();
                (owner, bases)
            })
            .collect::<Vec<_>>();
        let path = Path::new(&self.output_directory).join(SOURCES_FILE);
        let mut sink = BufWriter::new(std::fs::File::create(path)?);
        let columns = self
            .candidates
            .iter()
            .map(|set| format!("\t{0}_bin\t{0}_jaccard", set.name))
            .collect::<String>();
        writeln!(sink, "bin\tcontigs\tbp\torigin\torigin_jaccard{columns}")?;
        let mut credit = HashMap::<&str, (usize, usize)>::new();
        for (name, contigs) in scoring.bins() {
            let bp = contigs.iter().map(|contig| lengths[*contig]).sum::<usize>();
            let matches = owners
                .iter()
                .map(|(owner, bases)| {
                    let mut shared = HashMap::<u32, usize>::new();
                    for contig in contigs {
                        if owner[*contig] != u32::MAX {
                            *shared.entry(owner[*contig]).or_default() += lengths[*contig];
                        }
                    }
                    shared
                        .into_iter()
                        .max_by_key(|(bin, held)| (*held, std::cmp::Reverse(*bin)))
                        .map(|(bin, held)| Match {
                            bin: bin as usize,
                            jaccard: held as f64 / (bp + bases[bin as usize] - held) as f64,
                        })
                })
                .collect::<Vec<_>>();
            let best = matches
                .iter()
                .enumerate()
                .filter_map(|(set, found)| Some((set, found.as_ref()?.jaccard)))
                .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)));
            let (origin, jaccard) = match best {
                Some((set, jaccard)) if jaccard > 0.5 => {
                    (self.candidates[set].name.as_str(), jaccard)
                }
                Some((_, jaccard)) => (PARTITION, jaccard),
                None => (PARTITION, 0.0),
            };
            let entry = credit.entry(origin).or_default();
            entry.0 += 1;
            entry.1 += bp;
            write!(
                sink,
                "{name}\t{}\t{bp}\t{origin}\t{jaccard:.4}",
                contigs.len()
            )?;
            for (set, found) in self.candidates.iter().zip(&matches) {
                match found {
                    Some(found) => {
                        write!(sink, "\t{}\t{:.4}", set.bins[found.bin].0, found.jaccard)?
                    }
                    None => write!(sink, "\t\t0.0000")?,
                }
            }
            writeln!(sink)?;
        }
        sink.flush()?;
        self.write_ensemble(offered, &credit)
    }

    fn write_ensemble(
        &self,
        offered: &[SetRow],
        credit: &HashMap<&str, (usize, usize)>,
    ) -> Result<()> {
        let path = Path::new(&self.output_directory).join(ENSEMBLE_FILE);
        let mut sink = BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            sink,
            "set\toffered_bins\toffered_bp\trefined_bins\trefined_bp\ttaken_whole\ttaken_trimmed\
             \torigin_bins\torigin_bp"
        )?;
        for row in offered {
            let refined = row
                .refined
                .map_or_else(|| "\t".to_string(), |(bins, bp)| format!("{bins}\t{bp}"));
            let (bins, bp) = credit.get(row.name.as_str()).copied().unwrap_or_default();
            writeln!(
                sink,
                "{}\t{}\t{}\t{refined}\t{}\t{}\t{bins}\t{bp}",
                row.name, row.offered.0, row.offered.1, row.whole, row.trimmed
            )?;
        }
        sink.flush()?;
        Ok(())
    }
}
