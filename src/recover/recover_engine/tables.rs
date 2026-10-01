use std::path::Path;

use anyhow::Result;
use log::warn;

use super::RecoverEngine;
use super::bin_writer::{BIN_PREFIX, Published, REPLICON_PREFIX, Written};
use crate::quality::report::{Bin, Scored};

pub(super) struct Scoring {
    bins: Vec<(String, Vec<usize>)>,
    genome_contigs: Vec<usize>,
}

impl Scoring {
    pub(super) fn of(published: &Published) -> Self {
        let genome_contigs = published.bins.values().flatten().copied().collect();
        let bins = published
            .bins
            .iter()
            .map(|(bin, contigs)| (format!("{BIN_PREFIX}{bin}"), contigs.clone()))
            .chain(published.replicons.iter().enumerate().map(|(at, contig)| {
                (
                    format!("{BIN_PREFIX}{REPLICON_PREFIX}{}", at + 1),
                    vec![*contig],
                )
            }))
            .collect();
        Self {
            bins,
            genome_contigs,
        }
    }
}

impl RecoverEngine {
    pub(super) fn reports_markers(&self, contig: usize) -> bool {
        self.marker_report.is_some() && self.quality.hit_count(&[contig]) > 0
    }

    // Written from the bins that are written out, not from the refiner's last pass, so the
    // tables and the assignments never describe different partitions.
    pub(super) fn write_tables(&mut self, scoring: &Scoring, written: &Written) -> Result<()> {
        let names = &self.coverage_table.contig_names;
        if let Err(error) =
            self.annotator
                .complete_checkm(&mut self.quality, &scoring.genome_contigs, names)
        {
            warn!("Could not search the CheckM models: {error}");
        }
        let members = scoring
            .bins
            .iter()
            .map(|(_, contigs)| contigs.as_slice())
            .collect::<Vec<_>>();
        let strain = self
            .annotator
            .strain_heterogeneity(&self.quality, &members, names)
            .unwrap_or_else(|error| {
                warn!("Could not measure strain heterogeneity: {error}");
                vec![None; members.len()]
            });
        let bins = scoring
            .bins
            .iter()
            .zip(strain)
            .map(|((name, contigs), strain)| Bin {
                name: name.clone(),
                contigs,
                strain,
            })
            .collect::<Vec<_>>();
        let directory = Path::new(&self.output_directory);
        let scored = Scored {
            markers: &self.quality,
            names,
            lengths: &self.coverage_table.contig_lengths,
            bases: &written.bases,
        };
        if let Err(error) = scored.write(&bins, &directory.join(crate::defaults::QUALITY_FILE)) {
            warn!("Could not write the quality tables: {error}");
        }
        if let Err(error) = crate::recover::abundance::write(
            &scoring.bins,
            &self.coverage_table,
            &directory.join(crate::recover::abundance::ABUNDANCE_FILE),
        ) {
            warn!("Could not write the abundance table: {error}");
        }
        if let Some(path) = &self.marker_report {
            let mut placed = written
                .labels
                .iter()
                .map(|(contig, label)| (label.as_str(), *contig))
                .collect::<Vec<_>>();
            placed.sort_unstable_by_key(|(_, contig)| *contig);
            self.quality.report(placed, names, path)?;
        }
        Ok(())
    }
}
