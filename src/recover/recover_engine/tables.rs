use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use log::warn;

use super::RecoverEngine;
use super::bin_writer::{Published, Target, Written};
use crate::quality::report::{Bin, Scored};

pub const ORIGINS_FILE: &str = "replicon_origins.tsv";
pub const STRAYS_FILE: &str = "bridge_strays.tsv";

pub(super) struct Scoring {
    bins: Vec<(String, Vec<usize>)>,
}

impl Scoring {
    pub(super) fn bins(&self) -> &[(String, Vec<usize>)] {
        &self.bins
    }

    pub(super) fn of(published: &Published, prefix: &str) -> Self {
        let bins = published
            .bins
            .iter()
            .map(|(bin, contigs)| (published.genome(*bin).label(prefix), contigs.clone()))
            .chain(
                published
                    .replicon_bins()
                    .into_iter()
                    .enumerate()
                    .map(|(at, contigs)| (Target::Replicon(at + 1).label(prefix), contigs)),
            )
            .collect();
        Self { bins }
    }
}

impl RecoverEngine {
    pub(super) fn reports_markers(&self, contig: usize) -> bool {
        self.reports.marker.path.is_some() && self.quality.hit_count(&[contig]) > 0
    }

    // Written from the bins that are written out, not from the refiner's last pass, so the
    // tables and the assignments never describe different partitions.
    pub(super) fn write_tables(&self, scoring: &Scoring, written: &Written) -> Result<()> {
        let names = &self.tables.coverage.contig_names;
        let members = scoring
            .bins
            .iter()
            .map(|(_, contigs)| contigs.as_slice())
            .collect::<Vec<_>>();
        let strain = self
            .annotator
            .strain_heterogeneity(&self.quality, &members, names)
            .unwrap_or_else(|error| {
                warn!("Strain heterogeneity: {error}");
                vec![None; members.len()]
            });
        let lineages = self
            .annotator
            .checkm_lineages(&members, names)
            .unwrap_or_else(|error| {
                warn!("CheckM lineage sets: {error}");
                vec![None; members.len()]
            });
        let bins = scoring
            .bins
            .iter()
            .zip(strain)
            .zip(lineages)
            .map(|(((name, contigs), strain), checkm)| Bin {
                name: name.clone(),
                contigs,
                strain,
                checkm,
            })
            .collect::<Vec<_>>();
        let directory = Path::new(&self.output_directory);
        let scored = Scored {
            markers: &self.quality,
            names,
            lengths: &self.tables.coverage.contig_lengths,
            bases: &written.bases,
        };
        if let Err(error) = scored.write(&bins, &directory.join(crate::defaults::QUALITY_FILE)) {
            warn!("Quality tables: {error}");
        }
        if let Err(error) = crate::recover::abundance::write(
            &scoring.bins,
            &self.tables.coverage,
            &directory.join(crate::recover::abundance::ABUNDANCE_FILE),
        ) {
            warn!("Abundance table: {error}");
        }
        if let Some(path) = &self.reports.marker.path {
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

    // A replicon keeps the bin it left, so a scorer that counts plasmids with their host can fold
    // it back without the bins changing.
    pub(super) fn write_origins(&self, published: &Published) -> Result<()> {
        let path = Path::new(&self.output_directory).join(ORIGINS_FILE);
        let mut sink = BufWriter::new(std::fs::File::create(path)?);
        writeln!(sink, "replicon\tcontig\torigin")?;
        for (at, (contig, origin)) in published
            .replicons
            .iter()
            .zip(&published.origins)
            .enumerate()
        {
            let origin =
                origin.map_or_else(String::new, |bin| published.genome(bin).label(&self.prefix));
            writeln!(
                sink,
                "{}\t{}\t{origin}",
                Target::Replicon(published.groups[at] + 1).label(&self.prefix),
                self.tables.coverage.contig_names[*contig]
            )?;
        }
        sink.flush()?;
        Ok(())
    }

    // The bins before the move are the written bins with each stray handed back, so one run
    // measures the move without a second run that may take another path.
    pub(super) fn write_strays(&self, published: &Published) -> Result<()> {
        let path = Path::new(&self.output_directory).join(STRAYS_FILE);
        let mut sink = BufWriter::new(std::fs::File::create(path)?);
        writeln!(sink, "contig\tfrom\tto")?;
        let label = |bin: usize| published.genome(bin).label(&self.prefix);
        for (contig, from, to) in &published.strays {
            writeln!(
                sink,
                "{}\t{}\t{}",
                self.tables.coverage.contig_names[*contig],
                from.map_or_else(|| Target::Unbinned.name(), label),
                label(*to)
            )?;
        }
        sink.flush()?;
        Ok(())
    }
}
