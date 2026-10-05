use anyhow::{Result, bail};
use log::info;
use ndarray::s;

use crate::learned::{Learned, corpus::Inputs};
use crate::recover::recover_engine::RecoverEngine;

impl RecoverEngine {
    // A pass in k-mer space lends its near complete bins as labels, and the run then reads
    // composition from the space they teach.
    pub(super) fn learn_space(&mut self, contigs: &[usize]) -> Result<()> {
        if self.settings.features.distance.evidence.is_some() {
            bail!("--learned-space does not combine with a likelihood coverage distance");
        }
        let lengths = self.lengths();
        let cutoff = self.settings.lengths.cutoff;
        let long = contigs.partition_point(|contig| lengths[*contig] >= cutoff);
        let labels = self.first_labels(&contigs[..long])?;
        let tables = &self.tables;
        let learned = Learned::train(&Inputs {
            assembly: &tables.assembly,
            names: &tables.coverage.contig_names,
            lengths: &tables.coverage.contig_lengths,
            composition: tables.tnf.kmer_table.slice(s![..long, ..]),
            coverage: tables.coverage.table.slice(s![..long, ..]),
            labels: &labels,
            halves: &tables.halves,
            halves_floor: 2 * cutoff,
            kmer_size: tables.tnf.kmer_size(),
            seed: self.settings.seeds.seed,
        })?;
        let Some(learned) = learned else {
            return Ok(());
        };
        let rows = {
            let _timer = crate::timing::scope("learned_embed");
            learned.embed(tables.tnf.kmer_table.view(), tables.coverage.table.view())
        };
        info!(
            "Composition now reads {} learned dimensions for {} contigs.",
            rows.ncols(),
            rows.nrows()
        );
        self.tables.tnf.kmer_table = rows;
        self.tables.learned = Some(learned);
        Ok(())
    }
}
