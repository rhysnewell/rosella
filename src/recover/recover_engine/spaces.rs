use std::path::Path;

use anyhow::Result;
use log::debug;
use ndarray::Array2;

use crate::learned::Space;
use crate::recover::ladder::{Combined, Judge, combine_traced};
use crate::recover::recover_engine::RecoverEngine;
use crate::recover::recover_engine::weight::Settled;
use crate::tables::Tables;

pub(super) struct Aside {
    rows: Array2<f64>,
    space: Option<Space>,
}

impl Aside {
    fn swap(&mut self, tables: &mut Tables) {
        std::mem::swap(&mut self.rows, &mut tables.tnf.kmer_table);
        std::mem::swap(&mut self.space, &mut tables.learned);
    }
}

impl RecoverEngine {
    // The labelled space keeps sister strains apart that the windowed space fuses, and it loses
    // everywhere else, so it learns from the raw k-mer rows the windowed space replaced.
    pub(super) fn learn_aside(
        &mut self,
        raw: Array2<f64>,
        contigs: &[usize],
        dir: &Path,
    ) -> Result<()> {
        let mut aside = Aside {
            rows: raw,
            space: None,
        };
        aside.swap(&mut self.tables);
        let learned = self.learn_labelled(contigs, Some(dir));
        aside.swap(&mut self.tables);
        learned?;
        self.aside = aside.space.is_some().then_some(aside);
        Ok(())
    }

    // A bin fused from two sister strains carries every marker twice, so the per-bin combine
    // takes each genome from whichever space kept it whole.
    pub(super) fn partitioned_spaces(&mut self, contigs: &[usize]) -> Result<Settled> {
        let mut settled = self.partitioned(contigs)?;
        let Some(mut aside) = self.aside.take() else {
            return Ok(settled);
        };
        debug!("Labelled space: partitioning");
        aside.swap(&mut self.tables);
        let other = self.partitioned(contigs);
        aside.swap(&mut self.tables);
        drop(aside);
        let arms = [settled.partitioning, other?.partitioning];
        let long = contigs.len() - settled.parked.len();
        let features = self.features();
        let judge = Judge {
            quality: &self.quality,
            contigs: &contigs[..long],
            lengths: self.lengths(),
            bars: self.settings.bars,
            depth: Some(&features),
        };
        let Combined {
            partitioning,
            sources,
        } = combine_traced(&arms, &[0, 0], &judge, None);
        let [windowed, labelled, both] = sources.iter().fold([0; 3], |mut tally, (arms, _)| {
            tally[match arms.as_slice() {
                [0] => 0,
                [1] => 1,
                _ => 2,
            }] += 1;
            tally
        });
        debug!(
            "Spaces: {windowed} bins from the windowed space, {labelled} from the labelled, {both} from both"
        );
        if !settled.parked.is_empty() {
            settled.worth.kept = self.pass_worth(&partitioning, contigs);
        }
        settled.partitioning = partitioning;
        Ok(settled)
    }
}
