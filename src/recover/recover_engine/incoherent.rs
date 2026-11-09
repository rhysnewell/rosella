use std::collections::{BTreeMap, BTreeSet, HashSet};

use anyhow::Result;
use log::debug;
use ndarray::s;

use crate::markers::viruses::{Call, mobile};
use crate::quality::Scorer;
use crate::recover::coherence::Coherence;
use crate::recover::fold::Role;
use crate::recover::recover_engine::RecoverEngine;

impl RecoverEngine {
    // Markers cannot see a mix of genomes no marker reaches, so below a genome's piece a
    // marker-free bin stays only while its depths agree or its genes read phage or plasmid.
    pub(super) fn incoherent(
        &self,
        bins: &BTreeMap<usize, Vec<usize>>,
        whole: &HashSet<usize>,
    ) -> Result<Vec<usize>> {
        let Some(read_lengths) = self.tables.read_lengths.as_deref() else {
            return Ok(Vec::new());
        };
        if !self.settings.depth_unbin {
            return Ok(Vec::new());
        }
        let lengths = &self.tables.coverage.contig_lengths;
        let depths = self.tables.coverage.table.slice(s![.., ..;2]);
        let roles = bins
            .values()
            .map(|members| Role::of(self.quality.score(members), false))
            .collect::<Vec<_>>();
        let hosts = bins
            .values()
            .zip(&roles)
            .filter(|(_, role)| **role == Role::Host)
            .map(|(members, _)| members.as_slice())
            .collect::<Vec<_>>();
        let Some(coherence) = Coherence::of(&hosts, depths, lengths) else {
            return Ok(Vec::new());
        };
        let piece = self.settings.bars.piece();
        let unsure = bins
            .iter()
            .zip(&roles)
            .filter(|((bin, members), role)| {
                **role == Role::Donor
                    && members.len() > 1
                    && !whole.contains(bin)
                    && members.iter().map(|contig| lengths[*contig]).sum::<usize>() < piece
                    && !coherence.holds_one(members, depths, lengths, read_lengths)
            })
            .map(|((bin, _), _)| *bin)
            .collect::<Vec<_>>();
        let wanted = unsure
            .iter()
            .flat_map(|bin| bins[bin].iter().copied())
            .collect::<BTreeSet<_>>();
        let calls = self
            .annotator
            .virus_calls(&self.tables.coverage.contig_names, &wanted)?;
        let element = |bin: &usize| {
            mobile(
                bins[bin]
                    .iter()
                    .filter_map(|contig| calls.get(*contig).copied())
                    .fold(Call::default(), |total, call| total + call),
            )
        };
        let (kept, unbinned): (Vec<_>, Vec<_>) = unsure.into_iter().partition(element);
        debug!(
            "Marker-free bins unbinned on depth: {}, kept as phage or plasmid: {}",
            unbinned.len(),
            kept.len()
        );
        Ok(unbinned)
    }
}
