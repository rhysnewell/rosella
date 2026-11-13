use std::collections::{BTreeMap, BTreeSet, HashSet};

use anyhow::Result;
use log::debug;
use ndarray::s;

use crate::embedding::knn::KnnGraph;
use crate::markers::viruses::{Call, mobile};
use crate::quality::Scorer;
use crate::recover::coherence::Coherence;
use crate::recover::fold::Role;
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::adrift::scattered;
use crate::refine::owners::owners;

impl RecoverEngine {
    // A part genome reads clean on markers whatever it mixes, so below a genome's piece a bin
    // under medium quality stays only while its depths agree and its contigs are each other's
    // neighbours, or its genes read phage or plasmid.
    pub(super) fn incoherent(
        &self,
        bins: &BTreeMap<usize, Vec<usize>>,
        whole: &HashSet<usize>,
        knn: &KnnGraph,
    ) -> Result<Vec<usize>> {
        let Some(read_lengths) = self.tables.read_lengths.as_deref() else {
            return Ok(Vec::new());
        };
        let lengths = &self.tables.coverage.contig_lengths;
        let depths = self.tables.coverage.table.slice(s![.., ..;2]);
        let scores = bins
            .values()
            .map(|members| self.quality.score(members))
            .collect::<Vec<_>>();
        let roles = scores
            .iter()
            .map(|score| Role::of(*score, false))
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
        let owner = owners(bins.iter().map(|(bin, members)| (*bin, members)));
        let loose = |members: &[usize]| {
            !coherence.holds_one(members, depths, lengths, read_lengths)
                || scattered(members, &owner, knn)
        };
        let unsure = bins
            .iter()
            .zip(&scores)
            .filter(|((bin, members), score)| {
                Role::partial(**score)
                    && members.len() > 1
                    && !whole.contains(bin)
                    && members.iter().map(|contig| lengths[*contig]).sum::<usize>() < piece
                    && loose(members)
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
            "Part-genome bins unbinned on depth or neighbours: {}, kept as phage or plasmid: {}",
            unbinned.len(),
            kept.len()
        );
        Ok(unbinned)
    }

    // No eukaryote is as small as a prokaryote's piece, and depth and neighbours cannot tell
    // related eukaryotes apart, so a eukaryotic bin under a piece is a fragment of several.
    pub(super) fn fragments(
        &self,
        bins: &BTreeMap<usize, Vec<usize>>,
        eukaryotes: &HashSet<usize>,
        whole: &HashSet<usize>,
    ) -> Vec<usize> {
        let lengths = &self.tables.coverage.contig_lengths;
        let piece = self.settings.bars.piece();
        let fragments = bins
            .iter()
            .filter(|(bin, members)| {
                eukaryotes.contains(bin)
                    && !whole.contains(bin)
                    && members.iter().map(|contig| lengths[*contig]).sum::<usize>() < piece
            })
            .map(|(bin, _)| *bin)
            .collect::<Vec<_>>();
        debug!("Eukaryotic fragments unbinned: {}", fragments.len());
        fragments
    }
}
