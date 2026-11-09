use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use anyhow::Result;
use log::debug;
use ndarray::s;

use crate::markers::viruses::{Call, mobile};
use crate::quality::Scorer;
use crate::recover::coherence::Coherence;
use crate::recover::copies::Copies;
use crate::recover::fold::Role;
use crate::recover::recover_engine::RecoverEngine;

impl RecoverEngine {
    // Markers cannot see a mix of genomes no marker reaches, so below a genome's piece a
    // marker-free bin stays only while its depths agree or its genes read phage or plasmid.
    pub(super) fn incoherent(
        &self,
        bins: &BTreeMap<usize, Vec<usize>>,
        whole: &HashSet<usize>,
    ) -> Result<Vec<(usize, Vec<Vec<usize>>)>> {
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
        let copies = self.copies(bins, &unbinned, read_lengths)?;
        debug!(
            "Marker-free bins unbinned on depth: {}, kept as phage or plasmid: {}, copy sets \
             kept: {}",
            unbinned.len(),
            kept.len(),
            copies.iter().map(Vec::len).sum::<usize>()
        );
        Ok(unbinned.into_iter().zip(copies).collect())
    }

    fn copies(
        &self,
        bins: &BTreeMap<usize, Vec<usize>>,
        unbinned: &[usize],
        read_lengths: &[f64],
    ) -> Result<Vec<Vec<Vec<usize>>>> {
        let lengths = &self.tables.coverage.contig_lengths;
        let depths = self.tables.coverage.table.slice(s![.., ..;2]);
        let wanted = unbinned
            .iter()
            .flat_map(|bin| bins[bin].iter().copied())
            .collect::<BTreeSet<_>>();
        if wanted.is_empty() {
            return Ok(Vec::new());
        }
        let ties = self.tables.ties_among(&wanted)?;
        let reads = wanted
            .iter()
            .map(|contig| {
                let length = lengths[*contig] as f64;
                let reads = depths
                    .row(*contig)
                    .iter()
                    .zip(read_lengths)
                    .map(|(depth, read_length)| depth * length / read_length)
                    .sum::<f64>();
                (*contig, reads)
            })
            .collect::<HashMap<_, _>>();
        let copies = Copies {
            lengths,
            reads: &reads,
            ties: &ties,
            read_length: read_lengths.iter().copied().fold(0.0, f64::max),
        };
        let sets = unbinned
            .iter()
            .map(|bin| copies.of(&bins[bin]))
            .collect::<Vec<_>>();
        let held = sets
            .iter()
            .flatten()
            .flatten()
            .copied()
            .collect::<BTreeSet<_>>();
        let closed = match held.is_empty() {
            true => HashSet::new(),
            false => self.closed(&held)?,
        };
        let min_bin = self.settings.lengths.min_bin;
        Ok(unbinned
            .iter()
            .zip(sets)
            .map(|(bin, sets)| {
                // A contig that closes as a circle is written as a bin of its own.
                let sets = match sets.iter().flatten().any(|contig| closed.contains(contig)) {
                    true => copies.of(&bins[bin]
                        .iter()
                        .copied()
                        .filter(|contig| !closed.contains(contig))
                        .collect::<Vec<_>>()),
                    false => sets,
                };
                sets.into_iter()
                    .filter(|set| {
                        set.iter().map(|contig| lengths[*contig]).sum::<usize>() >= min_bin
                    })
                    .collect()
            })
            .collect())
    }

    fn closed(&self, wanted: &BTreeSet<usize>) -> Result<HashSet<usize>> {
        let names = &self.tables.coverage.contig_names;
        let rows = wanted
            .iter()
            .map(|contig| (names[*contig].as_str(), *contig))
            .collect::<HashMap<_, _>>();
        let mut closed = HashSet::new();
        crate::kmers::each_named(
            &self.tables.assembly,
            |name| rows.contains_key(name),
            |name, record| {
                if self.circular(name, record.id(), &record.seq()) {
                    closed.insert(rows[name]);
                }
                Ok(())
            },
        )?;
        Ok(closed)
    }
}
