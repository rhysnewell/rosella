use std::collections::{BTreeMap, BTreeSet, HashSet};

use anyhow::Result;
use log::debug;
use ndarray::{Array2, s};

use crate::embedding::knn::KnnGraph;
use crate::markers::viruses::{Call, mobile};
use crate::quality::Scorer;
use crate::recover::coherence::Coherence;
use crate::recover::fold::{FITS, Role};
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::adrift::scattered;
use crate::refine::owners::owners;

struct Wander {
    gc: Array2<f64>,
    coherence: Coherence,
}

fn measured(members: &[usize], gc: &Array2<f64>) -> Vec<usize> {
    members
        .iter()
        .copied()
        .filter(|contig| gc[[*contig, 0]].is_finite())
        .collect()
}

impl Wander {
    fn fit(&self, members: &[usize], lengths: &[usize]) -> f64 {
        let members = measured(members, &self.gc);
        self.coherence.fit_of(&members, self.gc.view(), lengths)
    }
}

impl RecoverEngine {
    // A part genome reads clean on markers whatever it mixes. Below a genome's piece it stays
    // only while its depths agree and its contigs are neighbours, and past the piece only while
    // its GC keeps to one genome's spread. Phage or plasmid genes keep either.
    pub(super) fn incoherent(
        &self,
        bins: &BTreeMap<usize, Vec<usize>>,
        whole: &HashSet<usize>,
        eukaryotes: &HashSet<usize>,
        knn: &KnnGraph,
    ) -> Result<Vec<usize>> {
        let names = &self.tables.coverage.contig_names;
        let lengths = &self.tables.coverage.contig_lengths;
        let depths = self.tables.coverage.table.slice(s![.., ..;2]);
        let scores = bins
            .values()
            .map(|members| self.quality.score(members))
            .collect::<Vec<_>>();
        let hosts = bins
            .values()
            .zip(&scores)
            .filter(|(_, score)| Role::of(**score, false) == Role::Host)
            .map(|(members, _)| members.as_slice())
            .collect::<Vec<_>>();
        let depth = self.tables.read_lengths.as_deref().and_then(|reads| {
            Coherence::of(&hosts, depths, lengths).map(|coherence| (coherence, reads))
        });
        let wander = match self.settings.gc_unbin {
            true => self.wander(bins, &hosts),
            false => None,
        };
        if depth.is_none() && wander.is_none() {
            return Ok(Vec::new());
        }
        let piece = self.settings.bars.piece();
        let owner = owners(bins.iter().map(|(bin, members)| (*bin, members)));
        let loose = |members: &[usize]| {
            depth.as_ref().is_some_and(|(coherence, reads)| {
                !coherence.holds_one(members, depths, lengths, reads)
                    || scattered(members, &owner, knn)
            })
        };
        let small = |members: &[usize]| {
            members.iter().map(|contig| lengths[*contig]).sum::<usize>() < piece
        };
        let tested = bins
            .iter()
            .zip(&scores)
            .filter(|((bin, members), score)| {
                Role::partial(**score) && members.len() > 1 && !whole.contains(bin)
            })
            .collect::<Vec<_>>();
        let fits = wander.as_ref().map_or_else(BTreeMap::new, |wander| {
            tested
                .iter()
                .filter(|((bin, members), _)| !small(members) && !eukaryotes.contains(bin))
                .map(|((bin, members), _)| (**bin, wander.fit(members, lengths)))
                .collect::<BTreeMap<_, _>>()
        });
        if let Some(report) = &self.stages {
            let rows = tested
                .iter()
                .filter_map(|((bin, members), score)| Some((bin, members, score, fits.get(bin)?)))
                .map(|(bin, members, score, fit)| {
                    format!(
                        "{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.3e}",
                        bin,
                        names[members[0]],
                        members.len(),
                        members.iter().map(|contig| lengths[*contig]).sum::<usize>(),
                        score.completeness,
                        score.contamination,
                        fit
                    )
                });
            report.table(
                "gc_spread",
                "bin\tfirst\tcontigs\tbases\tcompleteness\tcontamination\tfit",
                rows,
            );
        }
        let unsure = tested
            .iter()
            .filter(|((bin, members), _)| match small(members) {
                true => loose(members),
                false => fits.get(bin).is_some_and(|fit| *fit < FITS),
            })
            .map(|((bin, _), _)| **bin)
            .collect::<Vec<_>>();
        let wanted = unsure
            .iter()
            .flat_map(|bin| bins[bin].iter().copied())
            .collect::<BTreeSet<_>>();
        let calls = self.annotator.virus_calls(names, &wanted)?;
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
            "Part-genome bins unbinned on depth, neighbours or GC: {}, kept as phage or plasmid: {}",
            unbinned.len(),
            kept.len()
        );
        Ok(unbinned)
    }

    // A contig too short to be counted in halves carries no GC here and sits the test out.
    fn wander(&self, bins: &BTreeMap<usize, Vec<usize>>, hosts: &[&[usize]]) -> Option<Wander> {
        let names = &self.tables.coverage.contig_names;
        let mut gc = Array2::from_elem((names.len(), 1), f64::NAN);
        for contig in bins.values().flatten() {
            if let Some(share) = self.tables.halves.gc(&names[*contig]) {
                gc[[*contig, 0]] = share;
            }
        }
        let hosts = hosts
            .iter()
            .map(|members| measured(members, &gc))
            .collect::<Vec<_>>();
        let hosts = hosts.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let lengths = &self.tables.coverage.contig_lengths;
        let coherence = Coherence::widest(&hosts, gc.view(), lengths)?;
        Some(Wander { gc, coherence })
    }
}
