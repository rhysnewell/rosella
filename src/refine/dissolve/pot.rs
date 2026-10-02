use std::collections::{HashMap, HashSet};

use crate::embedding::features::ContigFeatures;
use crate::quality::{Quality, Scorer};
use crate::refine::owners::owners;
use crate::refine::ranking::{remaining, remaining_in};
use crate::refine::rung::{Rung, Verdict, verdict};

pub struct Pot<'a> {
    features: &'a ContigFeatures<'a>,
    quality: &'a dyn Scorer,
    worth: f64,
    floor: usize,
    origin: HashMap<usize, usize>,
    members: HashMap<usize, Vec<usize>>,
}

impl<'a> Pot<'a> {
    pub fn new(
        features: &'a ContigFeatures<'a>,
        quality: &'a dyn Scorer,
        worth: f64,
        floor: usize,
        dissolved: &[(usize, Vec<usize>)],
    ) -> Self {
        Self {
            features,
            quality,
            worth,
            floor,
            origin: owners(dissolved.iter().map(|(bin_id, contigs)| (*bin_id, contigs))),
            members: dissolved.iter().cloned().collect(),
        }
    }

    pub fn length(&self, contig: usize) -> usize {
        self.features.length(contig)
    }

    pub fn bases(&self, contigs: &[usize]) -> usize {
        self.features.bin_size(contigs)
    }

    pub fn quality_of(&self, contigs: &[usize]) -> Quality {
        self.quality.score(contigs)
    }

    pub fn worth(&self, contigs: &[usize]) -> f64 {
        self.quality.score(contigs).score(self.worth)
    }

    pub fn worth_of(&self, quality: Quality) -> f64 {
        quality.score(self.worth)
    }

    pub fn verdict(&self, contigs: &[usize], quality: Quality, rung: Rung) -> Verdict {
        verdict(self.bases(contigs), quality, rung)
    }

    pub fn origins(&self, contigs: &[usize]) -> Vec<(usize, usize)> {
        let mut taken: HashMap<usize, usize> = HashMap::new();
        for contig in contigs {
            if let Some(bin) = self.origin.get(contig) {
                *taken.entry(*bin).or_default() += self.features.length(*contig);
            }
        }
        let mut taken = taken.into_iter().collect::<Vec<_>>();
        taken.sort_unstable();
        taken
    }

    // What each bin a candidate draws on still holds in the pool, scored once for both checks.
    pub fn standing(
        &self,
        contigs: &[usize],
        pool: &HashSet<usize>,
        claimed: &HashSet<usize>,
    ) -> Vec<Standing> {
        self.origins(contigs)
            .into_iter()
            .filter_map(|(bin, taken)| {
                let members = self.members.get(&bin)?;
                let contigs = remaining(&remaining_in(members, pool), claimed);
                Some(Standing {
                    taken,
                    held: self.quality.score(&contigs),
                    contigs,
                })
            })
            .collect()
    }

    // A candidate outscores the contigs it drains almost by construction, so the loss sits in the
    // remainder. The bin left behind must be no worse, unless the candidate is at least that good.
    pub fn conserves(&self, contigs: &[usize], worth: f64, standing: &[Standing]) -> bool {
        let taking = contigs.iter().copied().collect::<HashSet<_>>();
        standing.iter().all(|bin| {
            let before = self.worth_of(bin.held);
            worth >= before || self.worth(&remaining(&bin.contigs, &taking)) >= before
        })
    }

    // A strain half reads complete and clean, so worth cannot tell a genome carved in two from an
    // organism pulled out of a bin holding two. Only the second leaves the duplication behind.
    pub fn unifies(&self, quality: Quality, standing: &[Standing]) -> bool {
        let [bin] = standing else {
            return true;
        };
        if self.bases(&bin.contigs).saturating_sub(bin.taken) < self.floor {
            return true;
        }
        let doubled = bin.held.contamination;
        doubled > 0.0 && 2.0 * quality.contamination <= doubled
    }
}

pub struct Standing {
    taken: usize,
    held: Quality,
    contigs: Vec<usize>,
}
