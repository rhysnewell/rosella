pub mod pairs;

use std::collections::{HashMap, HashSet};

use anyhow::Result;
use log::info;

pub use pairs::Pair;

use crate::cli::BinningParams;
use crate::external::skani_engine::SkaniEngine;

/// Two contigs that align to each other over most of both are the same locus in two organisms,
/// which is the one claim depth and composition cannot make. Two loci of one genome do not align.
#[derive(Debug, Clone, Copy)]
pub struct HomologySettings {
    pub min_identity: f64,
    pub min_aligned_fraction: f64,
    pub min_contig_length: usize,
    pub min_pair_length: usize,
}

impl Default for HomologySettings {
    fn default() -> Self {
        Self {
            min_identity: 90.0,
            min_aligned_fraction: 50.0,
            min_contig_length: 1_500,
            min_pair_length: 0,
        }
    }
}

impl HomologySettings {
    fn admits(&self, pair: &Pair, longer: usize) -> bool {
        pair.identity >= self.min_identity
            && pair.aligned_fraction() >= self.min_aligned_fraction
            && longer >= self.min_pair_length
    }
}

#[derive(Debug, Default)]
pub struct Homology {
    apart: HashSet<(u32, u32)>,
    neighbours: HashMap<u32, Vec<u32>>,
}

impl Homology {
    pub fn build(
        assembly: &str,
        threads: usize,
        settings: HomologySettings,
        names: &[String],
        lengths: &[usize],
    ) -> Result<Self> {
        if !SkaniEngine::is_installed() {
            bail!(
                "skani is not on PATH. rosella compares contigs to each other through it, so \
                 either install it or drop --homology"
            );
        }
        let pairs = SkaniEngine::new(assembly, threads, settings.min_contig_length)
            .pairs(names, lengths)?;
        let homology = Self::from_pairs(pairs, lengths, settings);
        info!("{} contig pairs are homologous", homology.len());
        Ok(homology)
    }

    pub fn from_pairs(
        pairs: impl IntoIterator<Item = Pair>,
        lengths: &[usize],
        settings: HomologySettings,
    ) -> Self {
        let apart = pairs
            .into_iter()
            .filter(|pair| {
                let longer = lengths[pair.one].max(lengths[pair.other]);
                settings.admits(pair, longer)
            })
            .map(|pair| key(pair.one, pair.other))
            .collect::<HashSet<_>>();
        let mut neighbours: HashMap<u32, Vec<u32>> = HashMap::new();
        for (one, other) in &apart {
            neighbours.entry(*one).or_default().push(*other);
            neighbours.entry(*other).or_default().push(*one);
        }
        Self { apart, neighbours }
    }

    /// Evidence that a bin holds two organisms, which is a reason to cluster it again rather
    /// than a reason to cut it in a particular place.
    pub fn holds_pair(&self, contigs: &[usize]) -> bool {
        let members = contigs.iter().map(|at| *at as u32).collect::<HashSet<_>>();
        members.iter().any(|contig| {
            self.neighbours
                .get(contig)
                .is_some_and(|near| near.iter().any(|other| members.contains(other)))
        })
    }

    pub fn cannot_link(&self, one: usize, other: usize) -> bool {
        self.apart.contains(&key(one, other))
    }

    pub fn len(&self) -> usize {
        self.apart.len()
    }

    /// Greedy colouring of the cannot-link graph. Contigs that never align stay together, which
    /// is what standing every genome-sized contig alone gets wrong.
    pub fn groups(&self, contigs: &[usize]) -> Vec<Vec<usize>> {
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for contig in contigs {
            match groups
                .iter_mut()
                .find(|group| group.iter().all(|held| !self.cannot_link(*held, *contig)))
            {
                Some(group) => group.push(*contig),
                None => groups.push(vec![*contig]),
            }
        }
        groups
    }
}

fn key(one: usize, other: usize) -> (u32, u32) {
    let (low, high) = if one <= other {
        (one, other)
    } else {
        (other, one)
    };
    (low as u32, high as u32)
}

pub fn homology_settings(
    binning: &BinningParams,
    min_contig_length: usize,
) -> Option<HomologySettings> {
    (binning.homology || binning.homology_trigger).then_some(HomologySettings {
        min_identity: binning.homology_identity,
        min_aligned_fraction: binning.homology_aligned_fraction,
        min_contig_length,
        min_pair_length: binning.homology_min_length,
    })
}
