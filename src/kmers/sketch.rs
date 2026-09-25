use std::collections::HashMap;

use anyhow::{Result, bail};
use needletail::Sequence;
use needletail::bitkmer::BitNuclKmer;
use rayon::prelude::*;

/// needletail's `extend_kmer` masks with `2^(2k) - 1`, which overflows a u64 at k=32.
pub const DEFAULT_KMER_SIZE: u8 = 31;
pub const DEFAULT_SCALE: u64 = 200;

const BATCH: usize = 512;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SketchParams {
    pub kmer_size: u8,
    pub scale: u64,
}

impl Default for SketchParams {
    fn default() -> Self {
        Self {
            kmer_size: DEFAULT_KMER_SIZE,
            scale: DEFAULT_SCALE,
        }
    }
}

pub struct ContigSketches {
    params: SketchParams,
    contig_names: Vec<String>,
    hashes: Vec<u64>,
    offsets: Vec<u32>,
    occurrences: Vec<u32>,
}

pub fn sketch_sequence(sequence: &[u8], params: SketchParams) -> (Vec<u64>, u32) {
    let bar = u64::MAX / params.scale;
    let mut kept = Vec::new();
    for (_, kmer, _) in BitNuclKmer::new(sequence, params.kmer_size, true) {
        let hashed = crate::seeds::mix(kmer.0);
        if hashed <= bar {
            kept.push(hashed);
        }
    }
    let occurrences = kept.len() as u32;
    kept.sort_unstable();
    kept.dedup();
    (kept, occurrences)
}

impl ContigSketches {
    pub fn build(assembly: &str) -> Result<Self> {
        let mut reader = needletail::parse_fastx_file(assembly)?;
        let mut built = Self {
            params: SketchParams::default(),
            contig_names: Vec::new(),
            hashes: Vec::new(),
            offsets: vec![0],
            occurrences: Vec::new(),
        };
        let mut batch: Vec<(String, Vec<u8>)> = Vec::with_capacity(BATCH);
        while let Some(record) = reader.next() {
            let seqrec = record?;
            batch.push((
                crate::contig_id(seqrec.id())?.to_string(),
                seqrec.normalize(false).into_owned(),
            ));
            if batch.len() == BATCH {
                built.absorb(&mut batch);
            }
        }
        built.absorb(&mut batch);
        Ok(built)
    }

    fn absorb(&mut self, batch: &mut Vec<(String, Vec<u8>)>) {
        let params = self.params;
        let sketched = batch
            .par_iter()
            .map(|(_, sequence)| sketch_sequence(sequence, params))
            .collect::<Vec<_>>();
        for ((name, _), (hashes, occurrences)) in batch.drain(..).zip(sketched) {
            self.contig_names.push(name);
            self.hashes.extend_from_slice(&hashes);
            self.offsets.push(self.hashes.len() as u32);
            self.occurrences.push(occurrences);
        }
    }

    /// `ContigFeatures` is indexed by the coverage table's row order, so the sketch is rebuilt
    /// in that order rather than trusting the assembly and the filter to have agreed.
    pub fn align_to(&mut self, wanted: &[String]) -> Result<()> {
        let at = self
            .contig_names
            .iter()
            .enumerate()
            .map(|(index, name)| (name.as_str(), index))
            .collect::<HashMap<_, _>>();
        let mut hashes = Vec::with_capacity(self.hashes.len());
        let mut offsets = vec![0u32];
        let mut occurrences = Vec::with_capacity(wanted.len());
        for name in wanted {
            let Some(index) = at.get(name.as_str()) else {
                bail!("the assembly holds no contig named {name} to sketch");
            };
            hashes.extend_from_slice(self.run(*index));
            offsets.push(hashes.len() as u32);
            occurrences.push(self.occurrences[*index]);
        }
        self.hashes = hashes;
        self.offsets = offsets;
        self.occurrences = occurrences;
        self.contig_names = wanted.to_vec();
        Ok(())
    }

    fn run(&self, index: usize) -> &[u64] {
        let start = self.offsets[index] as usize;
        let stop = self.offsets[index + 1] as usize;
        &self.hashes[start..stop]
    }

    fn unique(&self, contigs: &[usize]) -> usize {
        let mut every = contigs
            .iter()
            .flat_map(|c| self.run(*c).iter().copied())
            .collect::<Vec<_>>();
        every.sort_unstable();
        every.dedup();
        every.len()
    }

    /// The scale divides out of both halves, so a subsampled sketch estimates the same
    /// fraction an exact count would.
    pub fn duplication(&self, contigs: &[usize]) -> Option<f64> {
        let total: u64 = contigs.iter().map(|c| self.occurrences[*c] as u64).sum();
        if total == 0 {
            return None;
        }
        Some(1.0 - self.unique(contigs) as f64 / total as f64)
    }
}
