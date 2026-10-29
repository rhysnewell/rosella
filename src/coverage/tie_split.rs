use std::collections::{HashMap, HashSet};
use std::io::Read;

use anyhow::{Context, Result, bail};
use log::debug;
use noodles_bam as bam;
use noodles_sam::alignment::record::{Cigar as _, data::field::Tag};
use rayon::prelude::*;

use super::bam_depth::Pileup;
use super::depth_file::Contigs;
use super::read_stats::ReadStats;
use crate::digest::fold;

const ROUNDS: usize = 200;
const SETTLED: f64 = 1e-5;
// A fixed count, not the thread count, so the sums add in the same order on any machine.
const CHUNKS: usize = 16;

#[derive(Clone, Copy)]
struct Hit {
    contig: u32,
    start: u32,
    span: u32,
}

#[derive(Default)]
struct Mate {
    primary: bool,
    hits: Vec<(i64, Hit)>,
}

#[derive(Default)]
struct Ties {
    unique: Vec<Hit>,
    tied: Vec<Hit>,
    ends: Vec<usize>,
}

impl Ties {
    fn flush(&mut self, mate: &mut Mate, stats: &mut ReadStats) -> Result<()> {
        if mate.hits.is_empty() {
            return Ok(());
        }
        if !std::mem::take(&mut mate.primary) {
            bail!("a read's alignments are not next to each other. Keep the mapper's order");
        }
        let top = mate.hits.iter().map(|(score, _)| *score).max();
        let from = self.tied.len();
        for (score, hit) in mate.hits.drain(..) {
            if Some(score) == top
                && !self.tied[from..]
                    .iter()
                    .any(|held| held.contig == hit.contig)
            {
                self.tied.push(hit);
            }
        }
        if self.tied.len() - from == 1 {
            self.unique.extend(self.tied.pop());
        } else {
            stats.tie(self.tied[from..].iter().map(|hit| hit.contig));
            self.ends.push(self.tied.len());
        }
        Ok(())
    }

    fn settle(&mut self, mates: &mut [Mate; 2], stats: &mut ReadStats) -> Result<()> {
        pair(mates);
        mates
            .iter_mut()
            .try_for_each(|mate| self.flush(mate, stats))
    }

    fn read(&self, at: usize) -> &[Hit] {
        let start = at.checked_sub(1).map_or(0, |before| self.ends[before]);
        &self.tied[start..self.ends[at]]
    }
}

// A fragment lies on one contig, so a mate that ties is placed by the contigs its partner allows.
fn pair(mates: &mut [Mate; 2]) {
    if mates.iter().any(|mate| mate.hits.is_empty()) {
        return;
    }
    let best = |mate: &Mate, contig: u32| {
        mate.hits
            .iter()
            .filter(|(_, hit)| hit.contig == contig)
            .map(|(score, _)| *score)
            .max()
    };
    let mut joint: Vec<(u32, i64)> = Vec::new();
    for (_, hit) in &mates[0].hits {
        if joint.iter().any(|(contig, _)| *contig == hit.contig) {
            continue;
        }
        if let (Some(first), Some(second)) =
            (best(&mates[0], hit.contig), best(&mates[1], hit.contig))
        {
            joint.push((hit.contig, first + second));
        }
    }
    let Some(top) = joint.iter().map(|(_, score)| *score).max() else {
        return;
    };
    for mate in mates.iter_mut() {
        let kept = joint
            .iter()
            .filter(|(_, score)| *score == top)
            .filter_map(|&(contig, _)| {
                let own = best(mate, contig);
                mate.hits
                    .iter()
                    .find(|(score, hit)| hit.contig == contig && Some(*score) == own)
                    .map(|&(_, hit)| (top, hit))
            })
            .collect();
        mate.hits = kept;
    }
}

pub(super) fn tally<'a, R: Read>(
    mut reader: bam::io::Reader<R>,
    contigs: &'a Contigs,
    index: &'a HashMap<&'a str, u32>,
) -> Result<(Pileup<'a>, ReadStats<'a>)> {
    let mut ties = Ties::default();
    let mut stats = ReadStats::new(contigs, index);
    let mut mates = [Mate::default(), Mate::default()];
    let mut name = Vec::new();
    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        let held: &[u8] = record.name().map_or(b"", |name| name.as_ref());
        if held != name.as_slice() {
            ties.settle(&mut mates, &mut stats)?;
            name.clear();
            name.extend_from_slice(held);
        }
        stats.note(&record)?;
        let Some((score, hit)) = hit_of(&record)? else {
            continue;
        };
        let flags = record.flags();
        let mate = &mut mates[usize::from(flags.is_last_segment())];
        mate.primary |= !flags.is_secondary();
        mate.hits.push((score, hit));
    }
    ties.settle(&mut mates, &mut stats)?;
    let depth = fit(&ties, &contigs.lengths);
    Ok((lay(ties, &depth, contigs), stats))
}

// A sorted BAM scatters a read's hits, so a first pass finds the reads with a secondary hit and
// only theirs are gathered by name. A read without one has a single hit and is unique.
pub(super) fn sorted<'a, R: Read>(
    mut open: impl FnMut() -> Result<bam::io::Reader<R>>,
    contigs: &'a Contigs,
    index: &'a HashMap<&'a str, u32>,
) -> Result<(Pileup<'a>, ReadStats<'a>)> {
    let mut record = bam::Record::default();
    let mut tied = HashSet::new();
    let mut reader = open()?;
    while reader.read_record(&mut record)? != 0 {
        if record.flags().is_secondary() {
            tied.insert(named(&record));
        }
    }
    let mut ties = Ties::default();
    let mut stats = ReadStats::new(contigs, index);
    let mut held = Vec::new();
    let mut reader = open()?;
    while reader.read_record(&mut record)? != 0 {
        stats.note(&record)?;
        let Some((score, hit)) = hit_of(&record)? else {
            continue;
        };
        let (name, flags) = (named(&record), record.flags());
        match tied.contains(&name) {
            true => held.push((
                name,
                flags.is_last_segment(),
                !flags.is_secondary(),
                score,
                hit,
            )),
            false => ties.unique.push(hit),
        }
    }
    drop(tied);
    held.sort_by_key(|(name, ..)| *name);
    for read in held.chunk_by(|first, second| first.0 == second.0) {
        let mut mates = [Mate::default(), Mate::default()];
        for &(_, last, primary, score, hit) in read {
            let mate = &mut mates[usize::from(last)];
            mate.primary |= primary;
            mate.hits.push((score, hit));
        }
        ties.settle(&mut mates, &mut stats)?;
    }
    let depth = fit(&ties, &contigs.lengths);
    Ok((lay(ties, &depth, contigs), stats))
}

fn named(record: &bam::Record) -> u64 {
    fold(record.name().map_or(b"", |name| name.as_ref()))
}

fn hit_of(record: &bam::Record) -> Result<Option<(i64, Hit)>> {
    let flags = record.flags();
    if flags.is_unmapped() || flags.is_supplementary() {
        return Ok(None);
    }
    let (Some(contig), Some(start)) = (
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
    ) else {
        return Ok(None);
    };
    let score = record
        .data()
        .get(&Tag::ALIGNMENT_SCORE)
        .transpose()?
        .and_then(|value| value.as_int())
        .context("an alignment carries no AS score to rank a read's hits by")?;
    Ok(Some((
        score,
        Hit {
            contig: contig as u32,
            start: usize::from(start) as u32 - 1,
            span: record.cigar().alignment_span()? as u32,
        },
    )))
}

// A strain at 10x should take ten times the tied reads of one at 1x.
fn fit(ties: &Ties, lengths: &[usize]) -> Vec<f64> {
    let mut unique = vec![0.0; lengths.len()];
    for hit in &ties.unique {
        unique[hit.contig as usize] += f64::from(hit.span);
    }
    let per_length = |bp: Vec<f64>| {
        bp.into_iter()
            .zip(lengths)
            .map(|(bp, length)| bp / *length as f64)
            .collect::<Vec<_>>()
    };
    let held = Held::of(ties, lengths.len());
    let mut depth = per_length(spread(ties, &unique, None, &held));
    for round in 0..ROUNDS {
        let next = per_length(spread(ties, &unique, Some(&depth), &held));
        let moved = next
            .iter()
            .zip(&depth)
            .map(|(next, depth)| (next - depth).abs())
            .sum::<f64>()
            / next.iter().sum::<f64>();
        depth = next;
        if moved.is_nan() || moved < SETTLED {
            debug!("Tie split rounds: {}", round + 1);
            break;
        }
    }
    depth
}

fn shares<'a>(hits: &'a [Hit], depth: Option<&'a [f64]>) -> impl Iterator<Item = f64> + 'a {
    let total = depth.map_or(0.0, |depth| {
        hits.iter().map(|hit| depth[hit.contig as usize]).sum()
    });
    hits.iter().map(move |hit| match depth {
        Some(depth) if total > 0.0 => depth[hit.contig as usize] / total,
        _ => 1.0 / hits.len() as f64,
    })
}

// Tied reads land on few contigs, so each pass sums onto those alone.
struct Held {
    contigs: Vec<u32>,
    slot: Vec<u32>,
}

impl Held {
    fn of(ties: &Ties, contigs: usize) -> Self {
        let mut held = ties.tied.iter().map(|hit| hit.contig).collect::<Vec<_>>();
        held.sort_unstable();
        held.dedup();
        let mut slot = vec![u32::MAX; contigs];
        for (at, contig) in held.iter().enumerate() {
            slot[*contig as usize] = at as u32;
        }
        Self {
            contigs: held,
            slot,
        }
    }
}

fn spread(ties: &Ties, unique: &[f64], depth: Option<&[f64]>, held: &Held) -> Vec<f64> {
    let per = ties.ends.len().div_ceil(CHUNKS).max(1);
    let parts = (0..CHUNKS)
        .into_par_iter()
        .map(|chunk| {
            let mut bp = vec![0.0; held.contigs.len()];
            for at in (chunk * per).min(ties.ends.len())..((chunk + 1) * per).min(ties.ends.len()) {
                let hits = ties.read(at);
                for (hit, share) in hits.iter().zip(shares(hits, depth)) {
                    bp[held.slot[hit.contig as usize] as usize] += f64::from(hit.span) * share;
                }
            }
            bp
        })
        .collect::<Vec<_>>();
    let mut bp = unique.to_vec();
    for part in parts {
        for (contig, part) in held.contigs.iter().zip(part) {
            bp[*contig as usize] += part;
        }
    }
    bp
}

fn lay<'a>(mut ties: Ties, depth: &[f64], contigs: &'a Contigs) -> Pileup<'a> {
    let weights = (0..ties.ends.len())
        .flat_map(|at| shares(ties.read(at), Some(depth)))
        .collect::<Vec<_>>();
    let mut starts = vec![0usize; contigs.lengths.len() + 1];
    for hit in &ties.tied {
        starts[hit.contig as usize + 1] += 1;
    }
    for contig in 0..contigs.lengths.len() {
        starts[contig + 1] += starts[contig];
    }
    let mut order = vec![0u32; ties.tied.len()];
    let mut next = starts.clone();
    for (at, hit) in ties.tied.iter().enumerate() {
        order[next[hit.contig as usize]] = at as u32;
        next[hit.contig as usize] += 1;
    }
    ties.unique.sort_by_key(|hit| hit.contig);

    let mut pileup = Pileup::new(contigs);
    let mut unique = ties.unique.iter().peekable();
    for contig in 0..contigs.lengths.len() {
        while let Some(hit) = unique.next_if(|hit| hit.contig as usize == contig) {
            pileup.add(contig, hit.start as usize, hit.span as usize, 1.0);
        }
        for &at in &order[starts[contig]..starts[contig + 1]] {
            let (hit, weight) = (ties.tied[at as usize], weights[at as usize]);
            pileup.add(contig, hit.start as usize, hit.span as usize, weight);
        }
    }
    pileup
}
