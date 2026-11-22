use ndarray::{ArrayView2, Axis};

use crate::markers::replicon::DEPTH_FLOOR;
use crate::quality::FAMILY_ERROR;
use crate::recover::fold::{FITS, chi_square_above, median};

const FIRST_OCTAVE: usize = 1_000;
const LEAST: usize = 3;

// A short contig's depth scatters more than a long one's, so each length octave is judged on the
// spread near complete genomes show at that length.
pub struct Coherence {
    spreads: Vec<Option<Vec<f64>>>,
    scale: fn(f64) -> f64,
}

fn log_depth(depth: f64) -> f64 {
    (depth + DEPTH_FLOOR).ln()
}

fn widest(mut spreads: Vec<f64>) -> f64 {
    spreads.sort_unstable_by(|a, b| b.total_cmp(a));
    spreads[(FAMILY_ERROR * spreads.len() as f64) as usize]
}

fn octave(length: usize) -> usize {
    (length.max(FIRST_OCTAVE) / FIRST_OCTAVE).ilog2() as usize
}

fn scaled(
    contig: usize,
    values: ArrayView2<'_, f64>,
    scale: fn(f64) -> f64,
) -> impl Iterator<Item = f64> + '_ {
    values
        .index_axis_move(Axis(0), contig)
        .into_iter()
        .map(move |value| scale(*value))
}

fn centre(
    members: &[usize],
    values: ArrayView2<f64>,
    lengths: &[usize],
    scale: fn(f64) -> f64,
) -> Vec<f64> {
    let mut mean = vec![0.0; values.ncols()];
    let total = members.iter().map(|contig| lengths[*contig]).sum::<usize>() as f64;
    for contig in members {
        let weight = lengths[*contig] as f64 / total;
        for (sum, value) in mean.iter_mut().zip(scaled(*contig, values, scale)) {
            *sum += weight * value;
        }
    }
    mean
}

impl Coherence {
    pub fn of(hosts: &[&[usize]], depths: ArrayView2<f64>, lengths: &[usize]) -> Option<Self> {
        Self::fitted(hosts, depths, lengths, log_depth, median)
    }

    // No read count explains how far a genome's GC wanders, and some wander further than
    // others, so the spread is the one all but the family error of hosts sit under.
    pub fn widest(hosts: &[&[usize]], values: ArrayView2<f64>, lengths: &[usize]) -> Option<Self> {
        Self::fitted(hosts, values, lengths, std::convert::identity, widest)
    }

    fn fitted(
        hosts: &[&[usize]],
        values: ArrayView2<f64>,
        lengths: &[usize],
        scale: fn(f64) -> f64,
        typical: fn(Vec<f64>) -> f64,
    ) -> Option<Self> {
        let samples = values.ncols();
        let mut held: Vec<Vec<Vec<f64>>> = Vec::new();
        for members in hosts {
            let mean = centre(members, values, lengths, scale);
            let mut sums: Vec<(Vec<f64>, usize)> = Vec::new();
            for contig in *members {
                let at = octave(lengths[*contig]);
                if sums.len() <= at {
                    sums.resize(at + 1, (vec![0.0; samples], 0));
                }
                let (sum, count) = &mut sums[at];
                let read = scaled(*contig, values, scale);
                for ((total, value), centre) in sum.iter_mut().zip(read).zip(&mean) {
                    *total += (value - centre).powi(2);
                }
                *count += 1;
            }
            for (at, (sum, count)) in sums.into_iter().enumerate() {
                if count < LEAST {
                    continue;
                }
                if held.len() <= at {
                    held.resize(at + 1, Vec::new());
                }
                held[at].push(sum.into_iter().map(|total| total / count as f64).collect());
            }
        }
        let spreads = held
            .into_iter()
            .map(|bins| {
                (bins.len() >= LEAST).then(|| {
                    (0..samples)
                        .map(|sample| typical(bins.iter().map(|bin| bin[sample]).collect()))
                        .collect()
                })
            })
            .collect::<Vec<_>>();
        spreads
            .iter()
            .any(Option::is_some)
            .then_some(Self { spreads, scale })
    }

    pub(crate) fn spread(&self, length: usize) -> &[f64] {
        let at = octave(length) as isize;
        self.spreads
            .iter()
            .enumerate()
            .filter_map(|(octave, spread)| spread.as_ref().map(|spread| (octave as isize, spread)))
            .min_by_key(|(octave, _)| ((octave - at).abs(), *octave))
            .map(|(_, spread)| spread.as_slice())
            .unwrap_or_default()
    }

    // A contig at low depth holds few reads, so on top of the genome spread each term carries
    // the count noise its expected reads give: depth times read length over contig length.
    pub fn holds_one(
        &self,
        members: &[usize],
        depths: ArrayView2<f64>,
        lengths: &[usize],
        read_lengths: &[f64],
    ) -> bool {
        let total = members.iter().map(|contig| lengths[*contig]).sum::<usize>() as f64;
        let expected = (0..depths.ncols())
            .map(|sample| {
                members
                    .iter()
                    .map(|contig| depths[[*contig, sample]] * lengths[*contig] as f64)
                    .sum::<f64>()
                    / total
            })
            .collect::<Vec<_>>();
        let count = |contig: usize, sample: usize| {
            let depth = expected[sample];
            depth * read_lengths[sample] / lengths[contig] as f64 / (depth + DEPTH_FLOOR).powi(2)
        };
        self.fit(members, depths, lengths, count) >= FITS
    }

    pub fn fit_of(&self, members: &[usize], values: ArrayView2<f64>, lengths: &[usize]) -> f64 {
        self.fit(members, values, lengths, |_, _| 0.0)
    }

    fn fit(
        &self,
        members: &[usize],
        values: ArrayView2<f64>,
        lengths: &[usize],
        noise: impl Fn(usize, usize) -> f64,
    ) -> f64 {
        let mean = centre(members, values, lengths, self.scale);
        let (mut statistic, mut terms) = (0.0, 0usize);
        for contig in members {
            for (sample, ((value, centre), spread)) in scaled(*contig, values, self.scale)
                .zip(&mean)
                .zip(self.spread(lengths[*contig]))
                .enumerate()
            {
                if *spread > 0.0 {
                    statistic += (value - centre).powi(2) / (spread + noise(*contig, sample));
                    terms += 1;
                }
            }
        }
        let used = mean.len().min(terms);
        match terms <= used {
            true => 1.0,
            false => chi_square_above(statistic, terms - used),
        }
    }
}
