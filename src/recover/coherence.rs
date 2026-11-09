use ndarray::{ArrayView2, Axis};

use crate::markers::replicon::DEPTH_FLOOR;
use crate::recover::fold::{FITS, chi_square_above, median};

const FIRST_OCTAVE: usize = 1_000;
const LEAST: usize = 3;
const EVEN_GAP: f64 = 6.0;

// A short contig's depth scatters more than a long one's, so each length octave is judged on the
// spread near complete genomes show at that length.
pub struct Coherence {
    spreads: Vec<Option<Vec<f64>>>,
}

fn octave(length: usize) -> usize {
    (length.max(FIRST_OCTAVE) / FIRST_OCTAVE).ilog2() as usize
}

fn logs(contig: usize, depths: ArrayView2<'_, f64>) -> impl Iterator<Item = f64> + '_ {
    depths
        .index_axis_move(Axis(0), contig)
        .into_iter()
        .map(|depth| (depth + DEPTH_FLOOR).ln())
}

fn centre(
    members: &[usize],
    depths: ArrayView2<f64>,
    lengths: &[usize],
    share: impl Fn(usize, usize) -> f64,
) -> Vec<f64> {
    let weight = |contig: usize, sample: usize| lengths[contig] as f64 * share(contig, sample);
    let totals = (0..depths.ncols())
        .map(|sample| {
            members
                .iter()
                .map(|contig| weight(*contig, sample))
                .sum::<f64>()
        })
        .collect::<Vec<_>>();
    let mut mean = vec![0.0; depths.ncols()];
    for contig in members {
        for (sample, (sum, log)) in mean.iter_mut().zip(logs(*contig, depths)).enumerate() {
            *sum += weight(*contig, sample) / totals[sample] * log;
        }
    }
    mean
}

impl Coherence {
    pub fn of(hosts: &[&[usize]], depths: ArrayView2<f64>, lengths: &[usize]) -> Option<Self> {
        let samples = depths.ncols();
        let mut held: Vec<Vec<Vec<f64>>> = Vec::new();
        for members in hosts {
            let mean = centre(members, depths, lengths, |_, _| 1.0);
            let mut sums: Vec<(Vec<f64>, usize)> = Vec::new();
            for contig in *members {
                let at = octave(lengths[*contig]);
                if sums.len() <= at {
                    sums.resize(at + 1, (vec![0.0; samples], 0));
                }
                let (sum, count) = &mut sums[at];
                for ((total, log), centre) in sum.iter_mut().zip(logs(*contig, depths)).zip(&mean) {
                    *total += (log - centre).powi(2);
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
                        .map(|sample| median(bins.iter().map(|bin| bin[sample]).collect()))
                        .collect()
                })
            })
            .collect::<Vec<_>>();
        spreads
            .iter()
            .any(Option::is_some)
            .then_some(Self { spreads })
    }

    fn spread(&self, length: usize) -> &[f64] {
        let at = octave(length) as isize;
        self.spreads
            .iter()
            .enumerate()
            .filter_map(|(octave, spread)| spread.as_ref().map(|spread| (octave as isize, spread)))
            .min_by_key(|(octave, _)| ((octave - at).abs(), *octave))
            .map(|(_, spread)| spread.as_slice())
            .unwrap_or_default()
    }

    // A low depth contig holds few reads, so each term carries its count noise on top of the genome
    // spread. A contig whose depth is mostly tie split noise says little about where the bin sits.
    pub fn holds_one(
        &self,
        members: &[usize],
        depths: ArrayView2<f64>,
        lengths: &[usize],
        read_lengths: &[f64],
        tied: impl Fn(usize, usize) -> f64,
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
        let noise = |contig: usize, sample: usize| {
            let depth = expected[sample];
            let spread = self
                .spread(lengths[contig])
                .get(sample)
                .copied()
                .unwrap_or(0.0);
            let count = depth * read_lengths[sample] / lengths[contig] as f64;
            spread + count / (depth + DEPTH_FLOOR).powi(2)
        };
        let mean = centre(members, depths, lengths, |contig, sample| {
            match tied(contig, sample) {
                0.0 => 1.0,
                tie => noise(contig, sample) / (noise(contig, sample) + tie),
            }
        });
        let (mut statistic, mut terms) = (0.0, 0usize);
        for contig in members {
            for (sample, ((log, centre), spread)) in logs(*contig, depths)
                .zip(&mean)
                .zip(self.spread(lengths[*contig]))
                .enumerate()
            {
                if *spread > 0.0 {
                    let variance = noise(*contig, sample) + tied(*contig, sample);
                    statistic += (log - centre).powi(2) / variance;
                    terms += 1;
                }
            }
        }
        let used = mean.len().min(terms);
        terms <= used || chi_square_above(statistic, terms - used) >= FITS
    }
}

// The tie split and the truth can each sit anywhere from the depth a contig's own reads give to the
// depth every read it ties would give, so their gap carries two even draws' variance, width² / 6.
pub fn tie_spread(reads: u32, tied: u64, read_length: f64, length: usize) -> f64 {
    let depth = |count: f64| count * read_length / length as f64 + DEPTH_FLOOR;
    let own = (f64::from(reads) - tied as f64 / 2.0).max(0.0);
    (depth(own + tied as f64) / depth(own)).ln().powi(2) / EVEN_GAP
}
