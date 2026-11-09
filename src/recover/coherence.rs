use ndarray::{ArrayView2, Axis};

use crate::markers::replicon::DEPTH_FLOOR;
use crate::recover::fold::{FITS, chi_square_above, median};

const FIRST_OCTAVE: usize = 1_000;
const LEAST: usize = 3;

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

fn centre(members: &[usize], depths: ArrayView2<f64>, lengths: &[usize]) -> Vec<f64> {
    let mut mean = vec![0.0; depths.ncols()];
    let total = members.iter().map(|contig| lengths[*contig]).sum::<usize>() as f64;
    for contig in members {
        let weight = lengths[*contig] as f64 / total;
        for (sum, log) in mean.iter_mut().zip(logs(*contig, depths)) {
            *sum += weight * log;
        }
    }
    mean
}

impl Coherence {
    pub fn of(hosts: &[&[usize]], depths: ArrayView2<f64>, lengths: &[usize]) -> Option<Self> {
        let samples = depths.ncols();
        let mut held: Vec<Vec<Vec<f64>>> = Vec::new();
        for members in hosts {
            let mean = centre(members, depths, lengths);
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

    // With few samples the test has little power and keeps the bin, so a whole small genome
    // the depth cannot speak to stays binned.
    pub fn holds_one(&self, members: &[usize], depths: ArrayView2<f64>, lengths: &[usize]) -> bool {
        let mean = centre(members, depths, lengths);
        let (mut statistic, mut terms) = (0.0, 0usize);
        for contig in members {
            for ((log, centre), spread) in logs(*contig, depths)
                .zip(&mean)
                .zip(self.spread(lengths[*contig]))
            {
                if *spread > 0.0 {
                    statistic += (log - centre).powi(2) / spread;
                    terms += 1;
                }
            }
        }
        let used = mean.len().min(terms);
        terms <= used || chi_square_above(statistic, terms - used) >= FITS
    }
}
