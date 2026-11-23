use ndarray::{Array1, ArrayView2, Zip};

use crate::markers::replicon::DEPTH_FLOOR;
use crate::recover::coherence::{LEAST, octave};
use crate::recover::fold::{FITS, chi_square_above, median};

#[derive(Clone)]
pub(crate) struct Sums {
    depth: Array1<f64>,
    composition: Array1<f64>,
    weight: f64,
}

impl Sums {
    pub(crate) fn without(&self, part: &Sums) -> Sums {
        Sums {
            depth: &self.depth - &part.depth,
            composition: &self.composition - &part.composition,
            weight: self.weight - part.weight,
        }
    }

    fn centre(&self) -> (Array1<f64>, Array1<f64>) {
        (&self.depth / self.weight, &self.composition / self.weight)
    }
}

pub(crate) struct Statistics {
    pub(crate) depth: (f64, usize),
    pub(crate) composition: (f64, usize),
}

// Two sides of one genome sit apart by no more than two of its contigs do, so each signal is
// judged on the spread near complete bins show, never on how well each side's centre is known.
pub(crate) struct Agreement<'a> {
    depths: ArrayView2<'a, f64>,
    composition: ArrayView2<'a, f64>,
    lengths: &'a [usize],
    depth_spread: f64,
    composition_spread: Array1<f64>,
    present: f64,
    octave_depth: Vec<Option<f64>>,
    octave_composition: Vec<Option<Array1<f64>>>,
    pub(crate) bars: (f64, f64),
}

fn nearest<T>(spreads: &[Option<T>], length: usize) -> Option<&T> {
    let at = octave(length) as isize;
    spreads
        .iter()
        .enumerate()
        .filter_map(|(octave, spread)| spread.as_ref().map(|spread| (octave as isize, spread)))
        .min_by_key(|(octave, _)| ((octave - at).abs(), *octave))
        .map(|(_, spread)| spread)
}

fn quantile(mut values: Vec<f64>, level: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_unstable_by(f64::total_cmp);
    let at = (level * (values.len() - 1) as f64).round() as usize;
    values[at]
}

impl<'a> Agreement<'a> {
    pub(crate) fn of(
        hosts: &[Vec<usize>],
        depths: ArrayView2<'a, f64>,
        composition: ArrayView2<'a, f64>,
        lengths: &'a [usize],
    ) -> Option<Self> {
        let present = (2.0 * DEPTH_FLOOR).ln();
        let mut agreement = Self {
            depths,
            composition,
            lengths,
            depth_spread: 0.0,
            composition_spread: Array1::zeros(composition.ncols()),
            present,
            octave_depth: Vec::new(),
            octave_composition: Vec::new(),
            bars: (f64::INFINITY, f64::INFINITY),
        };
        let mut depth_spreads = Vec::new();
        let mut composition_spreads = vec![Vec::new(); composition.ncols()];
        let mut by_octave_depth: Vec<Vec<f64>> = Vec::new();
        let mut by_octave_composition: Vec<Vec<Array1<f64>>> = Vec::new();
        for members in hosts {
            let (depth_centre, composition_centre) =
                agreement.sums(members.iter().copied()).centre();
            let (depth_spread, composition_spread) =
                agreement.spreads(members, &depth_centre, &composition_centre);
            depth_spreads.extend(
                depth_centre
                    .iter()
                    .zip(&depth_spread)
                    .filter(|(centre, _)| **centre > present)
                    .map(|(_, spread)| *spread),
            );
            for (column, spread) in composition_spreads.iter_mut().zip(&composition_spread) {
                column.push(*spread);
            }
            for (at, (depth, composition)) in agreement
                .octave_sums(members, &depth_centre, &composition_centre)
                .into_iter()
                .enumerate()
            {
                if by_octave_depth.len() <= at {
                    by_octave_depth.resize(at + 1, Vec::new());
                    by_octave_composition.resize(at + 1, Vec::new());
                }
                by_octave_depth[at].extend(depth);
                by_octave_composition[at].extend(composition);
            }
        }
        agreement.depth_spread = median(depth_spreads);
        agreement.composition_spread = composition_spreads.into_iter().map(median).collect();
        agreement.octave_depth = by_octave_depth
            .into_iter()
            .map(|values| (values.len() >= LEAST).then(|| median(values)))
            .collect();
        agreement.octave_composition = by_octave_composition
            .into_iter()
            .map(|hosts| {
                (hosts.len() >= LEAST).then(|| {
                    (0..composition.ncols())
                        .map(|column| median(hosts.iter().map(|host| host[column]).collect()))
                        .collect()
                })
            })
            .collect();
        let scatters = hosts
            .iter()
            .map(|members| agreement.scatter(members))
            .collect::<Vec<_>>();
        agreement.bars = (
            quantile(scatters.iter().map(|(one, _)| *one).collect(), 1.0 - FITS),
            quantile(
                scatters.iter().map(|(_, other)| *other).collect(),
                1.0 - FITS,
            ),
        );
        (agreement.depth_spread > 0.0).then_some(agreement)
    }

    fn row_logs(&self, contig: usize) -> Array1<f64> {
        self.depths
            .row(contig)
            .mapv(|depth| (depth + DEPTH_FLOOR).ln())
    }

    pub(crate) fn sums(&self, members: impl IntoIterator<Item = usize>) -> Sums {
        let mut sums = Sums {
            depth: Array1::zeros(self.depths.ncols()),
            composition: Array1::zeros(self.composition.ncols()),
            weight: 0.0,
        };
        for contig in members {
            let weight = self.lengths[contig] as f64;
            sums.depth.scaled_add(weight, &self.row_logs(contig));
            sums.composition
                .scaled_add(weight, &self.composition.row(contig));
            sums.weight += weight;
        }
        sums
    }

    fn spreads(
        &self,
        members: &[usize],
        depth_centre: &Array1<f64>,
        composition_centre: &Array1<f64>,
    ) -> (Array1<f64>, Array1<f64>) {
        let mut depth = Array1::zeros(depth_centre.len());
        let mut composition = Array1::zeros(composition_centre.len());
        let total = members.iter().map(|at| self.lengths[*at]).sum::<usize>() as f64;
        for contig in members {
            let weight = self.lengths[*contig] as f64 / total;
            Zip::from(&mut depth)
                .and(&self.row_logs(*contig))
                .and(depth_centre)
                .for_each(|sum, log, centre| *sum += weight * (log - centre).powi(2));
            Zip::from(&mut composition)
                .and(self.composition.row(*contig))
                .and(composition_centre)
                .for_each(|sum, value, centre| *sum += weight * (value - centre).powi(2));
        }
        (depth, composition)
    }

    // A short contig scatters more than a long one, so each length octave holds its own spread.
    fn octave_sums(
        &self,
        members: &[usize],
        depth_centre: &Array1<f64>,
        composition_centre: &Array1<f64>,
    ) -> Vec<(Vec<f64>, Option<Array1<f64>>)> {
        let mut sums: Vec<(Array1<f64>, Array1<f64>, usize)> = Vec::new();
        for contig in members {
            let at = octave(self.lengths[*contig]);
            if sums.len() <= at {
                sums.resize(
                    at + 1,
                    (
                        Array1::zeros(depth_centre.len()),
                        Array1::zeros(composition_centre.len()),
                        0,
                    ),
                );
            }
            let (depth, composition, count) = &mut sums[at];
            Zip::from(depth)
                .and(&self.row_logs(*contig))
                .and(depth_centre)
                .for_each(|sum, log, centre| *sum += (log - centre).powi(2));
            Zip::from(composition)
                .and(self.composition.row(*contig))
                .and(composition_centre)
                .for_each(|sum, value, centre| *sum += (value - centre).powi(2));
            *count += 1;
        }
        sums.into_iter()
            .map(|(depth, composition, count)| {
                if count < LEAST {
                    return (Vec::new(), None);
                }
                let depth = depth
                    .iter()
                    .zip(depth_centre)
                    .filter(|(_, centre)| **centre > self.present)
                    .map(|(sum, _)| sum / count as f64)
                    .collect();
                (depth, Some(composition / count as f64))
            })
            .collect()
    }

    // The mean of each contig's squared distance from the centre over its octave's spread, so a
    // bin of two genomes reads wider than a genome and a large bin reads no stricter than a small.
    pub(crate) fn scatter(&self, members: &[usize]) -> (f64, f64) {
        let (depth_centre, composition_centre) = self.sums(members.iter().copied()).centre();
        let (mut composition, mut depth, mut total) = (0.0, 0.0, 0.0);
        for contig in members {
            let length = self.lengths[*contig];
            let weight = length as f64;
            if let Some(spread) = nearest(&self.octave_composition, length) {
                let (mut sum, mut terms) = (0.0, 0usize);
                for ((value, centre), spread) in self
                    .composition
                    .row(*contig)
                    .iter()
                    .zip(&composition_centre)
                    .zip(spread)
                {
                    if *spread > 0.0 {
                        sum += (value - centre).powi(2) / spread;
                        terms += 1;
                    }
                }
                composition += weight * sum / terms.max(1) as f64;
            }
            if let Some(spread) =
                nearest(&self.octave_depth, length).filter(|spread| **spread > 0.0)
            {
                let (mut sum, mut terms) = (0.0, 0usize);
                for (log, centre) in self.row_logs(*contig).iter().zip(&depth_centre) {
                    if *centre > self.present {
                        sum += (log - centre).powi(2) / spread;
                        terms += 1;
                    }
                }
                depth += weight * sum / terms.max(1) as f64;
            }
            total += weight;
        }
        (composition / total, depth / total)
    }

    pub(crate) fn holds_one(&self, members: &[usize]) -> bool {
        let (composition, depth) = self.scatter(members);
        composition <= self.bars.0 && depth <= self.bars.1
    }

    // A sample where neither side has reads says nothing about either, so it adds no term.
    pub(crate) fn statistics(&self, one: &Sums, other: &Sums) -> Statistics {
        let (one_depth, one_composition) = one.centre();
        let (other_depth, other_composition) = other.centre();
        let (mut depth, mut depth_terms) = (0.0, 0usize);
        for (mine, theirs) in one_depth.iter().zip(&other_depth) {
            if mine.max(*theirs) > self.present {
                depth += (mine - theirs).powi(2) / (2.0 * self.depth_spread);
                depth_terms += 1;
            }
        }
        let (mut composition, mut composition_terms) = (0.0, 0usize);
        for ((mine, theirs), spread) in one_composition
            .iter()
            .zip(&other_composition)
            .zip(&self.composition_spread)
        {
            if *spread > 0.0 {
                composition += (mine - theirs).powi(2) / (2.0 * spread);
                composition_terms += 1;
            }
        }
        Statistics {
            depth: (depth, depth_terms),
            composition: (composition, composition_terms),
        }
    }

    pub(crate) fn agree(&self, one: &Sums, other: &Sums) -> bool {
        if one.weight <= 0.0 || other.weight <= 0.0 {
            return false;
        }
        let found = self.statistics(one, other);
        let fits = |(statistic, terms): (f64, usize)| {
            terms > 0 && chi_square_above(statistic, terms) >= FITS
        };
        fits(found.depth) && fits(found.composition)
    }
}
