use ndarray::{Array1, ArrayView2, Zip};

use crate::markers::replicon::DEPTH_FLOOR;
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

// Two sides of one genome sit apart by no more than two of its contigs do, so each signal is
// judged on the spread near complete bins show, never on how well each side's centre is known.
pub(crate) struct Agreement<'a> {
    depths: ArrayView2<'a, f64>,
    composition: ArrayView2<'a, f64>,
    lengths: &'a [usize],
    depth_spread: f64,
    composition_spread: Array1<f64>,
    present: f64,
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
        };
        let mut depth_spreads = Vec::new();
        let mut composition_spreads = vec![Vec::new(); composition.ncols()];
        for members in hosts {
            let sums = agreement.sums(members.iter().copied());
            let (depth_centre, composition_centre) = sums.centre();
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
        }
        agreement.depth_spread = median(depth_spreads);
        agreement.composition_spread = composition_spreads.into_iter().map(median).collect();
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

    // A sample where neither side has reads says nothing about either, so it adds no term.
    pub(crate) fn agree(&self, one: &Sums, other: &Sums) -> bool {
        if one.weight <= 0.0 || other.weight <= 0.0 {
            return false;
        }
        let (one_depth, one_composition) = one.centre();
        let (other_depth, other_composition) = other.centre();
        let (mut statistic, mut terms) = (0.0, 0usize);
        for (mine, theirs) in one_depth.iter().zip(&other_depth) {
            if mine.max(*theirs) > self.present {
                statistic += (mine - theirs).powi(2) / (2.0 * self.depth_spread);
                terms += 1;
            }
        }
        if terms == 0 || chi_square_above(statistic, terms) < FITS {
            return false;
        }
        let (mut statistic, mut terms) = (0.0, 0usize);
        for ((mine, theirs), spread) in one_composition
            .iter()
            .zip(&other_composition)
            .zip(&self.composition_spread)
        {
            if *spread > 0.0 {
                statistic += (mine - theirs).powi(2) / (2.0 * spread);
                terms += 1;
            }
        }
        terms > 0 && chi_square_above(statistic, terms) >= FITS
    }
}
