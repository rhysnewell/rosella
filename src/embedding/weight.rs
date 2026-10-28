use rand::{SeedableRng, rngs::StdRng};
use rayon::prelude::*;

use crate::embedding::metrics::{
    Abundance, Centred, DistanceSettings, abundance_distance, combine, rho_between,
};

pub mod noise;

use noise::Partners;

pub const NEIGHBOURS: usize = 10;
pub const STEPS: usize = 20;
const STEP: f64 = 0.05;

pub struct Contigs<'a> {
    pub coverage: Vec<&'a [f64]>,
    pub whole: Vec<&'a [f64]>,
    pub first: Vec<&'a [f64]>,
    pub second: Vec<&'a [f64]>,
}

// Each contig's first half looks for its second half among every contig's second half. The weight
// that finds the most within NEIGHBOURS is the one that separates genomes on this assembly.
pub fn recall(contigs: &Contigs, settings: &DistanceSettings, seed: u64) -> Option<[f64; STEPS]> {
    let n = contigs.coverage.len();
    if n <= NEIGHBOURS {
        return None;
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let partners = noise::partners(&contigs.coverage, &contigs.whole, NEIGHBOURS, &mut rng)?;
    let halves = |rows: &[&[f64]]| -> Vec<Centred> {
        rows.par_iter().map(|row| Centred::new(row)).collect()
    };
    let prepared = Prepared {
        abundance: contigs
            .coverage
            .par_iter()
            .map(|row| Abundance::new(row, settings))
            .collect(),
        first: halves(&contigs.first),
        second: halves(&contigs.second),
        settings,
    };
    let found = (0..n)
        .into_par_iter()
        .map(|contig| recalled(&prepared, &partners[contig], contig))
        .collect::<Vec<_>>();
    Some(std::array::from_fn(|step| {
        found.iter().map(|held| held[step]).sum::<f64>() / n as f64
    }))
}

// A cut below the best makes the weight a step function of the pool. On marine 2 per cent more
// contigs moved it from 0.575 to 0.65. The peak of a least squares parabola moves smoothly.
pub fn centre(recall: &[f64; STEPS]) -> f64 {
    let [a, b, _] = parabola(recall);
    if a < 0.0 && b.is_finite() {
        return (-b / (2.0 * a)).clamp(0.0, weight_at(STEPS - 1));
    }
    (0..STEPS)
        .max_by(|left, right| recall[*left].total_cmp(&recall[*right]))
        .map_or(0.0, weight_at)
}

fn parabola(recall: &[f64; STEPS]) -> [f64; 3] {
    let mut moments = [0.0; 5];
    let mut against = [0.0; 3];
    for (step, value) in recall.iter().enumerate() {
        let mut power = 1.0;
        for (at, moment) in moments.iter_mut().enumerate() {
            *moment += power;
            if let Some(sum) = against.get_mut(at) {
                *sum += value * power;
            }
            power *= weight_at(step);
        }
    }
    let grid: [[f64; 3]; 3] =
        std::array::from_fn(|row| std::array::from_fn(|column| moments[4 - row - column]));
    let sums = [against[2], against[1], against[0]];
    let whole = determinant(&grid);
    std::array::from_fn(|column| {
        let mut swapped = grid;
        for (row, sum) in swapped.iter_mut().zip(sums) {
            row[column] = sum;
        }
        determinant(&swapped) / whole
    })
}

fn determinant(grid: &[[f64; 3]; 3]) -> f64 {
    grid[0][0] * (grid[1][1] * grid[2][2] - grid[1][2] * grid[2][1])
        - grid[0][1] * (grid[1][0] * grid[2][2] - grid[1][2] * grid[2][0])
        + grid[0][2] * (grid[1][0] * grid[2][1] - grid[1][1] * grid[2][0])
}

fn weight_at(step: usize) -> f64 {
    step as f64 * STEP
}

struct Prepared<'a> {
    abundance: Vec<Abundance>,
    first: Vec<Centred>,
    second: Vec<Centred>,
    settings: &'a DistanceSettings,
}

// The expectation over every candidate partner rather than one draw, which left the weight
// swinging twofold with the seed.
fn recalled(prepared: &Prepared, partners: &Partners, contig: usize) -> [f64; STEPS] {
    let mine = &prepared.abundance[contig];
    let own_composition = rho_between(&prepared.first[contig], &prepared.second[contig]);
    let own = partners
        .rows
        .iter()
        .map(|row| {
            let partner = Abundance::new(row, prepared.settings);
            let own_coverage = abundance_distance(mine, &partner).0;
            std::array::from_fn::<f64, STEPS, _>(|step| {
                combine(own_coverage, own_composition, weight_at(step))
            })
        })
        .collect::<Vec<_>>();
    // A partner is recalled when fewer than NEIGHBOURS others sit nearer, which is whether the
    // NEIGHBOURS-th nearest other sits no nearer than it, so one selection answers every partner.
    let others = prepared.abundance.len() - 1;
    let mut columns: [Vec<f64>; STEPS] = std::array::from_fn(|_| Vec::with_capacity(others));
    for other in (0..prepared.abundance.len()).filter(|other| *other != contig) {
        let coverage = abundance_distance(mine, &prepared.abundance[other]).0;
        let composition = rho_between(&prepared.first[contig], &prepared.second[other]);
        for (step, column) in columns.iter_mut().enumerate() {
            column.push(combine(coverage, composition, weight_at(step)));
        }
    }
    let nearest = columns
        .iter_mut()
        .map(|column| {
            *column
                .select_nth_unstable_by(NEIGHBOURS - 1, f64::total_cmp)
                .1
        })
        .collect::<Vec<_>>();
    std::array::from_fn(|step| {
        own.iter()
            .zip(&partners.weights)
            .map(|(bars, weight)| weight * f64::from(u8::from(nearest[step] >= bars[step])))
            .sum()
    })
}
