use anyhow::{Result, bail};
use log::debug;
use ndarray::{Array2, Axis};
use rayon::prelude::*;

const START: usize = 40;
const STEP: usize = 2;
const CAP: usize = 60;
pub const TARGET: f64 = 0.75;
const SWEEPS: usize = 60;
const TOLERANCE: f64 = 1e-12;

pub fn project(table: &mut Array2<f64>, target: f64) -> Result<()> {
    let n_rows = table.nrows();
    let n_cols = table.ncols();
    if n_rows < 2 || n_cols < 2 {
        bail!("a projection needs at least two contigs and two columns, got {n_rows}x{n_cols}");
    }

    let means = table.mean_axis(Axis(0)).expect("the table has rows");
    table
        .axis_iter_mut(Axis(0))
        .into_par_iter()
        .for_each(|mut row| row.iter_mut().zip(&means).for_each(|(v, m)| *v -= m));

    let covariance = covariance_of(table);
    let (values, vectors) = jacobi_eigen(&covariance);
    let kept = components_for(&values, target).min(n_cols);
    debug!(
        "Composition projected from {n_cols} to {kept} components at {:.3} explained variance",
        explained(&values, kept)
    );

    let basis = vectors.slice(ndarray::s![.., ..kept]).to_owned();
    *table = table.dot(&basis);
    Ok(())
}

fn covariance_of(table: &Array2<f64>) -> Array2<f64> {
    let n_cols = table.ncols();
    let scale = (table.nrows() - 1) as f64;
    let mut covariance = table.t().dot(table);
    covariance.iter_mut().for_each(|value| *value /= scale);
    debug_assert_eq!(covariance.dim(), (n_cols, n_cols));
    covariance
}

fn jacobi_eigen(input: &Array2<f64>) -> (Vec<f64>, Array2<f64>) {
    let n = input.nrows();
    let mut a = input.clone();
    let mut v = Array2::<f64>::eye(n);

    for _ in 0..SWEEPS {
        let off = (0..n)
            .flat_map(|p| ((p + 1)..n).map(move |q| (p, q)))
            .map(|(p, q)| a[[p, q]] * a[[p, q]])
            .sum::<f64>();
        if off <= TOLERANCE {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[[p, q]];
                if apq.abs() <= TOLERANCE {
                    continue;
                }
                let theta = (a[[q, q]] - a[[p, p]]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let akp = a[[k, p]];
                    let akq = a[[k, q]];
                    a[[k, p]] = c * akp - s * akq;
                    a[[k, q]] = s * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[[p, k]];
                    let aqk = a[[q, k]];
                    a[[p, k]] = c * apk - s * aqk;
                    a[[q, k]] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let vkp = v[[k, p]];
                    let vkq = v[[k, q]];
                    v[[k, p]] = c * vkp - s * vkq;
                    v[[k, q]] = s * vkp + c * vkq;
                }
            }
        }
    }

    let mut order = (0..n).collect::<Vec<_>>();
    order.sort_by(|left, right| a[[*right, *right]].total_cmp(&a[[*left, *left]]));
    let values = order.iter().map(|index| a[[*index, *index]]).collect();
    let mut sorted = Array2::<f64>::zeros((n, n));
    for (column, index) in order.into_iter().enumerate() {
        sorted.column_mut(column).assign(&v.column(index));
    }
    (values, sorted)
}

fn explained(values: &[f64], kept: usize) -> f64 {
    let total = values.iter().filter(|value| **value > 0.0).sum::<f64>();
    if total <= 0.0 {
        return 1.0;
    }
    values.iter().take(kept).sum::<f64>() / total
}

fn components_for(values: &[f64], target: f64) -> usize {
    let mut kept = START.min(values.len());
    while kept < CAP.min(values.len()) && explained(values, kept) < target {
        kept += STEP;
    }
    kept
}
