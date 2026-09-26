use std::collections::HashSet;

use anyhow::Result;
use ndarray::{Array, Array2, Axis};

/// The coverage and composition tables are row-aligned with the contig name list, so a contig
/// dropped from one has to leave all of them at the same index.
pub fn keep_rows(table: &Array2<f64>, dropped: &HashSet<usize>) -> Result<Array2<f64>> {
    let kept = table
        .axis_iter(Axis(0))
        .enumerate()
        .filter(|(index, _)| !dropped.contains(index))
        .flat_map(|(_, row)| row.to_vec());
    let rows = table.nrows() - dropped.len();
    Ok(Array::from_iter(kept).into_shape_with_order((rows, table.ncols()))?)
}

pub fn keep<T: Clone>(values: &[T], dropped: &HashSet<usize>) -> Vec<T> {
    values
        .iter()
        .enumerate()
        .filter(|(index, _)| !dropped.contains(index))
        .map(|(_, value)| value.clone())
        .collect()
}

pub fn reorder<T: Clone>(values: &[T], order: &[usize]) -> Vec<T> {
    order.iter().map(|at| values[*at].clone()).collect()
}

pub fn dropped_names(names: &[String], dropped: &HashSet<usize>) -> HashSet<String> {
    names
        .iter()
        .enumerate()
        .filter(|(index, _)| dropped.contains(index))
        .map(|(_, name)| name.clone())
        .collect()
}
