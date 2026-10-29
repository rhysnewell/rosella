use std::cmp::Ordering;
use std::collections::HashMap;
use std::path::Path;

use anyhow::Result;
use log::info;
use rand::{Rng, rngs::StdRng};
use rayon::prelude::*;

use crate::coverage::depth_file::{self, Contigs};

// Near-clone strains share the sequence a read ties across. A genome's own contigs share none.
// So tied contigs train in one step as each other's negatives.
pub struct Partners {
    starts: Vec<usize>,
    rows: Vec<u32>,
    weights: Vec<f64>,
}

impl Partners {
    pub fn read(dir: &Path, samples: &[String], names: &[String]) -> Result<Self> {
        let at = names
            .iter()
            .enumerate()
            .map(|(row, name)| (name.as_str(), row as u32))
            .collect::<HashMap<_, _>>();
        let row_of = Contigs::read(dir)?
            .names
            .iter()
            .map(|name| at.get(name.as_str()).copied())
            .collect::<Vec<_>>();
        let pairs = samples
            .par_iter()
            .map(|sample| {
                let mut pairs = depth_file::triples(&depth_file::ties_path(dir, sample))?
                    .into_iter()
                    .filter_map(|[first, second, count]| {
                        let first = row_of.get(first as usize).copied().flatten()?;
                        let second = row_of.get(second as usize).copied().flatten()?;
                        Some((key(first, second), f64::from(count)))
                    })
                    .collect::<Vec<_>>();
                pairs.sort_unstable_by_key(|(key, _)| *key);
                Ok::<_, anyhow::Error>(merged(pairs, Vec::new()))
            })
            .try_reduce(Vec::new, |left, right| Ok(merged(left, right)))?;
        let partners = Self::of(names.len(), &pairs);
        info!(
            "Tie partners: {} pairs over {} of {} training contigs",
            pairs.len(),
            partners
                .starts
                .windows(2)
                .filter(|span| span[1] > span[0])
                .count(),
            names.len()
        );
        Ok(partners)
    }

    fn of(contigs: usize, pairs: &[(u64, f64)]) -> Self {
        let mut starts = vec![0; contigs + 1];
        for (key, _) in pairs {
            let (first, second) = split(*key);
            starts[first as usize + 1] += 1;
            starts[second as usize + 1] += 1;
        }
        for at in 0..contigs {
            starts[at + 1] += starts[at];
        }
        let mut next = starts.clone();
        let (mut rows, mut weights) = (vec![0; 2 * pairs.len()], vec![0.0; 2 * pairs.len()]);
        for (key, weight) in pairs {
            let (first, second) = split(*key);
            for (from, to) in [(first, second), (second, first)] {
                let slot = &mut next[from as usize];
                (rows[*slot], weights[*slot]) = (to, *weight);
                *slot += 1;
            }
        }
        Self {
            starts,
            rows,
            weights,
        }
    }

    pub fn paired(&self, order: &[u32], rng: &mut StdRng) -> Vec<u32> {
        let mut open = vec![false; self.starts.len() - 1];
        order.iter().for_each(|row| open[*row as usize] = true);
        let mut out = Vec::with_capacity(order.len());
        for row in order {
            if !std::mem::take(&mut open[*row as usize]) {
                continue;
            }
            out.push(*row);
            if let Some(partner) = self.draw(*row, &open, rng) {
                open[partner as usize] = false;
                out.push(partner);
            }
        }
        out
    }

    fn draw(&self, row: u32, open: &[bool], rng: &mut StdRng) -> Option<u32> {
        let span = self.starts[row as usize]..self.starts[row as usize + 1];
        let free = || span.clone().filter(|at| open[self.rows[*at] as usize]);
        let total = free().map(|at| self.weights[at]).sum::<f64>();
        if total <= 0.0 {
            return None;
        }
        let mut left = rng.random::<f64>() * total;
        let mut last = None;
        for at in free() {
            last = Some(self.rows[at]);
            left -= self.weights[at];
            if left <= 0.0 {
                break;
            }
        }
        last
    }
}

fn key(first: u32, second: u32) -> u64 {
    (u64::from(first.min(second)) << 32) | u64::from(first.max(second))
}

fn split(key: u64) -> (u32, u32) {
    ((key >> 32) as u32, key as u32)
}

fn merged(left: Vec<(u64, f64)>, right: Vec<(u64, f64)>) -> Vec<(u64, f64)> {
    let mut out = Vec::with_capacity(left.len().max(right.len()));
    let (mut left, mut right) = (left.into_iter().peekable(), right.into_iter().peekable());
    while let (Some(&(first, _)), Some(&(second, _))) = (left.peek(), right.peek()) {
        out.extend(match first.cmp(&second) {
            Ordering::Less => left.next(),
            Ordering::Greater => right.next(),
            Ordering::Equal => left.next().zip(right.next()).map(|(a, b)| (a.0, a.1 + b.1)),
        });
    }
    out.extend(left);
    out.extend(right);
    out
}
