use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Result, bail};
use rayon::prelude::*;

use crate::external::hmmer_engine;
use crate::markers::strain::masked;

const GAP: u8 = b'-';

pub struct Profile {
    pub name: String,
    pub length: usize,
    pub text: String,
}

impl Profile {
    pub fn read(path: &Path) -> Result<Vec<Self>> {
        let mut profiles = Vec::new();
        let reader = BufReader::new(std::fs::File::open(path)?);
        hmmer_engine::each_model(reader, |name, text| {
            let length = text
                .lines()
                .find_map(|line| line.strip_prefix("LENG"))
                .and_then(|rest| rest.trim().parse().ok());
            let Some(length) = length else {
                bail!("{name} in {} carries no length", path.display());
            };
            let name = name.to_string();
            profiles.push(Self { name, length, text });
            Ok(())
        })?;
        profiles.sort_by(|a, b| a.name.cmp(&b.name));
        Ok(profiles)
    }
}

pub struct References {
    names: Vec<String>,
    width: usize,
    rows: Vec<u8>,
}

fn residue(byte: u8) -> bool {
    !matches!(byte, GAP | b'.' | b'X' | b'x')
}

impl References {
    pub fn read(path: &Path) -> Result<Self> {
        let mut references = Self {
            names: Vec::new(),
            width: 0,
            rows: Vec::new(),
        };
        let reader = BufReader::new(std::fs::File::open(path)?);
        for line in reader.lines() {
            let line = line?;
            match line.strip_prefix('>') {
                Some(name) => references.names.push(name.trim().to_string()),
                None => references.rows.extend(line.trim().bytes()),
            }
        }
        let count = references.names.len();
        if count == 0 || !references.rows.len().is_multiple_of(count) {
            bail!("{} holds rows of unequal width", path.display());
        }
        references.width = references.rows.len() / count;
        Ok(references)
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn nearest(&self, row: &[u8]) -> Option<&str> {
        let columns = (0..row.len().min(self.width))
            .filter(|at| residue(row[*at]))
            .collect::<Vec<_>>();
        let tally = |reference: usize| {
            let held = &self.rows[reference * self.width..(reference + 1) * self.width];
            let (mut same, mut compared) = (0u64, 0u64);
            for at in &columns {
                if residue(held[*at]) {
                    compared += 1;
                    same += u64::from(held[*at] == row[*at]);
                }
            }
            (reference, same, compared)
        };
        let closer = |a: &(usize, u64, u64), b: &(usize, u64, u64)| {
            (a.1 * b.2)
                .cmp(&(b.1 * a.2))
                .then(a.2.cmp(&b.2))
                .then(b.0.cmp(&a.0))
        };
        let best = crate::pool::install(|| {
            (0..self.names.len())
                .into_par_iter()
                .map(tally)
                .filter(|(_, _, compared)| *compared > 0)
                .max_by(closer)
        })?;
        Some(self.names[best.0].as_str())
    }
}

pub fn concatenate(
    profiles: &[Profile],
    picked: &[Vec<(usize, String)>],
    slots: usize,
) -> Result<Vec<Option<Vec<u8>>>> {
    let width = profiles.iter().map(|profile| profile.length).sum();
    let aligned = crate::pool::install(|| {
        profiles
            .par_iter()
            .zip(picked)
            .map(|(profile, held)| {
                if held.is_empty() {
                    return Ok(Vec::new());
                }
                let proteins = held
                    .iter()
                    .map(|(_, protein)| protein.as_str())
                    .collect::<Vec<_>>();
                Ok(masked(
                    &hmmer_engine::align(&profile.text, &proteins)?,
                    held.len(),
                ))
            })
            .collect::<Result<Vec<_>>>()
    })?;
    let mut rows = vec![None; slots];
    let mut offset = 0;
    for ((profile, held), aligned) in profiles.iter().zip(picked).zip(aligned) {
        for ((slot, _), row) in held.iter().zip(aligned) {
            if row.len() != profile.length {
                continue;
            }
            let target = rows[*slot].get_or_insert_with(|| vec![GAP; width]);
            target[offset..offset + profile.length].copy_from_slice(&row);
        }
        offset += profile.length;
    }
    Ok(rows)
}
