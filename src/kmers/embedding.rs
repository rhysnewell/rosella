use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};
use ndarray::Array2;

pub struct Embedding {
    width: usize,
    whole: HashMap<String, usize>,
    rows: Vec<f32>,
    halves: HashMap<String, usize>,
    first: Vec<f32>,
    second: Vec<f32>,
}

fn floats(path: &Path) -> Result<Vec<f32>> {
    let bytes = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
    if bytes.len() % 4 != 0 {
        bail!("{} is not a whole number of floats", path.display());
    }
    Ok(bytes
        .chunks_exact(4)
        .map(|chunk| f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
        .collect())
}

fn names(path: &Path) -> Result<HashMap<String, usize>> {
    Ok(std::fs::read_to_string(path)
        .with_context(|| format!("reading {}", path.display()))?
        .lines()
        .enumerate()
        .map(|(row, name)| (name.to_string(), row))
        .collect())
}

fn gather<'a>(
    index: &HashMap<String, usize>,
    values: &[f32],
    width: usize,
    wanted: impl Iterator<Item = &'a str>,
) -> Result<Array2<f64>> {
    let mut out = Vec::new();
    for name in wanted {
        let row = index
            .get(name)
            .ok_or_else(|| anyhow!("{name} has no composition embedding"))?;
        out.extend(
            values[row * width..(row + 1) * width]
                .iter()
                .map(|value| *value as f64),
        );
    }
    Ok(Array2::from_shape_vec((out.len() / width, width), out)?)
}

impl Embedding {
    pub fn read(directory: &Path) -> Result<Self> {
        let whole = names(&directory.join("whole.names"))?;
        let rows = floats(&directory.join("whole.f32"))?;
        let halves = names(&directory.join("halves.names"))?;
        let first = floats(&directory.join("first.f32"))?;
        let second = floats(&directory.join("second.f32"))?;
        let width = rows.len() / whole.len().max(1);
        if width == 0
            || rows.len() != width * whole.len()
            || first.len() != width * halves.len()
            || second.len() != first.len()
        {
            bail!(
                "the embedding files in {} disagree on their width",
                directory.display()
            );
        }
        Ok(Self {
            width,
            whole,
            rows,
            halves,
            first,
            second,
        })
    }

    pub fn rows<'a>(&self, wanted: impl Iterator<Item = &'a str>) -> Result<Array2<f64>> {
        gather(&self.whole, &self.rows, self.width, wanted)
    }

    pub fn halves(&self, wanted: &[&str]) -> Result<[Array2<f64>; 2]> {
        Ok([
            gather(
                &self.halves,
                &self.first,
                self.width,
                wanted.iter().copied(),
            )?,
            gather(
                &self.halves,
                &self.second,
                self.width,
                wanted.iter().copied(),
            )?,
        ])
    }
}
