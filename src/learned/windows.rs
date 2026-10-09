use std::collections::HashMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use rayon::prelude::*;

use crate::learned::pieces::{VIEWS, view_key, views};

const WINDOW: usize = 100;

pub struct Windows {
    dir: PathBuf,
    samples: Vec<String>,
    first: HashMap<String, usize>,
}

// Each partitioned contig's views in row order, the whole first, so a view reads the depth of
// its own span.
pub struct SpanDepth {
    samples: usize,
    values: Vec<f32>,
}

impl SpanDepth {
    pub fn view(&self, row: usize, view: usize) -> &[f32] {
        let at = (row * VIEWS + view) * self.samples;
        &self.values[at..at + self.samples]
    }
}

impl Windows {
    pub fn read(dir: &Path, samples: &[String]) -> Result<Self> {
        let mut first = HashMap::new();
        let mut total = 0;
        for line in crate::get_file_reader(dir.join("contigs.tsv"))?
            .lines()
            .skip(1)
        {
            let line = line?;
            let (name, length) = line
                .split_once('\t')
                .context("a contig row has no length")?;
            first.insert(name.to_string(), total);
            total += length.parse::<usize>()?.div_ceil(WINDOW);
        }
        for sample in samples {
            let size = std::fs::metadata(dir.join(format!("{sample}.win.f16")))?.len() as usize;
            if size != 2 * total {
                bail!(
                    "{sample}.win.f16 holds {size} bytes where {total} windows need {}",
                    2 * total
                );
            }
        }
        Ok(Self {
            dir: dir.to_path_buf(),
            samples: samples.to_vec(),
            first,
        })
    }

    pub fn spans(
        &self,
        names: &[String],
        lengths: &[usize],
        long: usize,
        seed: u64,
    ) -> Result<SpanDepth> {
        let requests = names
            .iter()
            .zip(lengths)
            .take(long)
            .map(|(name, length)| {
                let window = *self
                    .first
                    .get(name)
                    .with_context(|| format!("{name} has no windows"))?;
                Ok((window, views(view_key(seed, name), *length)))
            })
            .collect::<Result<Vec<_>>>()?;
        let per_sample = self
            .samples
            .par_iter()
            .map(|sample| {
                let bytes = std::fs::read(self.dir.join(format!("{sample}.win.f16")))?;
                let window = |at: usize| {
                    f64::from(f16_value(u16::from_le_bytes([
                        bytes[2 * at],
                        bytes[2 * at + 1],
                    ])))
                };
                Ok(requests
                    .iter()
                    .flat_map(|&(first, ref spans)| {
                        spans.iter().map(move |&(start, end)| {
                            let sum = (start / WINDOW..end.div_ceil(WINDOW))
                                .map(|at| {
                                    let overlap =
                                        end.min((at + 1) * WINDOW) - start.max(at * WINDOW);
                                    overlap as f64 * window(first + at)
                                })
                                .sum::<f64>();
                            (sum / (end - start).max(1) as f64) as f32
                        })
                    })
                    .collect::<Vec<_>>())
            })
            .collect::<Result<Vec<_>>>()?;
        let samples = self.samples.len();
        let mut values = vec![0.0; requests.len() * VIEWS * samples];
        for (sample, depths) in per_sample.iter().enumerate() {
            for (at, depth) in depths.iter().enumerate() {
                values[at * samples + sample] = *depth;
            }
        }
        Ok(SpanDepth { samples, values })
    }
}

fn f16_value(bits: u16) -> f32 {
    let exponent = u32::from((bits >> 10) & 0x1f);
    let mantissa = u32::from(bits & 0x3ff);
    match exponent {
        0 => mantissa as f32 * (2.0f32).powi(-24),
        _ => f32::from_bits(((exponent + 112) << 23) | (mantissa << 13)),
    }
}
