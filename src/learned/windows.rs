use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use rayon::prelude::*;

use crate::coverage::depth_file::{self, Contigs};
use crate::learned::pieces::{view_key, views};

pub const WINDOW: usize = 100;

pub struct Windows {
    dir: PathBuf,
    samples: Vec<String>,
    contigs: Contigs,
    first: HashMap<String, usize>,
}

// Each contig's spans in row order, so a view reads the depth of its own span.
pub struct SpanDepth {
    samples: usize,
    starts: Vec<usize>,
    values: Vec<f32>,
}

impl SpanDepth {
    pub fn view(&self, row: usize, view: usize) -> &[f32] {
        let at = (self.starts[row] + view) * self.samples;
        &self.values[at..at + self.samples]
    }
}

impl Windows {
    pub fn read(dir: &Path, samples: &[String]) -> Result<Self> {
        let contigs = Contigs::read(dir)?;
        let mut first = HashMap::new();
        let mut total = 0;
        for (name, length) in contigs.names.iter().zip(&contigs.lengths) {
            first.insert(name.clone(), total);
            total += length.div_ceil(WINDOW);
        }
        if let Some(sample) = samples
            .iter()
            .find(|sample| !depth_file::path(dir, sample).is_file())
        {
            bail!("{} holds no depth for {sample}", dir.display());
        }
        Ok(Self {
            dir: dir.to_path_buf(),
            samples: samples.to_vec(),
            contigs,
            first,
        })
    }

    pub fn spans(&self, names: &[String], lengths: &[usize], seed: u64) -> Result<SpanDepth> {
        self.cut(names, |at| views(view_key(seed, &names[at]), lengths[at]))
    }

    pub fn cut(
        &self,
        names: &[String],
        cut: impl Fn(usize) -> Vec<(usize, usize)>,
    ) -> Result<SpanDepth> {
        let requests = names
            .iter()
            .enumerate()
            .map(|(at, name)| {
                let window = *self
                    .first
                    .get(name)
                    .with_context(|| format!("{name} has no windows"))?;
                Ok((window, cut(at)))
            })
            .collect::<Result<Vec<_>>>()?;
        let starts = requests
            .iter()
            .scan(0, |total, (_, spans)| {
                let start = *total;
                *total += spans.len();
                Some(start)
            })
            .collect::<Vec<_>>();
        let spans = requests.iter().map(|(_, spans)| spans.len()).sum::<usize>();
        let per_sample = self
            .samples
            .par_iter()
            .map(|sample| {
                let bytes =
                    depth_file::windows(&depth_file::path(&self.dir, sample), &self.contigs)?;
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
        let mut values = vec![0.0; spans * samples];
        for (sample, depths) in per_sample.iter().enumerate() {
            for (at, depth) in depths.iter().enumerate() {
                values[at * samples + sample] = *depth;
            }
        }
        Ok(SpanDepth {
            samples,
            starts,
            values,
        })
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

// The inverse of f16_value, rounding half up and flushing what half precision cannot hold as a
// normal number to zero, as the window files were first written.
pub fn f16_bits(value: f32) -> u16 {
    let bits = value.to_bits();
    let exponent = ((bits >> 23) & 0xff) as i32 - 112;
    if exponent <= 0 {
        return 0;
    }
    let mut rounded = (bits & 0x7f_ffff) + 0x1000;
    let mut exponent = exponent as u32;
    if rounded & 0x80_0000 != 0 {
        rounded = 0;
        exponent += 1;
    }
    if exponent >= 31 {
        return 0x7bff;
    }
    ((exponent << 10) | (rounded >> 13)) as u16
}
