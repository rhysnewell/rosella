use anyhow::Result;

use crate::learned::device::{Device, Ptr, Table, current, spread};

// Enough queries that their warps fill the card, and a product tile of 512 MB.
const QUERY_TILE: usize = 32768;
const BASE_TILE: usize = 4096;

pub const NOBODY: u32 = u32::MAX;

// A warp holds a query's list four slots to a lane.
const MOST_NEAREST: usize = 128;

pub struct Rows<'a> {
    pub values: &'a [f32],
    pub variance: &'a [f32],
}

// `rows` names each query's row of `moments` and `presence`, and a base row is its own.
pub struct Coverage<'a> {
    pub rows: &'a [u32],
    pub moments: &'a [f64],
    pub presence: &'a [f64],
    pub samples: usize,
    pub terms: &'a [f64],
}

#[derive(Default)]
struct Mixed {
    moments: Ptr,
    presence: Ptr,
    stride: i32,
    samples: i32,
    terms: Ptr,
}

// Every query's k nearest base rows, measured exactly. A query's own base row, named in `selves`,
// is never its neighbour.
pub fn nearest(
    base: Rows,
    queries: Rows,
    selves: &[u32],
    coverage: Option<Coverage>,
    k: usize,
) -> Result<Option<(Vec<u32>, Vec<f32>)>> {
    let Some(device) = current().filter(|_| k <= MOST_NEAREST) else {
        return Ok(None);
    };
    let mark = device.mark();
    let found = measure(&device, base, queries, selves, coverage, k);
    device.release(mark)?;
    found.map(Some)
}

fn reals(device: &Device, values: &[f64]) -> Result<Ptr> {
    let wide = device.wide();
    let at = device.alloc(values.len().max(1) * if wide { 8 } else { 4 })?;
    match wide {
        true => device.upload(at, values)?,
        false => device.upload(at, &values.iter().map(|v| *v as f32).collect::<Vec<_>>())?,
    }
    Ok(at)
}

fn measure(
    device: &Device,
    base: Rows,
    queries: Rows,
    selves: &[u32],
    coverage: Option<Coverage>,
    k: usize,
) -> Result<(Vec<u32>, Vec<f32>)> {
    let (bases, rows) = (base.variance.len(), queries.variance.len());
    let width = base.values.len() / bases.max(1);
    let floats = |count: usize| device.alloc(count.max(1) * size_of::<f32>());
    let (base_at, base_variance) = (floats(bases * width)?, floats(bases)?);
    device.upload(base_at, base.values)?;
    device.upload(base_variance, base.variance)?;
    let (query_at, query_variance, own) = (floats(rows * width)?, floats(rows)?, floats(rows)?);
    device.upload(query_at, queries.values)?;
    device.upload(query_variance, queries.variance)?;
    device.upload(own, selves)?;
    let query_rows = floats(rows)?;
    let mixed = match coverage {
        Some(coverage) => {
            device.upload(query_rows, coverage.rows)?;
            Mixed {
                moments: reals(device, coverage.moments)?,
                presence: reals(device, coverage.presence)?,
                stride: coverage.presence.len() as i32,
                samples: coverage.samples as i32,
                terms: reals(device, coverage.terms)?,
            }
        }
        None => Mixed::default(),
    };
    let (distances, indices) = (floats(rows * k)?, floats(rows * k)?);
    device.upload(distances, &vec![f32::INFINITY; rows * k])?;
    device.upload(indices, &vec![NOBODY; rows * k])?;
    let dots = floats(BASE_TILE * QUERY_TILE)?;
    for start in (0..rows).step_by(QUERY_TILE) {
        let tile = QUERY_TILE.min(rows - start);
        let query = Table::new(query_at.after(start * width), tile, width);
        for first in (0..bases).step_by(BASE_TILE) {
            let span = BASE_TILE.min(bases - first);
            let base = Table::new(base_at.after(first * width), span, width);
            let measured = Table::new(dots, tile, span);
            device.product(1.0, (query, false), (base, true), 0.0, measured)?;
            device.nearest(
                spread(tile * 32),
                dots,
                tile as i32,
                span as i32,
                first as i32,
                own.after(start),
                query_variance.after(start),
                base_variance,
                k as i32,
                distances.after(start * k),
                indices.after(start * k),
                query_rows.after(start),
                mixed.moments,
                mixed.presence,
                mixed.stride,
                mixed.samples,
                mixed.terms,
            )?;
        }
    }
    Ok((
        device.download(indices, rows * k)?,
        device.download(distances, rows * k)?,
    ))
}
