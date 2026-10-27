use anyhow::Result;

use crate::learned::device::{Device, Table, current, spread};

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

// Every query's k nearest base rows by the composition distance, measured exactly. A query's own
// base row, named in `selves`, is never its neighbour.
pub fn nearest(
    base: Rows,
    queries: Rows,
    selves: &[u32],
    k: usize,
) -> Result<Option<(Vec<u32>, Vec<f32>)>> {
    let Some(device) = current().filter(|_| k <= MOST_NEAREST) else {
        return Ok(None);
    };
    let mark = device.mark();
    let found = measure(&device, base, queries, selves, k);
    device.release(mark)?;
    found.map(Some)
}

fn measure(
    device: &Device,
    base: Rows,
    queries: Rows,
    selves: &[u32],
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
            )?;
        }
    }
    Ok((
        device.download(indices, rows * k)?,
        device.download(distances, rows * k)?,
    ))
}
