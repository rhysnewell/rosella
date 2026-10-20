use std::rc::Rc;

use anyhow::Result;

use crate::learned::device::{Device, Ptr, Table, columns, spread};
use crate::learned::net::Net;
use crate::learned::norm::{BLOCK, EPSILON, MOMENTUM, dropout_bar, dropout_scale};
use crate::learned::optim::Clock;

const FLOAT: usize = size_of::<f32>();

struct Layer {
    units: usize,
    from: usize,
    weights: usize,
    scale: usize,
    shift: usize,
    running: usize,
    normed: Ptr,
    active: Ptr,
    inverse: Ptr,
    mean: Ptr,
}

// Every tensor sits in one flat table in the order the host's slices run, so the optimiser
// is one launch and the weights come back in one copy.
pub(crate) struct DeviceNet {
    device: Rc<Device>,
    layers: Vec<Layer>,
    width: usize,
    last: usize,
    last_bias: usize,
    sizes: Vec<usize>,
    params: Ptr,
    grads: Ptr,
    first: Ptr,
    second: Ptr,
    running: Ptr,
    output: Ptr,
    lengths: Ptr,
    pub(crate) gradient: Ptr,
    pub(crate) input_gradient: Ptr,
    scratch: [Ptr; 2],
    partial: Ptr,
    input: Table,
    rows: usize,
    kept: f32,
    clock: Clock,
}

impl DeviceNet {
    pub(crate) fn new(
        device: Rc<Device>,
        net: &Net,
        rows: usize,
        decay: f64,
        input_from: Option<usize>,
    ) -> Result<Self> {
        let floats = |count: usize| device.alloc(count * FLOAT);
        let tensors = net.params.slices();
        let sizes = tensors
            .iter()
            .map(|tensor| tensor.len())
            .collect::<Vec<_>>();
        let offsets = sizes
            .iter()
            .scan(0, |at, size| Some(std::mem::replace(at, *at + size)))
            .collect::<Vec<_>>();
        let count = (sizes.len() - 2) / 3;
        let mut running_at = 0;
        let mut layers = Vec::with_capacity(count);
        for at in 0..count {
            let units = sizes[3 * at + 1];
            layers.push(Layer {
                units,
                from: sizes[3 * at] / units,
                weights: offsets[3 * at],
                scale: offsets[3 * at + 1],
                shift: offsets[3 * at + 2],
                running: running_at,
                normed: floats(rows * units)?,
                active: floats(rows * units)?,
                inverse: floats(units)?,
                mean: floats(units)?,
            });
            running_at += 2 * units;
        }
        let width = sizes[3 * count + 1];
        let widest = layers
            .iter()
            .map(|layer| layer.units)
            .fold(width, usize::max);
        let total = offsets[3 * count + 1] + width;
        let params = floats(total)?;
        device.upload(params, &tensors.concat())?;
        let stats = net
            .running()
            .iter()
            .flat_map(|running| running.parts().concat())
            .collect::<Vec<_>>();
        let running = floats(stats.len())?;
        device.upload(running, &stats)?;
        let input_width = layers
            .first()
            .map_or(sizes[3 * count] / width, |layer| layer.from);
        let input_gradient = rows * input_from.map_or(0, |from| input_width - from);
        Ok(Self {
            width,
            last: offsets[3 * count],
            last_bias: offsets[3 * count + 1],
            sizes,
            grads: floats(total)?,
            first: floats(total)?,
            second: floats(total)?,
            params,
            running,
            output: floats(rows * width)?,
            lengths: floats(rows)?,
            gradient: floats(rows * width)?,
            input_gradient: floats(input_gradient)?,
            scratch: [floats(rows * widest)?, floats(rows * widest)?],
            partial: floats(rows.div_ceil(BLOCK) * widest)?,
            input: Table::new(Ptr::default(), 0, input_width),
            rows: 0,
            kept: 1.0,
            clock: Clock::new(decay),
            layers,
            device,
        })
    }

    pub(crate) fn input_width(&self) -> usize {
        self.input.cols
    }

    fn param(&self, offset: usize, rows: usize, cols: usize) -> Table {
        Table::new(self.params.after(offset), rows, cols)
    }

    fn grad(&self, offset: usize, rows: usize, cols: usize) -> Table {
        Table::new(self.grads.after(offset), rows, cols)
    }

    fn below(&self, at: usize) -> Table {
        match at {
            0 => self.input,
            _ => {
                let layer = &self.layers[at - 1];
                Table::new(layer.active, self.rows, layer.units)
            }
        }
    }

    pub(crate) fn output(&self) -> Table {
        Table::new(self.output, self.rows, self.width)
    }

    pub(crate) fn gradient(&self) -> Table {
        Table::new(self.gradient, self.rows, self.width)
    }

    pub(crate) fn forward(&mut self, input: Table, drop: Option<(f32, u64)>) -> Result<()> {
        let device = self.device.clone();
        (self.input, self.rows, self.kept) = (input, input.rows, dropout_scale(drop));
        for at in 0..self.layers.len() {
            let layer = &self.layers[at];
            let normed = Table::new(layer.normed, input.rows, layer.units);
            let weights = self.param(layer.weights, layer.units, layer.from);
            device.product(1.0, (self.below(at), false), (weights, true), 0.0, normed)?;
            self.normalise(at, drop.filter(|_| at == 0))?;
        }
        let (output, top) = (self.output(), self.below(self.layers.len()));
        let (n, width) = (input.rows as i32, self.width as i32);
        let bias = self.params.after(self.last_bias);
        device.broadcast(spread(input.rows * self.width), output.at, bias, n, width)?;
        let last = self.param(self.last, self.width, top.cols);
        device.product(1.0, (top, false), (last, true), 1.0, output)?;
        device.unit(spread(input.rows), output.at, self.lengths, n, width)
    }

    fn normalise(&self, at: usize, drop: Option<(f32, u64)>) -> Result<()> {
        let device = &self.device;
        let layer = &self.layers[at];
        let (rows, units) = (self.rows, layer.units);
        let blocks = rows.div_ceil(BLOCK);
        let (n, u, b, block) = (rows as i32, units as i32, blocks as i32, BLOCK as i32);
        let (normed, partial, mean, none) =
            (layer.normed, self.partial, layer.mean, Ptr::default());
        let grid = columns(units, blocks);
        device.block_sums(grid, normed, none, none, n, u, block, 0, partial)?;
        device.fold(spread(units), partial, b, u, rows as f32, mean)?;
        device.block_sums(grid, normed, none, mean, n, u, block, 1, partial)?;
        let running = self.running.after(layer.running);
        let variance = running.after(units);
        let epsilon = EPSILON as f32;
        device.finish(
            spread(units),
            partial,
            b,
            u,
            n,
            mean,
            MOMENTUM,
            epsilon,
            running,
            variance,
            layer.inverse,
        )?;
        let (dropping, bar, key) = dropout_bar(drop).map_or((0, 0, 0), |(bar, key)| (1, bar, key));
        let (scale, shift) = (
            self.params.after(layer.scale),
            self.params.after(layer.shift),
        );
        device.normalise(
            spread(rows * units),
            normed,
            layer.active,
            mean,
            layer.inverse,
            scale,
            shift,
            n,
            u,
            dropping,
            bar,
            key,
            dropout_scale(drop),
        )
    }

    fn sums(&self, table: Table, times: Option<Ptr>, out: Ptr) -> Result<()> {
        let blocks = table.rows.div_ceil(BLOCK);
        let (n, cols, block) = (table.rows as i32, table.cols as i32, BLOCK as i32);
        let mode = if times.is_some() { 2 } else { 0 };
        let (times, none) = (times.unwrap_or_default(), Ptr::default());
        let grid = columns(table.cols, blocks);
        let device = &self.device;
        device.block_sums(
            grid,
            table.at,
            times,
            none,
            n,
            cols,
            block,
            mode,
            self.partial,
        )?;
        device.fold(
            spread(table.cols),
            self.partial,
            blocks as i32,
            cols,
            1.0,
            out,
        )
    }

    pub(crate) fn backward(&mut self, input_from: Option<usize>) -> Result<()> {
        let device = self.device.clone();
        let (rows, width) = (self.rows, self.width);
        let gradient = self.gradient();
        let (n, w) = (rows as i32, width as i32);
        device.unit_backward(spread(rows), gradient.at, self.output, self.lengths, n, w)?;
        let top = self.below(self.layers.len());
        let last_grad = self.grad(self.last, width, top.cols);
        device.product(1.0, (gradient, true), (top, false), 0.0, last_grad)?;
        self.sums(gradient, None, self.grads.after(self.last_bias))?;
        let mut below = Table::new(self.scratch[0], rows, top.cols);
        let last = self.param(self.last, width, top.cols);
        device.product(1.0, (gradient, false), (last, false), 0.0, below)?;
        let mut spare = 1;
        for at in (0..self.layers.len()).rev() {
            let layer = &self.layers[at];
            let (units, count) = (layer.units, rows * layer.units);
            let kept = if at == 0 { self.kept } else { 1.0 };
            device.mask(spread(count), below.at, layer.active, count as i32, kept)?;
            let (shift, scaled) = (self.grads.after(layer.shift), self.grads.after(layer.scale));
            self.sums(below, None, shift)?;
            self.sums(below, Some(layer.normed), scaled)?;
            let scale = self.params.after(layer.scale);
            let (normed, inverse, u) = (layer.normed, layer.inverse, units as i32);
            device.rewrite(
                spread(count),
                below.at,
                normed,
                scale,
                inverse,
                shift,
                scaled,
                n,
                u,
            )?;
            let weights = self.param(layer.weights, units, layer.from);
            let next = match (at, input_from) {
                (0, None) => None,
                (0, Some(from)) => {
                    let out = Table::new(self.input_gradient, rows, layer.from - from);
                    let rest = weights.cols_from(from);
                    device.product(1.0, (below, false), (rest, false), 0.0, out)?;
                    None
                }
                _ => {
                    let out = Table::new(self.scratch[spare], rows, layer.from);
                    device.product(1.0, (below, false), (weights, false), 0.0, out)?;
                    Some(out)
                }
            };
            let weight_grad = self.grad(layer.weights, units, layer.from);
            device.product(
                1.0,
                (below, true),
                (self.below(at), false),
                0.0,
                weight_grad,
            )?;
            if let Some(next) = next {
                (below, spare) = (next, 1 - spare);
            }
        }
        Ok(())
    }

    pub(crate) fn update(&mut self, rate: f64, momentum: f64) -> Result<()> {
        let factors = self.clock.tick(rate, momentum);
        let count = self.sizes.iter().sum::<usize>();
        self.device.adamw(
            spread(count),
            self.params,
            self.grads,
            self.first,
            self.second,
            count as i32,
            factors.shrink,
            factors.step,
            factors.momentum,
            factors.keep,
            factors.correction,
            factors.epsilon,
        )
    }

    fn split(&self, flat: &[f32]) -> Vec<Vec<f32>> {
        let mut at = 0;
        self.sizes
            .iter()
            .map(|size| {
                at += size;
                flat[at - size..at].to_vec()
            })
            .collect()
    }

    pub(crate) fn grads(&self) -> Result<Vec<Vec<f32>>> {
        let total = self.sizes.iter().sum();
        Ok(self.split(&self.device.download(self.grads, total)?))
    }

    pub(crate) fn copy_into(&self, net: &mut Net) -> Result<()> {
        let total = self.sizes.iter().sum();
        let params = self.split(&self.device.download(self.params, total)?);
        for (tensor, values) in net.params.slices_mut().into_iter().zip(params) {
            tensor.copy_from_slice(&values);
        }
        let units = self.layers.iter().map(|layer| 2 * layer.units).sum();
        let stats = self.device.download::<f32>(self.running, units)?;
        let mut at = 0;
        for running in net.running_mut() {
            for part in running.parts_mut() {
                part.copy_from_slice(&stats[at..at + part.len()]);
                at += part.len();
            }
        }
        Ok(())
    }
}
