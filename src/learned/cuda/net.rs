use std::rc::Rc;

use anyhow::Result;
use cudarc::driver::CudaSlice;

use crate::learned::cuda::Device;
use crate::learned::cuda::kernels::{Table, columns, spread};
use crate::learned::net::Net;
use crate::learned::norm::{BLOCK, EPSILON, MOMENTUM, dropout_bar, dropout_scale};
use crate::learned::optim::Clock;

const FLOAT: u64 = size_of::<f32>() as u64;

struct Layer {
    units: usize,
    from: usize,
    weights: usize,
    scale: usize,
    shift: usize,
    running: usize,
    normed: CudaSlice<f32>,
    active: CudaSlice<f32>,
    inverse: CudaSlice<f32>,
    mean: CudaSlice<f32>,
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
    params: CudaSlice<f32>,
    grads: CudaSlice<f32>,
    first: CudaSlice<f32>,
    second: CudaSlice<f32>,
    running: CudaSlice<f32>,
    output: CudaSlice<f32>,
    lengths: CudaSlice<f32>,
    pub(crate) gradient: CudaSlice<f32>,
    pub(crate) input_gradient: CudaSlice<f32>,
    scratch: [CudaSlice<f32>; 2],
    partial: CudaSlice<f64>,
    mean: CudaSlice<f64>,
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
                normed: device.zeros(rows * units)?,
                active: device.zeros(rows * units)?,
                inverse: device.zeros(units)?,
                mean: device.zeros(units)?,
            });
            running_at += 2 * units;
        }
        let width = sizes[3 * count + 1];
        let widest = layers.iter().map(|layer| layer.units).max().unwrap_or(0);
        let widest = widest.max(width);
        let total = offsets[3 * count + 1] + width;
        let mut params = device.zeros(total)?;
        device.upload(&tensors.concat(), &mut params)?;
        let stats = net
            .running()
            .iter()
            .flat_map(|running| running.parts().concat())
            .collect::<Vec<_>>();
        let mut running = device.zeros(stats.len())?;
        device.upload(&stats, &mut running)?;
        let input_width = layers
            .first()
            .map_or(sizes[3 * count] / width, |layer| layer.from);
        let input_gradient = rows * input_from.map_or(0, |from| input_width - from);
        Ok(Self {
            width,
            last: offsets[3 * count],
            last_bias: offsets[3 * count + 1],
            sizes,
            grads: device.zeros(total)?,
            first: device.zeros(total)?,
            second: device.zeros(total)?,
            params,
            running,
            output: device.zeros(rows * width)?,
            lengths: device.zeros(rows)?,
            gradient: device.zeros(rows * width)?,
            input_gradient: device.zeros(input_gradient)?,
            scratch: [device.zeros(rows * widest)?, device.zeros(rows * widest)?],
            partial: device.zeros(rows.div_ceil(BLOCK) * widest)?,
            mean: device.zeros(widest)?,
            input: Table::new(0, 0, input_width),
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
        Table::new(self.param_at(offset), rows, cols)
    }

    fn param_at(&self, offset: usize) -> u64 {
        self.device.at(&self.params) + offset as u64 * FLOAT
    }

    fn grad(&self, offset: usize, rows: usize, cols: usize) -> Table {
        Table::new(self.grad_at(offset), rows, cols)
    }

    fn grad_at(&self, offset: usize) -> u64 {
        self.device.at(&self.grads) + offset as u64 * FLOAT
    }

    fn below(&self, at: usize) -> Table {
        match at {
            0 => self.input,
            _ => {
                let layer = &self.layers[at - 1];
                Table::new(self.device.at(&layer.active), self.rows, layer.units)
            }
        }
    }

    pub(crate) fn output(&self) -> Table {
        Table::new(self.device.at(&self.output), self.rows, self.width)
    }

    pub(crate) fn gradient(&self) -> Table {
        Table::new(self.device.at(&self.gradient), self.rows, self.width)
    }

    pub(crate) fn forward(&mut self, input: Table, drop: Option<(f32, u64)>) -> Result<()> {
        let device = self.device.clone();
        (self.input, self.rows, self.kept) = (input, input.rows, dropout_scale(drop));
        for at in 0..self.layers.len() {
            let layer = &self.layers[at];
            let normed = Table::new(device.at(&layer.normed), input.rows, layer.units);
            let weights = self.param(layer.weights, layer.units, layer.from);
            device.product(1.0, (self.below(at), false), (weights, true), 0.0, normed)?;
            self.normalise(at, drop.filter(|_| at == 0))?;
        }
        let (output, top) = (self.output(), self.below(self.layers.len()));
        let (n, width) = (input.rows as i32, self.width as i32);
        let bias = self.param_at(self.last_bias);
        device.broadcast(spread(input.rows * self.width), output.at, bias, n, width)?;
        let last = self.param(self.last, self.width, top.cols);
        device.product(1.0, (top, false), (last, true), 1.0, output)?;
        let lengths = device.at(&self.lengths);
        device.unit(spread(input.rows), output.at, lengths, n, width)
    }

    fn normalise(&self, at: usize, drop: Option<(f32, u64)>) -> Result<()> {
        let device = &self.device;
        let layer = &self.layers[at];
        let (rows, units) = (self.rows, layer.units);
        let blocks = rows.div_ceil(BLOCK);
        let (n, u, b, block) = (rows as i32, units as i32, blocks as i32, BLOCK as i32);
        let (normed, partial, mean) = (
            device.at(&layer.normed),
            device.at(&self.partial),
            device.at(&self.mean),
        );
        let grid = columns(units, blocks);
        device.block_sums(grid, normed, 0, 0, n, u, block, 0, partial)?;
        device.fold(spread(units), partial, b, u, rows as f64, mean, 0)?;
        device.block_sums(grid, normed, 0, mean, n, u, block, 1, partial)?;
        let running = device.at(&self.running) + layer.running as u64 * FLOAT;
        let variance = running + units as u64 * FLOAT;
        let (inverse, mean32) = (device.at(&layer.inverse), device.at(&layer.mean));
        device.finish(
            spread(units),
            partial,
            b,
            u,
            n,
            mean,
            MOMENTUM,
            EPSILON,
            running,
            variance,
            inverse,
            mean32,
        )?;
        let (dropping, bar, key) = dropout_bar(drop).map_or((0, 0, 0), |(bar, key)| (1, bar, key));
        let active = device.at(&layer.active);
        let (scale, shift) = (self.param_at(layer.scale), self.param_at(layer.shift));
        device.normalise(
            spread(rows * units),
            normed,
            active,
            mean32,
            inverse,
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

    fn sums(&self, table: Table, times: Option<u64>, out: u64) -> Result<()> {
        let blocks = table.rows.div_ceil(BLOCK);
        let partial = self.device.at(&self.partial);
        let (n, cols) = (table.rows as i32, table.cols as i32);
        let mode = if times.is_some() { 2 } else { 0 };
        let grid = columns(table.cols, blocks);
        let block = BLOCK as i32;
        let times = times.unwrap_or(0);
        self.device
            .block_sums(grid, table.at, times, 0, n, cols, block, mode, partial)?;
        let grid = spread(table.cols);
        self.device
            .fold(grid, partial, blocks as i32, cols, 1.0, 0, out)
    }

    pub(crate) fn backward(&mut self, input_from: Option<usize>) -> Result<()> {
        let device = self.device.clone();
        let (rows, width) = (self.rows, self.width);
        let gradient = self.gradient();
        let (output, lengths) = (device.at(&self.output), device.at(&self.lengths));
        device.unit_backward(
            spread(rows),
            gradient.at,
            output,
            lengths,
            rows as i32,
            width as i32,
        )?;
        let top = self.below(self.layers.len());
        let last_grad = self.grad(self.last, width, top.cols);
        device.product(1.0, (gradient, true), (top, false), 0.0, last_grad)?;
        self.sums(gradient, None, self.grad_at(self.last_bias))?;
        let mut below = Table::new(device.at(&self.scratch[0]), rows, top.cols);
        let last = self.param(self.last, width, top.cols);
        device.product(1.0, (gradient, false), (last, false), 0.0, below)?;
        let mut spare = 1;
        for at in (0..self.layers.len()).rev() {
            let layer = &self.layers[at];
            let (units, count) = (layer.units, rows * layer.units);
            let kept = if at == 0 { self.kept } else { 1.0 };
            let (normed, active) = (device.at(&layer.normed), device.at(&layer.active));
            device.mask(spread(count), below.at, active, count as i64, kept)?;
            let (shift, scaled) = (self.grad_at(layer.shift), self.grad_at(layer.scale));
            self.sums(below, None, shift)?;
            self.sums(below, Some(normed), scaled)?;
            let (scale, inverse) = (self.param_at(layer.scale), device.at(&layer.inverse));
            let (n, u) = (rows as i32, units as i32);
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
                    let out = Table::new(device.at(&self.input_gradient), rows, layer.from - from);
                    device.product(
                        1.0,
                        (below, false),
                        (weights.from_column(from), false),
                        0.0,
                        out,
                    )?;
                    None
                }
                _ => {
                    let out = Table::new(device.at(&self.scratch[spare]), rows, layer.from);
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
        let device = &self.device;
        device.adamw(
            spread(count),
            device.at(&self.params),
            device.at(&self.grads),
            device.at(&self.first),
            device.at(&self.second),
            count as i64,
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
        Ok(self.split(&self.device.download(&self.grads, total)?))
    }

    pub(crate) fn copy_into(&self, net: &mut Net) -> Result<()> {
        let total = self.sizes.iter().sum();
        let params = self.split(&self.device.download(&self.params, total)?);
        for (tensor, values) in net.params.slices_mut().into_iter().zip(params) {
            tensor.copy_from_slice(&values);
        }
        let units = self.layers.iter().map(|layer| 2 * layer.units).sum();
        let stats = self.device.download(&self.running, units)?;
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
