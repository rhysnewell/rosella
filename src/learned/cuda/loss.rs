use std::collections::HashMap;
use std::rc::Rc;

use anyhow::{Result, ensure};
use cudarc::driver::CudaSlice;

use crate::learned::cuda::Device;
use crate::learned::cuda::kernels::{Table, per_row};
use crate::learned::loss::LANES;

const MOST_MATES: usize = 32;

pub(crate) struct Loss {
    device: Rc<Device>,
    weights: CudaSlice<f32>,
    start: CudaSlice<i32>,
    size: CudaSlice<i32>,
    members: CudaSlice<i32>,
    row_loss: CudaSlice<f64>,
}

impl Loss {
    pub(crate) fn new(device: Rc<Device>, rows: usize) -> Result<Self> {
        Ok(Self {
            weights: device.zeros(rows * rows)?,
            start: device.zeros(rows)?,
            size: device.zeros(rows)?,
            members: device.zeros(rows)?,
            row_loss: device.zeros(rows)?,
            device,
        })
    }

    // The row losses stay on the device until read, so the host waits on the device once a step.
    pub(crate) fn launch(
        &mut self,
        rows: Table,
        labels: &[u32],
        temperature: f32,
        gradient: Table,
    ) -> Result<Vec<usize>> {
        let n = labels.len();
        let mut groups = HashMap::<u32, Vec<i32>>::new();
        for (at, label) in labels.iter().enumerate() {
            groups.entry(*label).or_default().push(at as i32);
        }
        let mut starts = HashMap::with_capacity(groups.len());
        let mut members = Vec::with_capacity(n);
        for label in labels {
            starts.entry(*label).or_insert_with(|| {
                let start = members.len() as i32;
                members.extend(&groups[label]);
                start
            });
        }
        let sizes = labels
            .iter()
            .map(|label| groups[label].len() as i32)
            .collect::<Vec<_>>();
        ensure!(
            sizes.iter().all(|size| *size as usize <= MOST_MATES),
            "a contig has more views than the device loss holds"
        );
        let start = labels.iter().map(|label| starts[label]).collect::<Vec<_>>();
        let device = &self.device;
        device.upload(&start, &mut self.start)?;
        device.upload(&sizes, &mut self.size)?;
        device.upload(&members, &mut self.members)?;
        let weights = Table::new(device.at(&self.weights), n, n);
        device.product(1.0 / temperature, (rows, false), (rows, true), 0.0, weights)?;
        device.per_positive(
            per_row(n),
            weights.at,
            n as i32,
            device.at(&self.start),
            device.at(&self.size),
            device.at(&self.members),
            LANES as i32,
            device.at(&self.row_loss),
        )?;
        let owns = sizes
            .iter()
            .map(|size| *size as usize - 1)
            .collect::<Vec<_>>();
        let pairs = owns.iter().sum::<usize>().max(1);
        let scale = 1.0 / (temperature * pairs as f32);
        device.product(scale, (weights, false), (rows, false), 0.0, gradient)?;
        device.product(scale, (weights, true), (rows, false), 1.0, gradient)?;
        Ok(owns)
    }

    pub(crate) fn read(&self, owns: Vec<usize>) -> Result<(f64, Vec<(f64, usize)>)> {
        let losses = self.device.download(&self.row_loss, owns.len())?;
        let pairs = owns.iter().sum::<usize>().max(1);
        let loss = losses.iter().sum::<f64>() / pairs as f64;
        Ok((loss, losses.into_iter().zip(owns).collect()))
    }
}
