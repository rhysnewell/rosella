use std::collections::HashMap;
use std::rc::Rc;

use anyhow::{Result, ensure};

use crate::learned::device::{Device, Ptr, Table, per_row};

const MOST_MATES: usize = 32;
const INT: usize = size_of::<i32>();

#[derive(Clone, Copy)]
pub(crate) enum Kind {
    PerPositive,
    Supervised,
}

pub(crate) struct Loss {
    device: Rc<Device>,
    kind: Kind,
    weights: Ptr,
    start: Ptr,
    size: Ptr,
    members: Ptr,
    row_loss: Ptr,
    owns: Vec<usize>,
}

impl Loss {
    pub(crate) fn new(device: Rc<Device>, rows: usize, kind: Kind) -> Result<Self> {
        Ok(Self {
            kind,
            weights: device.alloc(rows * rows * size_of::<f32>())?,
            start: device.alloc(rows * INT)?,
            size: device.alloc(rows * INT)?,
            members: device.alloc(rows * INT)?,
            row_loss: device.alloc(rows * size_of::<f32>())?,
            owns: Vec::new(),
            device,
        })
    }

    // The groups go up before any of the step is queued, so no upload waits on the device.
    pub(crate) fn prepare(&mut self, labels: &[u32]) -> Result<()> {
        let mut groups = HashMap::<u32, Vec<i32>>::new();
        for (at, label) in labels.iter().enumerate() {
            groups.entry(*label).or_default().push(at as i32);
        }
        let mut starts = HashMap::with_capacity(groups.len());
        let mut members = Vec::<i32>::with_capacity(labels.len());
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
            matches!(self.kind, Kind::Supervised)
                || sizes.iter().all(|size| *size as usize <= MOST_MATES),
            "a contig has more views than the device loss holds"
        );
        let start = labels.iter().map(|label| starts[label]).collect::<Vec<_>>();
        self.device.upload(self.start, &start)?;
        self.device.upload(self.size, &sizes)?;
        self.device.upload(self.members, &members)?;
        self.owns = sizes
            .iter()
            .map(|size| match self.kind {
                Kind::PerPositive => *size as usize - 1,
                Kind::Supervised => usize::from(*size > 1),
            })
            .collect();
        Ok(())
    }

    pub(crate) fn launch(&self, rows: Table, temperature: f32, gradient: Table) -> Result<()> {
        let n = self.owns.len();
        let device = &self.device;
        let weights = Table::new(self.weights, n, n);
        device.product(1.0 / temperature, (rows, false), (rows, true), 0.0, weights)?;
        let pairs = self.owns.iter().sum::<usize>().max(1);
        let scale = 1.0 / (temperature * pairs as f32);
        let kernel = match self.kind {
            Kind::PerPositive => Device::per_positive,
            Kind::Supervised => Device::contrastive,
        };
        kernel(
            device,
            per_row(n),
            weights.at,
            n as i32,
            self.start,
            self.size,
            self.members,
            1.0 / temperature,
            scale,
            self.row_loss,
        )?;
        device.product(1.0, (weights, false), (rows, false), 0.0, gradient)?;
        device.product(1.0, (weights, true), (rows, false), 1.0, gradient)
    }

    pub(crate) fn read(&self) -> Result<(f64, Vec<(f64, usize)>)> {
        let losses = self
            .device
            .download::<f32>(self.row_loss, self.owns.len())?;
        let pairs = self.owns.iter().sum::<usize>().max(1);
        let rows = losses
            .iter()
            .map(|loss| f64::from(*loss))
            .zip(self.owns.iter().copied())
            .collect::<Vec<_>>();
        let loss = rows.iter().map(|(loss, _)| loss).sum::<f64>() / pairs as f64;
        Ok((loss, rows))
    }
}
