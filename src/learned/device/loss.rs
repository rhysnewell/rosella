use std::collections::HashMap;
use std::rc::Rc;

use anyhow::{Result, ensure};

use crate::learned::device::{Device, Ptr, Table, per_row};
use crate::learned::loss::Mates;

pub(crate) const MOST_MATES: usize = 32;
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
    pulled: Ptr,
    members: Ptr,
    held: usize,
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
            pulled: device.alloc(rows * INT)?,
            members: device.alloc(rows * INT)?,
            held: rows,
            row_loss: device.alloc(rows * size_of::<f64>())?,
            owns: Vec::new(),
            device,
        })
    }

    // The groups go up before any of the step is queued, so no upload waits on the device.
    pub(crate) fn prepare(&mut self, labels: &[u32], mates: &Mates) -> Result<()> {
        let groups = mates.groups(labels);
        let mut starts = HashMap::with_capacity(groups.len());
        let mut members = Vec::<i32>::with_capacity(labels.len());
        for label in labels {
            starts.entry(*label).or_insert_with(|| {
                let start = members.len() as i32;
                members.extend(groups[label].0.iter().map(|row| *row as i32));
                start
            });
        }
        let (sizes, pulled): (Vec<_>, Vec<_>) = labels
            .iter()
            .map(|label| {
                let (rows, pulled) = &groups[label];
                (rows.len() as i32, *pulled as i32)
            })
            .unzip();
        ensure!(
            matches!(self.kind, Kind::Supervised)
                || pulled.iter().all(|pulled| *pulled as usize <= MOST_MATES),
            "a contig has more positives than the device loss holds"
        );
        if members.len() > self.held {
            self.held = members.len().next_power_of_two();
            self.members = self.device.alloc(self.held * INT)?;
        }
        let start = labels.iter().map(|label| starts[label]).collect::<Vec<_>>();
        self.device.upload(self.start, &start)?;
        self.device.upload(self.size, &sizes)?;
        self.device.upload(self.pulled, &pulled)?;
        self.device.upload(self.members, &members)?;
        self.owns = pulled
            .iter()
            .map(|pulled| match self.kind {
                Kind::PerPositive => *pulled as usize - 1,
                Kind::Supervised => usize::from(*pulled > 1),
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
        let (launch, peak) = (per_row(n), 1.0 / temperature);
        match self.kind {
            Kind::PerPositive => device.per_positive(
                launch,
                weights.at,
                n as i32,
                self.start,
                self.size,
                self.pulled,
                self.members,
                peak,
                scale,
                self.row_loss,
            ),
            Kind::Supervised => device.contrastive(
                launch,
                weights.at,
                n as i32,
                self.start,
                self.size,
                self.members,
                peak,
                scale,
                self.row_loss,
            ),
        }?;
        device.product(1.0, (weights, false), (rows, false), 0.0, gradient)?;
        device.product(1.0, (weights, true), (rows, false), 1.0, gradient)
    }

    pub(crate) fn read(&self) -> Result<(f64, Vec<(f64, usize)>)> {
        let losses = self.device.sums(self.row_loss, self.owns.len())?;
        let pairs = self.owns.iter().sum::<usize>().max(1);
        let rows = losses
            .into_iter()
            .zip(self.owns.iter().copied())
            .collect::<Vec<_>>();
        let loss = rows.iter().map(|(loss, _)| loss).sum::<f64>() / pairs as f64;
        Ok((loss, rows))
    }
}
