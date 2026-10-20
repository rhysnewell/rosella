use std::rc::Rc;

use anyhow::{Result, ensure};
use cudarc::driver::CudaSlice;
use ndarray::Array2;

use crate::learned::VIEWS_PER_STEP;
use crate::learned::cuda::kernels::{Table, spread};
use crate::learned::cuda::loss::Loss;
use crate::learned::cuda::net::DeviceNet;
use crate::learned::cuda::{Device, current};
use crate::learned::towers::{Towers, lift};

pub struct Trainer {
    device: Rc<Device>,
    joint: DeviceNet,
    depth: DeviceNet,
    loss: Loss,
    split: usize,
    width: usize,
    features: CudaSlice<f32>,
    joined: CudaSlice<f32>,
}

impl Trainer {
    pub fn new(towers: &Towers, decay: f64) -> Result<Option<Self>> {
        let Some(device) = current() else {
            return Ok(None);
        };
        let (rows, split) = (VIEWS_PER_STEP, towers.split);
        let joint = DeviceNet::new(device.clone(), &towers.joint, rows, decay, Some(split))?;
        let depth = DeviceNet::new(device.clone(), &towers.depth, rows, decay, None)?;
        let width = split + depth.input_width();
        Ok(Some(Self {
            features: device.zeros(rows * width)?,
            joined: device.zeros(rows * joint.input_width())?,
            loss: Loss::new(device.clone(), rows)?,
            joint,
            depth,
            split,
            width,
            device,
        }))
    }

    fn launch(
        &mut self,
        features: &Array2<f32>,
        labels: &[u32],
        drop: Option<(f32, u64)>,
        temperature: f32,
    ) -> Result<Vec<usize>> {
        let rows = features.nrows();
        ensure!(
            rows <= VIEWS_PER_STEP && features.ncols() == self.width,
            "a batch of {rows} by {} does not fit the device tables",
            features.ncols()
        );
        let device = self.device.clone();
        let values = features.as_slice().expect("a batch is contiguous");
        device.upload(values, &mut self.features)?;
        let table = Table::new(device.at(&self.features), rows, self.width);
        let below = drop.map(|(share, key)| (share, crate::digest::mix(key ^ u64::MAX)));
        self.depth.forward(table.from_column(self.split), below)?;
        let depth = self.depth.output();
        let joined = Table::new(device.at(&self.joined), rows, self.split + depth.cols);
        device.join(
            spread(rows * joined.cols),
            table.at,
            self.width as i32,
            self.split as i32,
            depth.at,
            depth.cols as i32,
            lift(),
            joined.at,
            rows as i32,
        )?;
        self.joint.forward(joined, drop)?;
        let gradient = self.joint.gradient();
        let owns = self
            .loss
            .launch(self.joint.output(), labels, temperature, gradient)?;
        self.joint.backward(Some(self.split))?;
        let count = rows * depth.cols;
        let (from, to) = (
            device.at(&self.joint.input_gradient),
            device.at(&self.depth.gradient),
        );
        device.scaled_copy(spread(count), from, to, count as i64, lift())?;
        self.depth.backward(None)?;
        Ok(owns)
    }

    pub fn gradients(
        &mut self,
        features: &Array2<f32>,
        labels: &[u32],
        drop: Option<(f32, u64)>,
        temperature: f32,
    ) -> Result<(f64, Vec<(f64, usize)>)> {
        let owns = self.launch(features, labels, drop, temperature)?;
        self.loss.read(owns)
    }

    pub fn step(
        &mut self,
        features: &Array2<f32>,
        labels: &[u32],
        drop: Option<(f32, u64)>,
        (rate, momentum): (f64, f64),
        temperature: f32,
    ) -> Result<(f64, Vec<(f64, usize)>)> {
        let owns = self.launch(features, labels, drop, temperature)?;
        self.joint.update(rate, momentum)?;
        self.depth.update(rate, momentum)?;
        self.loss.read(owns)
    }

    pub fn grads(&self) -> Result<[Vec<Vec<f32>>; 2]> {
        Ok([self.joint.grads()?, self.depth.grads()?])
    }

    pub fn copy_into(&self, towers: &mut Towers) -> Result<()> {
        self.joint.copy_into(&mut towers.joint)?;
        self.depth.copy_into(&mut towers.depth)
    }
}
