use rand::seq::SliceRandom;
use rand::{Rng, rngs::StdRng};

const BANDS: usize = 10;

// A fast run draws each contig with chance its last loss over the mean, so the contigs the space
// already places are skipped and the steps go where it still learns (Katharopoulos & Fleuret 2018).
// The mean is taken over contigs of like length, since a short contig's views are noisier.
pub struct Draw {
    losses: Option<Vec<f64>>,
    bands: Vec<u8>,
    order: Vec<u32>,
    drawn: Vec<u32>,
}

fn bands(lengths: &[usize]) -> Vec<u8> {
    let mut sorted = lengths.to_vec();
    sorted.sort_unstable();
    let cuts = (1..BANDS)
        .map(|band| sorted[band * sorted.len() / BANDS])
        .collect::<Vec<_>>();
    lengths
        .iter()
        .map(|length| cuts.partition_point(|cut| cut <= length) as u8)
        .collect()
}

impl Draw {
    pub fn new(lengths: &[usize], fast: bool) -> Self {
        let items = lengths.len();
        Self {
            losses: fast.then(|| vec![f64::NAN; items]),
            bands: if fast { bands(lengths) } else { Vec::new() },
            order: (0..items as u32).collect(),
            drawn: Vec::with_capacity(items),
        }
    }

    pub fn fast(&self) -> bool {
        self.losses.is_some()
    }

    pub fn note(&self, drawn: usize, all: usize) -> String {
        match self.fast() {
            true => format!(", {:.0}% drawn", 100.0 * drawn as f64 / all.max(1) as f64),
            false => String::new(),
        }
    }

    pub fn epoch(&mut self, rng: &mut StdRng) -> &[u32] {
        self.order.shuffle(rng);
        let Some(losses) = &self.losses else {
            return &self.order;
        };
        let mut sums = [(0.0, 0); BANDS];
        for (loss, band) in losses.iter().zip(&self.bands) {
            if loss.is_finite() {
                let sum = &mut sums[*band as usize];
                (sum.0, sum.1) = (sum.0 + loss, sum.1 + 1);
            }
        }
        let means = sums.map(|(sum, known)| sum / f64::from(known.max(1)));
        self.drawn.clear();
        for item in &self.order {
            let loss = losses[*item as usize];
            let mean = means[self.bands[*item as usize] as usize];
            let chance = match loss.is_finite() && mean > 0.0 {
                true => (loss / mean).min(1.0),
                false => 1.0,
            };
            if rng.random::<f64>() < chance {
                self.drawn.push(*item);
            }
        }
        &self.drawn
    }

    pub fn record(&mut self, labels: &[u32], rows: &[(f64, usize)], offset: u32) {
        let Some(losses) = &mut self.losses else {
            return;
        };
        let mut sums = std::collections::HashMap::<u32, (f64, usize)>::new();
        for (label, (loss, pairs)) in labels.iter().zip(rows) {
            if *label >= offset && *pairs > 0 {
                let entry = sums.entry(label - offset).or_default();
                entry.0 += loss;
                entry.1 += pairs;
            }
        }
        for (item, (loss, pairs)) in sums {
            losses[item as usize] = loss / pairs as f64;
        }
    }
}
