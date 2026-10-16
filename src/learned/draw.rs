use rand::seq::SliceRandom;
use rand::{Rng, rngs::StdRng};

// A fast run draws each contig with chance its last loss over the mean, so the contigs the space
// already places are skipped and the steps go where it still learns (Katharopoulos & Fleuret 2018).
pub struct Draw {
    losses: Option<Vec<f64>>,
    order: Vec<u32>,
    drawn: Vec<u32>,
}

impl Draw {
    pub fn new(items: usize, fast: bool) -> Self {
        Self {
            losses: fast.then(|| vec![f64::NAN; items]),
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
        let (sum, known) = losses
            .iter()
            .filter(|loss| loss.is_finite())
            .fold((0.0, 0), |(sum, known), loss| (sum + loss, known + 1));
        let mean = sum / f64::from(known.max(1));
        self.drawn.clear();
        for item in &self.order {
            let loss = losses[*item as usize];
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
