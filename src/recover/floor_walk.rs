use std::ops::Range;

// Each band holds as many contigs as are already in play, so a walk that stops costs at most
// what it kept and the search never needs a length grid.
pub fn bands(lengths: &[usize], in_play: usize) -> Vec<Range<usize>> {
    let mut bands = Vec::new();
    let (mut start, mut in_play) = (0, in_play.max(1));
    while start < lengths.len() {
        let mut end = (start + in_play).min(lengths.len());
        while end < lengths.len() && lengths[end] == lengths[end - 1] {
            end += 1;
        }
        bands.push(start..end);
        in_play += end - start;
        start = end;
    }
    bands
}

pub struct Reach {
    kept: f64,
    bar: f64,
    long_hits: f64,
    carried: f64,
}

impl Reach {
    pub fn new(kept: f64, bar: f64, long_hits: usize) -> Self {
        Self {
            kept,
            bar,
            long_hits: long_hits.max(1) as f64,
            carried: 0.0,
        }
    }

    pub fn add(&mut self, hits: usize) -> f64 {
        let share = hits as f64 / self.long_hits;
        let before = 1.0 + self.carried;
        self.carried += share;
        self.kept * ((before + share).powi(2) - before.powi(2))
    }

    pub fn bar(&self) -> f64 {
        self.bar
    }
}

// The foreign share pools from the cutoff down because a band under 750 bp holds a few dozen
// marker contigs, too few to read alone. One half is where attaching stops paying.
#[derive(Default)]
pub struct Foreign {
    repeats: f64,
    complete: f64,
}

impl Foreign {
    pub fn admits(&mut self, repeats: usize, complete: f64) -> bool {
        self.repeats += repeats as f64;
        self.complete += complete;
        self.share().is_some_and(|share| share < 0.5)
    }

    pub fn share(&self) -> Option<f64> {
        (self.complete > 0.0).then(|| self.repeats / self.complete)
    }
}
