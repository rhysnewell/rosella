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

// The share a join needs rises until the marker joins above it read under half foreign. Bands pool
// from the cutoff down, since one under 750 bp holds too few marker contigs to read alone.
#[derive(Default)]
pub struct Bar {
    seen: Vec<(f32, bool, f64)>,
}

impl Bar {
    pub fn add(&mut self, seen: impl IntoIterator<Item = (f32, bool, f64)>) -> Option<f32> {
        self.seen.extend(seen);
        self.seen.sort_by(|a, b| b.0.total_cmp(&a.0));
        let (mut repeats, mut complete, mut bar) = (0.0, 0.0, None);
        for (at, (share, repeated, filled)) in self.seen.iter().enumerate() {
            repeats += f64::from(u8::from(*repeated));
            complete += filled;
            let next = self.seen.get(at + 1).map(|next| next.0);
            if next.is_none_or(|next| next < *share) && complete > 0.0 && repeats < complete / 2.0 {
                bar = Some(next.unwrap_or(0.5));
            }
        }
        bar
    }
}
