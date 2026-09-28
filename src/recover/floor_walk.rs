use std::collections::HashMap;
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

pub struct Seen {
    pub bin: usize,
    pub band: usize,
    pub repeats: bool,
    pub complete: f64,
}

// Past the band where the assembly's short contigs turn foreign, a bin keeps walking only
// while its own markers still vouch for what it takes.
pub fn deepest(
    seen: &[Seen],
    stop: usize,
    bands: usize,
    worth: f64,
) -> (f64, HashMap<usize, Option<usize>>) {
    let judged = seen.iter().filter(|seen| seen.band <= stop);
    let repeats = judged.clone().filter(|seen| seen.repeats).count() as f64;
    let complete = judged.clone().map(|seen| seen.complete).sum::<f64>();
    let foreign = (repeats / complete).min(1.0);
    let trade = |seen: &Seen| match seen.repeats {
        true => -worth,
        false => 1.0 - foreign * (1.0 - seen.complete),
    };
    let mut gain = HashMap::<usize, f64>::new();
    for seen in judged {
        *gain.entry(seen.bin).or_default() += trade(seen);
    }
    let mut deepest = gain
        .iter()
        .map(|(bin, gain)| (*bin, (*gain > 0.0).then_some(stop)))
        .collect::<HashMap<_, _>>();
    for band in stop + 1..bands {
        let mut step = HashMap::<usize, f64>::new();
        for seen in seen.iter().filter(|seen| seen.band == band) {
            *step.entry(seen.bin).or_default() += trade(seen);
        }
        for (bin, at) in deepest.iter_mut() {
            if *at != Some(band - 1) {
                continue;
            }
            let total = gain.get_mut(bin).expect("every judged bin has a gain");
            *total += step.get(bin).copied().unwrap_or(0.0);
            if *total > 0.0 {
                *at = Some(band);
            }
        }
    }
    (foreign, deepest)
}
