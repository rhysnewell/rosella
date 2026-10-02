use std::collections::HashMap;

// CheckM1 leaves contigs this short out of the GC spread (Parks et al. 2015).
const GC_SPREAD_FROM: usize = 1000;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Bases {
    pub gc: u64,
    pub at: u64,
    pub ambiguous: u64,
}

impl Bases {
    // Setting the 0x20 bit lowercases a letter without mapping any other byte onto one.
    pub fn count(sequence: &[u8]) -> Self {
        let of = |wanted: [u8; 2]| {
            sequence
                .chunks(usize::from(u16::MAX))
                .map(|chunk| {
                    let found = chunk
                        .iter()
                        .map(|base| u16::from(wanted.contains(&(base | 0x20))))
                        .sum::<u16>();
                    u64::from(found)
                })
                .sum()
        };
        Self {
            gc: of(*b"gc"),
            at: of(*b"at"),
            ambiguous: of(*b"nn"),
        }
    }

    fn share(&self) -> f64 {
        match self.gc + self.at {
            0 => 0.0,
            called => self.gc as f64 / called as f64,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Composition {
    pub gc: f64,
    pub gc_spread: f64,
    pub ambiguous: u64,
}

// The spread is taken about the bin's pooled GC, not the mean of its contigs, as CheckM1 takes it.
pub fn composition(
    bases: &HashMap<usize, Bases>,
    lengths: &[usize],
    contigs: &[usize],
) -> Option<Composition> {
    let held = contigs
        .iter()
        .map(|contig| Some((bases.get(contig)?, lengths[*contig])))
        .collect::<Option<Vec<_>>>()?;
    let total = held.iter().fold(Bases::default(), |sum, (bases, _)| Bases {
        gc: sum.gc + bases.gc,
        at: sum.at + bases.at,
        ambiguous: sum.ambiguous + bases.ambiguous,
    });
    let gc = total.share();
    let spread = held
        .iter()
        .filter(|(_, length)| *length > GC_SPREAD_FROM)
        .map(|(bases, _)| (bases.share() - gc).powi(2))
        .collect::<Vec<_>>();
    let gc_spread = match spread.len() {
        0 | 1 => 0.0,
        count => (spread.iter().sum::<f64>() / count as f64).sqrt(),
    };
    Some(Composition {
        gc: 100.0 * gc,
        gc_spread: 100.0 * gc_spread,
        ambiguous: total.ambiguous,
    })
}

pub fn n50(lengths: &mut [usize]) -> usize {
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let half = lengths.iter().sum::<usize>().div_ceil(2);
    let mut running = 0;
    for length in lengths.iter() {
        running += length;
        if running >= half {
            return *length;
        }
    }
    0
}
