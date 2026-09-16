use std::collections::HashMap;
use std::collections::hash_map::Entry;

pub const DOMAIN_TARGET: usize = 0;
pub const DOMAIN_MODEL: usize = 3;
pub const DOMAIN_MODEL_LENGTH: usize = 5;
pub const DOMAIN_SEQUENCE_SCORE: usize = 7;
pub const DOMAIN_SCORE: usize = 13;
pub const DOMAIN_HMM_FROM: usize = 15;
pub const DOMAIN_HMM_TO: usize = 16;

pub type Hits = HashMap<usize, (String, f64)>;

/// The table is read left to right, so the wanted columns come off one pass over the line
/// rather than collecting every field of every row.
pub struct Columns<'a> {
    fields: std::str::SplitWhitespace<'a>,
    at: usize,
}

impl<'a> Columns<'a> {
    pub fn new(line: &'a str) -> Self {
        Self {
            fields: line.split_whitespace(),
            at: 0,
        }
    }

    pub fn at(&mut self, column: usize) -> Option<&'a str> {
        while self.at < column {
            self.fields.next()?;
            self.at += 1;
        }
        self.at += 1;
        self.fields.next()
    }
}

pub fn rows(table: &str) -> impl Iterator<Item = Columns<'_>> {
    table
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(Columns::new)
}

/// A gene that trips two models is one gene, so counting it under both would inflate presence
/// and duplication at once. The name settles a tie, since the order hmmsearch lists its rows in
/// is not something the answer should depend on.
pub fn keep_best(best: &mut Hits, protein: usize, model: &str, score: f64) {
    match best.entry(protein) {
        Entry::Occupied(mut held) => {
            let (kept, top) = held.get();
            if score > *top || (score == *top && model < kept.as_str()) {
                held.insert((model.to_string(), score));
            }
        }
        Entry::Vacant(slot) => {
            slot.insert((model.to_string(), score));
        }
    }
}
