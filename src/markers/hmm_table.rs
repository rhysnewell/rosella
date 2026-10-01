use std::collections::HashMap;
use std::collections::hash_map::Entry;

pub const DOMAIN_TARGET: usize = 0;
pub const DOMAIN_MODEL: usize = 3;
pub const DOMAIN_MODEL_LENGTH: usize = 5;
pub const DOMAIN_SEQUENCE_E_VALUE: usize = 6;
pub const DOMAIN_SEQUENCE_SCORE: usize = 7;
pub const DOMAIN_I_E_VALUE: usize = 12;
pub const DOMAIN_SCORE: usize = 13;
pub const DOMAIN_HMM_FROM: usize = 15;
pub const DOMAIN_HMM_TO: usize = 16;
pub const DOMAIN_ALI_FROM: usize = 17;
pub const DOMAIN_ALI_TO: usize = 18;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Reach {
    pub model_from: u32,
    pub model_to: u32,
    pub model_length: u32,
    pub protein_from: u32,
    pub protein_to: u32,
}

impl Reach {
    fn widen(self, other: Self) -> Self {
        Self {
            model_from: self.model_from.min(other.model_from),
            model_to: self.model_to.max(other.model_to),
            protein_from: self.protein_from.min(other.protein_from),
            protein_to: self.protein_to.max(other.protein_to),
            ..self
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Best {
    pub model: String,
    pub score: f64,
    pub reach: Reach,
}

pub type Hits = HashMap<usize, Best>;

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

pub struct Domain<'a> {
    pub protein: usize,
    pub model: &'a str,
    pub sequence_e_value: Option<f64>,
    pub sequence_score: f64,
    pub i_e_value: Option<f64>,
    pub score: f64,
    pub reach: Reach,
}

impl Domain<'_> {
    pub fn span(&self) -> f64 {
        let covered = (self.reach.model_to + 1).saturating_sub(self.reach.model_from);
        (f64::from(covered) / f64::from(self.reach.model_length)).clamp(0.0, 1.0)
    }
}

// The E-values are optional because only CheckM's reading settles on them.
pub fn domains(table: &str) -> impl Iterator<Item = Domain<'_>> {
    rows(table).filter_map(|mut fields| {
        let protein = fields.at(DOMAIN_TARGET)?.parse::<usize>().ok()?;
        let model = fields.at(DOMAIN_MODEL)?;
        let model_length = fields.at(DOMAIN_MODEL_LENGTH)?.parse::<u32>().ok()?;
        let sequence_e_value = fields.at(DOMAIN_SEQUENCE_E_VALUE)?.parse().ok();
        let sequence_score = fields.at(DOMAIN_SEQUENCE_SCORE)?.parse::<f64>().ok()?;
        let i_e_value = fields.at(DOMAIN_I_E_VALUE)?.parse().ok();
        let score = fields.at(DOMAIN_SCORE)?.parse::<f64>().ok()?;
        let mut position = |column| fields.at(column)?.parse::<u32>().ok();
        let reach = Reach {
            model_length,
            model_from: position(DOMAIN_HMM_FROM)?,
            model_to: position(DOMAIN_HMM_TO)?,
            protein_from: position(DOMAIN_ALI_FROM)?,
            protein_to: position(DOMAIN_ALI_TO)?,
        };
        (model_length > 0).then_some(Domain {
            protein,
            model,
            sequence_e_value,
            sequence_score,
            i_e_value,
            score,
            reach,
        })
    })
}

pub fn floor(bars: impl Iterator<Item = (f64, f64)>, scale: f64) -> f64 {
    let lowest = bars
        .map(|(sequence, domain)| sequence.min(domain))
        .fold(f64::INFINITY, f64::min);
    match lowest.is_finite() {
        true => (lowest * scale * 100.0).floor() / 100.0,
        false => 0.0,
    }
}

/// A gene that trips two models is one gene, so counting it under both would inflate presence
/// and duplication at once. The name settles a tie, since the order hmmsearch lists its rows in
/// is not something the answer should depend on.
pub fn keep_best(best: &mut Hits, protein: usize, model: &str, score: f64, reach: Reach) {
    match best.entry(protein) {
        Entry::Occupied(mut held) => {
            let held = held.get_mut();
            if held.model == model {
                held.score = held.score.max(score);
                held.reach = held.reach.widen(reach);
            } else if score > held.score || (score == held.score && model < held.model.as_str()) {
                *held = Best {
                    model: model.to_string(),
                    score,
                    reach,
                };
            }
        }
        Entry::Vacant(slot) => {
            slot.insert(Best {
                model: model.to_string(),
                score,
                reach,
            });
        }
    }
}
