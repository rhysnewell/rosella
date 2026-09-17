use std::collections::HashMap;
use std::path::Path;

use anyhow::Result;

use crate::markers::hmm_table::{self, Hits};

pub const DEFAULT_SPAN: f64 = 0.3;

#[derive(Clone, Copy, Debug)]
pub struct Bar {
    pub sequence: f64,
    pub domain: f64,
}

pub type Bars = HashMap<String, Bar>;

pub fn gathering(hmm: &Path) -> Result<Bars> {
    let text = std::fs::read_to_string(hmm)?;
    let mut bars = Bars::new();
    let mut name: Option<String> = None;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("NAME") {
            name = Some(rest.trim().to_string());
        } else if let Some(rest) = line.strip_prefix("GA") {
            let Some(model) = name.take() else {
                continue;
            };
            let mut cutoffs = rest
                .trim()
                .trim_end_matches(';')
                .split_whitespace()
                .filter_map(|value| value.trim_end_matches(';').parse::<f64>().ok());
            if let Some(sequence) = cutoffs.next() {
                let domain = cutoffs.next().unwrap_or(sequence);
                bars.insert(model, Bar { sequence, domain });
            }
        }
    }
    Ok(bars)
}

/// One search serves both readings, so its reporting floor has to sit under the lowest score
/// either can accept, which is the smallest gathering cutoff scaled by the span.
pub fn floor(bars: &Bars, min_span: f64) -> String {
    let lowest = bars
        .values()
        .map(|bar| bar.sequence.min(bar.domain))
        .fold(f64::INFINITY, f64::min);
    match lowest.is_finite() {
        true => format!("{:.2}", (lowest * min_span * 100.0).floor() / 100.0),
        false => "0".to_string(),
    }
}

pub struct Domain<'a> {
    pub protein: usize,
    pub model: &'a str,
    pub sequence_score: f64,
    pub score: f64,
    pub span: f64,
}

/// A gene cut by a contig end can only align to the part of the model it still carries, so the
/// full length gathering threshold is scaled to the span that could have matched at all.
pub fn accepted(table: &str, bars: &Bars, min_span: f64, keep: impl Fn(usize) -> bool) -> Hits {
    let mut best = Hits::new();
    for domain in parse(table) {
        let Some(bar) = bars.get(domain.model) else {
            continue;
        };
        if !keep(domain.protein) {
            continue;
        }
        if domain.span < min_span || domain.score < bar.sequence * domain.span {
            continue;
        }
        hmm_table::keep_best(&mut best, domain.protein, domain.model, domain.score);
    }
    best
}

/// What a gathering-cutoff search reports. The decision is the whole protein's score against
/// the sequence cutoff, not any one alignment's, since a model can be reached by domains that
/// each fall short of it.
pub fn complete(table: &str, bars: &Bars) -> Hits {
    let mut best = Hits::new();
    for domain in parse(table) {
        let Some(bar) = bars.get(domain.model) else {
            continue;
        };
        if domain.sequence_score < bar.sequence {
            continue;
        }
        hmm_table::keep_best(
            &mut best,
            domain.protein,
            domain.model,
            domain.sequence_score,
        );
    }
    best
}

fn parse(table: &str) -> impl Iterator<Item = Domain<'_>> {
    hmm_table::rows(table).filter_map(|mut fields| {
        let protein = fields.at(hmm_table::DOMAIN_TARGET)?.parse::<usize>().ok()?;
        let model = fields.at(hmm_table::DOMAIN_MODEL)?;
        let length = fields
            .at(hmm_table::DOMAIN_MODEL_LENGTH)?
            .parse::<f64>()
            .ok()?;
        let sequence_score = fields
            .at(hmm_table::DOMAIN_SEQUENCE_SCORE)?
            .parse::<f64>()
            .ok()?;
        let score = fields.at(hmm_table::DOMAIN_SCORE)?.parse::<f64>().ok()?;
        let from = fields.at(hmm_table::DOMAIN_HMM_FROM)?.parse::<f64>().ok()?;
        let to = fields.at(hmm_table::DOMAIN_HMM_TO)?.parse::<f64>().ok()?;
        (length > 0.0).then(|| Domain {
            protein,
            model,
            sequence_score,
            score,
            span: ((to - from + 1.0) / length).clamp(0.0, 1.0),
        })
    })
}
