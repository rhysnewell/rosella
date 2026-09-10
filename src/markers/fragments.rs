use std::collections::HashMap;
use std::path::Path;

use anyhow::Result;

pub const DEFAULT_SPAN: f64 = 0.3;
pub const DOMAIN_FLOOR: &str = "10";

pub fn gathering(hmm: &Path) -> Result<HashMap<String, f64>> {
    let text = std::fs::read_to_string(hmm)?;
    let mut bars = HashMap::new();
    let mut name: Option<String> = None;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("NAME") {
            name = Some(rest.trim().to_string());
        } else if let Some(rest) = line.strip_prefix("GA") {
            let Some(model) = name.take() else {
                continue;
            };
            if let Some(bar) = rest
                .trim()
                .trim_end_matches(';')
                .split_whitespace()
                .next()
                .and_then(|value| value.parse::<f64>().ok())
            {
                bars.insert(model, bar);
            }
        }
    }
    Ok(bars)
}

pub struct Domain {
    pub protein: String,
    pub model: String,
    pub score: f64,
    pub span: f64,
}

/// A gene cut by a contig end can only align to the part of the model it still carries, so the
/// full length gathering threshold is scaled to the span that could have matched at all.
pub fn accepted(
    table: &str,
    bars: &HashMap<String, f64>,
    min_span: f64,
) -> HashMap<String, (String, f64)> {
    let mut best: HashMap<String, (String, f64)> = HashMap::new();
    for domain in parse(table) {
        let Some(bar) = bars.get(&domain.model) else {
            continue;
        };
        if domain.span < min_span || domain.score < bar * domain.span {
            continue;
        }
        match best.get(&domain.protein) {
            Some((_, held)) if *held >= domain.score => {}
            _ => {
                best.insert(domain.protein, (domain.model, domain.score));
            }
        }
    }
    best
}

fn parse(table: &str) -> impl Iterator<Item = Domain> + '_ {
    table
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            let fields = line.split_whitespace().collect::<Vec<_>>();
            let protein = (*fields.first()?).to_string();
            let model = (*fields.get(3)?).to_string();
            let length = fields.get(5)?.parse::<f64>().ok()?;
            let score = fields.get(13)?.parse::<f64>().ok()?;
            let from = fields.get(15)?.parse::<f64>().ok()?;
            let to = fields.get(16)?.parse::<f64>().ok()?;
            (length > 0.0).then(|| Domain {
                protein,
                model,
                score,
                span: ((to - from + 1.0) / length).clamp(0.0, 1.0),
            })
        })
}
