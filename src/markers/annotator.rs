use std::collections::HashMap;
use std::ops::Range;
use std::path::PathBuf;

use anyhow::Result;

use super::{ContigMarkers, MarkerAnnotation, MarkerRules};

pub struct Annotator {
    pub assembly: String,
    pub threads: usize,
    pub shards: Option<usize>,
    pub rules: MarkerRules,
    pub cache: Option<PathBuf>,
}

impl Annotator {
    pub fn annotate(&self, band: Range<usize>) -> Result<MarkerAnnotation> {
        MarkerAnnotation::build(
            &self.assembly,
            band,
            self.threads,
            self.shards,
            self.rules,
            self.cache.as_deref(),
        )
    }
}

impl ContigMarkers {
    pub fn fill(&mut self, annotation: MarkerAnnotation, names: &[String], contigs: &[usize]) {
        let mut rows = annotation.rows;
        let index = rows
            .names
            .iter()
            .enumerate()
            .map(|(at, name)| (name.as_str(), at))
            .collect::<HashMap<_, _>>();
        let found = contigs
            .iter()
            .filter_map(|contig| Some((*contig, *index.get(names[*contig].as_str())?)))
            .collect::<Vec<_>>();
        for (contig, at) in found {
            self.per_contig[contig] = std::mem::take(&mut rows.hits[at]);
            if let Some(shape) = self.shapes.get_mut(contig) {
                *shape = rows.shapes[at];
            }
        }
    }
}
