use std::collections::HashSet;
use std::io::{BufWriter, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::Result;
use needletail::parser::{LineEnding, write_fasta};

use super::checkm::Counted;
use super::panels::Panel;
use super::{ContigMarkers, MarkerAnnotation, MarkerRules, cache};

#[derive(Clone)]
pub struct Annotator {
    pub assembly: String,
    pub rules: MarkerRules,
    pub cache: Option<PathBuf>,
    pub content: Option<u64>,
    pub proteins: Option<PathBuf>,
    pub panel: Option<Arc<Panel>>,
    pub checkm_tree: Option<Arc<super::genome_tree::Data>>,
    pub elements: Option<Arc<super::elements::Elements>>,
}

impl Annotator {
    pub(super) fn key(&self, floor: usize) -> Result<String> {
        cache::key(
            &self.assembly,
            self.content,
            floor,
            self.rules.fragment_span,
            self.panel.as_deref().map(Panel::fingerprint),
        )
    }

    pub fn annotate(&self, band: Range<usize>) -> Result<MarkerAnnotation> {
        MarkerAnnotation::build(self, band, false)
    }

    pub fn annotate_named(&self, names: &[String], file: &str) -> Result<MarkerAnnotation> {
        self.named(names, file, false)
    }

    pub(super) fn checkm_proteins(
        &self,
        names: &[String],
        file: &str,
    ) -> Result<Vec<(String, Counted, String)>> {
        let MarkerAnnotation {
            rows, sequences, ..
        } = self.named(names, file, true)?;
        Ok(sequences
            .into_iter()
            .map(|(copy, protein)| (rows.names[copy.contig].clone(), copy, protein))
            .collect())
    }

    fn named(&self, names: &[String], file: &str, sequences: bool) -> Result<MarkerAnnotation> {
        let directory = tempfile::tempdir()?;
        let subset = directory.path().join(file);
        write_contigs(&self.assembly, names, &subset)?;
        let annotator = Self {
            assembly: subset.to_string_lossy().into_owned(),
            cache: None,
            proteins: None,
            ..self.clone()
        };
        MarkerAnnotation::build(&annotator, 0..usize::MAX, sequences)
    }
}

pub(super) fn write_contigs(assembly: &str, names: &[String], path: &Path) -> Result<()> {
    let wanted = names.iter().map(String::as_str).collect::<HashSet<_>>();
    let mut sink = BufWriter::new(std::fs::File::create(path)?);
    crate::kmers::each_named(
        assembly,
        |name| wanted.contains(name),
        |name, record| {
            write_fasta(name.as_bytes(), &record.seq(), &mut sink, LineEnding::Unix)?;
            Ok(())
        },
    )?;
    sink.flush()?;
    Ok(())
}

impl ContigMarkers {
    pub fn fill(&mut self, annotation: MarkerAnnotation, names: &[String], contigs: &[usize]) {
        let mut rows = annotation.rows;
        let index = crate::rows::positions(rows.names.iter().map(String::as_str));
        let found = contigs
            .iter()
            .filter_map(|contig| Some((*contig, *index.get(names[*contig].as_str())?)))
            .collect::<Vec<_>>();
        for (contig, at) in found {
            self.per_contig[contig] = std::mem::take(&mut rows.hits[at]);
            if let Some(shape) = self.shapes.get_mut(contig) {
                *shape = rows.shapes[at];
                self.annotated[contig] = true;
            }
            if let Some(copies) = self.checkm.get_mut(contig) {
                *copies = std::mem::take(&mut rows.checkm[at]);
            }
        }
    }
}
