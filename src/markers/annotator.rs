use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};

use anyhow::Result;
use log::info;

use super::{ContigMarkers, MarkerAnnotation, MarkerRules, cache};

#[derive(Clone)]
pub struct Annotator {
    pub assembly: String,
    pub threads: usize,
    pub shards: Option<usize>,
    pub rules: MarkerRules,
    pub cache: Option<PathBuf>,
    pub checkm: bool,
}

impl Annotator {
    pub fn annotate(&self, band: Range<usize>) -> Result<MarkerAnnotation> {
        MarkerAnnotation::build(self, band)
    }

    pub fn complete_checkm(
        &self,
        markers: &mut ContigMarkers,
        contigs: &[usize],
        names: &[String],
    ) -> Result<()> {
        if !self.checkm {
            return Ok(());
        }
        // A contig outside these bins reads NA even when the cache knows it from another run.
        let mut binned = vec![false; markers.checkm.len()];
        for contig in contigs {
            binned[*contig] = true;
        }
        for (copies, binned) in markers.checkm.iter_mut().zip(binned) {
            copies.take_if(|_| !binned);
        }
        let mut searched = markers.search_checkm(contigs)?;
        let mut missing = contigs
            .iter()
            .copied()
            .filter(|contig| markers.checkm[*contig].is_none())
            .collect::<Vec<_>>();
        missing.sort_unstable();
        missing.dedup();
        if !missing.is_empty() {
            info!(
                "Calling genes again on {} binned contigs no CheckM search has seen.",
                missing.len()
            );
            let directory = tempfile::tempdir()?;
            let subset = directory.path().join("unsearched.fna");
            let wanted = missing
                .iter()
                .map(|contig| names[*contig].clone())
                .collect::<Vec<_>>();
            write_contigs(&self.assembly, &wanted, &subset)?;
            let annotator = Self {
                assembly: subset.to_string_lossy().into_owned(),
                cache: None,
                ..self.clone()
            };
            let mut extra = annotator.annotate(0..usize::MAX)?.select(&wanted)?;
            extra.search_checkm(&(0..wanted.len()).collect::<Vec<_>>())?;
            for (at, contig) in missing.iter().enumerate() {
                markers.checkm[*contig] = extra.checkm[at].take();
            }
            searched.extend(missing);
        }
        let Some(directory) = &self.cache else {
            return Ok(());
        };
        let Some(shortest) = searched
            .iter()
            .filter_map(|contig| markers.lengths.get(*contig))
            .min()
        else {
            return Ok(());
        };
        let found = searched
            .iter()
            .filter_map(|contig| {
                let copies = markers.checkm[*contig].as_deref()?;
                Some((names[*contig].as_str(), copies))
            })
            .collect::<HashMap<_, _>>();
        let key = cache::key(&self.assembly, *shortest, self.rules.fragment_span)?;
        cache::record(directory, &key, &markers.set, &found)
    }
}

fn write_contigs(assembly: &str, names: &[String], path: &Path) -> Result<()> {
    let wanted = names.iter().map(String::as_str).collect::<HashSet<_>>();
    let mut reader = needletail::parse_fastx_file(assembly)?;
    let mut sink = BufWriter::new(std::fs::File::create(path)?);
    while let Some(record) = reader.next() {
        let record = record?;
        let name = crate::contig_id(record.id())?;
        if wanted.contains(name) {
            writeln!(sink, ">{name}")?;
            sink.write_all(&record.seq())?;
            writeln!(sink)?;
        }
    }
    sink.flush()?;
    Ok(())
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
        let mut placed = vec![None; rows.names.len()];
        for (contig, at) in found {
            placed[at] = Some(contig);
            self.per_contig[contig] = std::mem::take(&mut rows.hits[at]);
            if let Some(shape) = self.shapes.get_mut(contig) {
                *shape = rows.shapes[at];
            }
            if let Some(copies) = self.checkm.get_mut(contig) {
                *copies = rows.checkm[at].take();
            }
        }
        self.proteins
            .extend(annotation.proteins.map(|kept| kept.renumbered(&placed)));
    }
}
