use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use log::{debug, warn};

use super::{
    Annotator, ContigMarkers, HMM_GZ, Hit, MarkerSet, cache, fragments, in_marker_order, proteins,
    replicon,
};
use crate::external::hmmer_engine::HmmerEngine;
use crate::quality::orfs;

// hmmsearch refuses a target past this, and an ORF this long is an uncovered N span in a gold
// standard assembly or a scaffold gap, never a marker gene.
const MAX_SEARCH_RESIDUES: usize = 100_000;

pub struct MarkerAnnotation {
    pub(super) rows: cache::Rows,
    set: MarkerSet,
    pub(super) proteins: Option<proteins::Proteins>,
    checkm: bool,
}

impl MarkerAnnotation {
    pub fn build(annotator: &Annotator, band: std::ops::Range<usize>) -> Result<Self> {
        let (assembly, rules) = (annotator.assembly.as_str(), annotator.rules);
        let set = MarkerSet::embedded();
        let cached = annotator
            .cache
            .as_deref()
            .map(|directory| {
                cache::key(assembly, annotator.content, band.start, rules.fragment_span)
                    .map(|key| (directory, key))
            })
            .transpose()?;
        let held = cached
            .as_ref()
            .and_then(|(directory, key)| cache::find(directory, key))
            .and_then(|(path, floor)| match cache::read(&path, &set) {
                Ok(rows) => {
                    debug!("Read the marker annotation from {}", path.display());
                    Some((floor, (path, rows)))
                }
                Err(error) => {
                    warn!("Ignoring {}: {error}", path.display());
                    None
                }
            });
        let (ceiling, held) = match held {
            Some((floor, (path, rows))) if floor <= band.start => {
                if let Some((directory, key)) = &cached {
                    cache::upgrade(directory, key, &path, floor, &set, &rows);
                }
                return Ok(Self {
                    rows: rows.at_least(band.start),
                    set,
                    proteins: None,
                    checkm: annotator.checkm,
                });
            }
            Some((floor, held)) if floor <= band.end => (floor, Some(held)),
            _ => (band.end, None),
        };
        let (mut rows, proteins) = annotate(annotator, band.start..ceiling, &set)?;
        let whole = held.is_some() || band.end == usize::MAX;
        let superseded = match held {
            Some((path, held)) => {
                rows = rows.fill_from(held, ceiling)?;
                Some(path)
            }
            None => None,
        };
        if let Some((directory, key)) = cached.as_ref().filter(|_| whole) {
            let path = cache::write_path(directory, key);
            match cache::write(&path, key, &set, &rows) {
                Ok(()) => {
                    debug!("Wrote the marker annotation to {}", path.display());
                    if let Some(old) = superseded.filter(|old| *old != path) {
                        let _ = std::fs::remove_file(old);
                    }
                }
                Err(error) => warn!("Could not write {}: {error}", path.display()),
            }
        }
        Ok(Self {
            rows,
            set,
            proteins,
            checkm: annotator.checkm,
        })
    }

    // Contigs under the annotated band are filled later or not at all, so a missing contig
    // reads as featureless rather than refusing the whole table.
    pub fn select_present(self, names: &[String]) -> Result<ContigMarkers> {
        self.selected(names, true)
    }

    pub fn select(self, names: &[String]) -> Result<ContigMarkers> {
        self.selected(names, false)
    }

    fn selected(mut self, names: &[String], tolerate_missing: bool) -> Result<ContigMarkers> {
        let index = self
            .rows
            .names
            .iter()
            .enumerate()
            .map(|(position, name)| (name.as_str(), position))
            .collect::<HashMap<_, _>>();
        let mut per_contig = Vec::with_capacity(names.len());
        let mut shapes = Vec::with_capacity(names.len());
        let mut checkm = Vec::with_capacity(names.len());
        let mut placed = vec![None; self.rows.names.len()];
        for (contig, name) in names.iter().enumerate() {
            match index.get(name.as_str()) {
                Some(position) => {
                    placed[*position] = Some(contig);
                    per_contig.push(std::mem::take(&mut self.rows.hits[*position]));
                    shapes.push(Some(self.rows.shapes[*position]));
                    checkm.push(self.rows.checkm[*position].take());
                }
                None if tolerate_missing => {
                    per_contig.push(Default::default());
                    shapes.push(None);
                    checkm.push(None);
                }
                None => anyhow::bail!("the marker table does not hold {name}"),
            }
        }
        let markers = ContigMarkers::new(per_contig, self.set).with_shapes(shapes);
        Ok(match self.checkm {
            true => markers
                .with_checkm(checkm)
                .with_proteins(self.proteins.map(|kept| kept.renumbered(&placed))),
            false => markers,
        })
    }
}

// Only contigs in `band` are called and searched. Those at or above its end come back unscored,
// for the caller to fill from a cached annotation.
fn annotate(
    annotator: &Annotator,
    band: std::ops::Range<usize>,
    set: &MarkerSet,
) -> Result<(cache::Rows, Option<proteins::Proteins>)> {
    let engine = HmmerEngine::new(annotator.threads, annotator.shards);
    let rules = annotator.rules;

    let directory = tempfile::tempdir()?;
    let hmm = directory.path().join("markers.hmm");
    inflate(HMM_GZ, &hmm)?;

    let (walked, called, pieces, shapes) = {
        let _timer = crate::timing::scope("genes");
        let mut sink = engine.protein_shards(directory.path())?;
        let mut called: Vec<orfs::Orf> = Vec::new();
        let mut shapes: Vec<replicon::Shape> = Vec::new();
        let walked = orfs::call_over(&annotator.assembly, band.clone(), |batch| {
            for mut orf in batch {
                if shapes.len() <= orf.contig {
                    shapes.resize(orf.contig + 1, replicon::Shape::default());
                }
                shapes[orf.contig].add(orf.bases);
                if searchable(&orf.protein) {
                    sink.write(called.len(), &orf.protein)?;
                }
                if !orf.partial() {
                    orf.protein = String::new();
                }
                called.push(orf);
            }
            Ok(())
        })?;
        let pieces = sink.finish()?;
        shapes.resize(walked.names.len(), replicon::Shape::default());
        let banded = walked.lengths.iter().filter(|length| **length < band.end);
        debug!(
            "Called {} genes over {} contigs",
            called.len(),
            banded.count()
        );
        (walked, called, pieces, shapes)
    };

    let bars = fragments::gathering(&hmm)?;
    let table = {
        let _timer = crate::timing::scope("search");
        // CheckM reads this table for the models both panels share, so it reports down to
        // CheckM's lowest bar as well.
        let floor = fragments::floor(&bars, rules.fragment_span).min(set.checkm.floor());
        engine.search(&hmm, &pieces, directory.path(), floor, "markers")?
    };
    let mut hits = fragments::complete(&table, &bars);
    {
        let _timer = crate::timing::scope("fragments");
        let cut = |protein: usize| {
            called.get(protein).is_some_and(|orf| orf.partial()) && !hits.contains_key(&protein)
        };
        let rescued = fragments::accepted(&table, &bars, rules.fragment_span, cut);
        debug!("{} markers rescued from cut genes", rescued.len());
        hits.extend(rescued);
    }

    let mut per_contig = vec![Vec::new(); walked.names.len()];
    for (protein, best) in hits {
        let Some(marker) = set.id(&best.model) else {
            continue;
        };
        let Some(orf) = called.get(protein) else {
            continue;
        };
        per_contig[orf.contig].push(Hit::called(marker, orf, &best));
    }
    in_marker_order(&mut per_contig);
    let carriers = per_contig.iter().filter(|hits| !hits.is_empty()).count();
    debug!(
        "{carriers} of {} contigs carry a single copy marker",
        walked.names.len()
    );
    let proteins = annotator
        .checkm
        .then(|| {
            let covered = walked.lengths.iter().map(|length| *length < band.end);
            proteins::Proteins::keep(directory, pieces, &table, &called, covered, engine)
        })
        .transpose()?;
    Ok((
        cache::Rows {
            checkm: vec![None; walked.names.len()],
            names: walked.names,
            lengths: walked.lengths,
            hits: per_contig,
            shapes,
        },
        proteins,
    ))
}

fn searchable(protein: &str) -> bool {
    !protein.is_empty() && protein.len() <= MAX_SEARCH_RESIDUES
}

pub(super) fn inflate(compressed: &[u8], target: &Path) -> Result<()> {
    let mut decoder = flate2::read::GzDecoder::new(compressed);
    let mut sink = BufWriter::new(std::fs::File::create(target)?);
    std::io::copy(&mut decoder, &mut sink)?;
    sink.flush()?;
    Ok(())
}
