use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use log::{debug, warn};

use super::checkm::{self, Counted};
use super::{
    Annotator, ContigMarkers, HMM_GZ, Hit, MarkerSet, cache, fragments, in_marker_order, replicon,
};
use crate::external::hmmer_engine::{self, ProteinFile};
use crate::quality::orfs;

// HMMER refuses a target past this, and an ORF this long is an uncovered N span in a gold
// standard assembly or a scaffold gap, never a marker gene.
const MAX_SEARCH_RESIDUES: usize = 100_000;

type Sequenced = Vec<(Counted, String)>;

pub struct MarkerAnnotation {
    pub(super) rows: cache::Rows,
    set: MarkerSet,
    pub(super) sequences: Sequenced,
}

impl MarkerAnnotation {
    // `sequences` keeps each CheckM copy's protein, which only the strain reading needs.
    pub fn build(
        annotator: &Annotator,
        band: std::ops::Range<usize>,
        sequences: bool,
    ) -> Result<Self> {
        let set = MarkerSet::joined(annotator.panel.as_deref());
        let cached = annotator
            .cache
            .as_deref()
            .map(|directory| annotator.key(band.start).map(|key| (directory, key)))
            .transpose()?;
        let held = cached
            .as_ref()
            .and_then(|(directory, key)| cache::find(directory, key))
            .and_then(|(path, floor)| match cache::read(&path, &set) {
                Ok(rows) => {
                    debug!("Marker cache: read {}", path.display());
                    Some((floor, (path, rows)))
                }
                Err(error) => {
                    warn!("Marker cache: ignored {} ({error})", path.display());
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
                    sequences: Vec::new(),
                });
            }
            Some((floor, held)) if floor <= band.end => (floor, Some(held)),
            _ => (band.end, None),
        };
        let (mut rows, kept) = annotate(annotator, band.start..ceiling, &set, sequences)?;
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
                    debug!("Marker cache: wrote {}", path.display());
                    if let Some(old) = superseded.filter(|old| *old != path) {
                        let _ = std::fs::remove_file(old);
                    }
                }
                Err(error) => warn!("Marker cache: {error} ({})", path.display()),
            }
        }
        Ok(Self {
            rows,
            set,
            sequences: kept,
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
        let index = crate::rows::positions(self.rows.names.iter().map(String::as_str));
        let mut per_contig = Vec::with_capacity(names.len());
        let mut shapes = Vec::with_capacity(names.len());
        let mut checkm = Vec::with_capacity(names.len());
        for name in names {
            match index.get(name.as_str()) {
                Some(position) => {
                    per_contig.push(std::mem::take(&mut self.rows.hits[*position]));
                    shapes.push(Some(self.rows.shapes[*position]));
                    checkm.push(std::mem::take(&mut self.rows.checkm[*position]));
                }
                None if tolerate_missing => {
                    per_contig.push(Default::default());
                    shapes.push(None);
                    checkm.push(Vec::new());
                }
                None => anyhow::bail!("the marker table does not hold {name}"),
            }
        }
        Ok(ContigMarkers::new(per_contig, self.set)
            .with_shapes(shapes)
            .with_checkm(checkm))
    }
}

// Only contigs in `band` are called and searched. Those at or above its end come back unscored,
// for the caller to fill from a cached annotation.
fn annotate(
    annotator: &Annotator,
    band: std::ops::Range<usize>,
    set: &MarkerSet,
    sequences: bool,
) -> Result<(cache::Rows, Sequenced)> {
    let rules = annotator.rules;

    let directory = tempfile::tempdir()?;
    let hmm = directory.path().join("markers.hmm");
    inflate(HMM_GZ, &hmm)?;

    let (walked, called, proteins, shapes) = {
        let _timer = crate::timing::scope("genes");
        let mut sink = ProteinFile::create(directory.path().join("proteins.faa"))?;
        let mut called: Vec<orfs::Orf> = Vec::new();
        let mut shapes: Vec<replicon::Shape> = Vec::new();
        let mut written = annotator
            .proteins
            .as_deref()
            .map(GeneFasta::create)
            .transpose()?;
        let walked = orfs::call_over(&annotator.assembly, band.clone(), |batch, names| {
            for mut orf in batch {
                if let Some(written) = written.as_mut() {
                    written.write(&orf, names)?;
                }
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
        let proteins = sink.finish()?;
        if let Some(written) = written {
            written.finish()?;
        }
        shapes.resize(walked.names.len(), replicon::Shape::default());
        let banded = walked.lengths.iter().filter(|length| **length < band.end);
        debug!("Genes: {} over {} contigs", called.len(), banded.count());
        (walked, called, proteins, shapes)
    };

    let bars = fragments::gathering(&hmm)?;
    let table = {
        let _timer = crate::timing::scope("search");
        // CheckM reads this table for the models both panels share, so it reports down to
        // CheckM's lowest bar as well.
        let floor = fragments::floor(&bars, rules.fragment_span).min(set.checkm.floor());
        match &proteins {
            Some(proteins) => hmmer_engine::search(&hmm, proteins, floor)?,
            None => String::new(),
        }
    };
    let hits = {
        let _timer = crate::timing::scope("fragments");
        let cut = |protein: usize| called.get(protein).is_some_and(orfs::Orf::partial);
        let (hits, rescued) = fragments::called(&table, &bars, rules.fragment_span, cut);
        debug!("Markers from cut genes: {rescued}");
        hits
    };

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
    let carriers = per_contig.iter().filter(|hits| !hits.is_empty()).count();
    debug!("Marker contigs: {carriers} of {}", walked.names.len());
    if let (Some(panel), Some(proteins)) = (&annotator.panel, &proteins) {
        let _timer = crate::timing::scope("eukaryote panels");
        let found = panel.search(proteins, &called, set, rules.fragment_span)?;
        debug!("Eukaryote panels: {} hits", found.len());
        for (contig, hit) in found {
            per_contig[contig].push(hit);
        }
    }
    in_marker_order(&mut per_contig);
    let (checkm, kept) = match &proteins {
        Some(proteins) => {
            let _timer = crate::timing::scope("checkm");
            search_checkm(set, directory.path(), proteins, &table, &called, sequences)?
        }
        None => (Vec::new(), Vec::new()),
    };
    let mut checkm = set.checkm.tally(&checkm, walked.names.len());
    checkm.resize(walked.names.len(), Vec::new());
    Ok((
        cache::Rows {
            checkm,
            names: walked.names,
            lengths: walked.lengths,
            hits: per_contig,
            shapes,
        },
        kept,
    ))
}

fn search_checkm(
    set: &MarkerSet,
    directory: &Path,
    proteins: &Path,
    table: &str,
    called: &[orfs::Orf],
    sequences: bool,
) -> Result<(Vec<Counted>, Sequenced)> {
    let hmm = directory.join("checkm.hmm");
    inflate(checkm::HMM_GZ, &hmm)?;
    let panel = &set.checkm;
    let joined = table.to_string() + &hmmer_engine::search(&hmm, proteins, panel.floor())?;
    let counted = panel.counted(&joined, |protein| {
        called.get(protein).map(|orf| (orf.contig, orf.partial()))
    });
    if !sequences {
        return Ok((counted, Vec::new()));
    }
    let wanted = counted
        .iter()
        .flat_map(|copy| [Some(copy.protein), copy.partner])
        .flatten()
        .collect::<HashSet<_>>();
    let mut held = HashMap::new();
    hmmer_engine::each_protein(proteins, |id, protein| {
        if wanted.contains(&id) {
            held.insert(id, protein.to_string());
        }
        Ok(())
    })?;
    let kept = counted
        .iter()
        .map(|copy| {
            let joined = [Some(copy.protein), copy.partner]
                .into_iter()
                .flatten()
                .filter_map(|id| held.get(&id).map(String::as_str))
                .collect::<String>();
            (*copy, joined)
        })
        .collect();
    Ok((counted, kept))
}

fn searchable(protein: &str) -> bool {
    !protein.is_empty() && protein.len() <= MAX_SEARCH_RESIDUES
}

struct GeneFasta {
    sink: BufWriter<std::fs::File>,
    contig: usize,
    gene: usize,
}

impl GeneFasta {
    fn create(path: &Path) -> Result<Self> {
        Ok(Self {
            sink: BufWriter::new(std::fs::File::create(path)?),
            contig: usize::MAX,
            gene: 0,
        })
    }

    fn write(&mut self, orf: &orfs::Orf, names: &[String]) -> Result<()> {
        if orf.contig != self.contig {
            (self.contig, self.gene) = (orf.contig, 0);
        }
        self.gene += 1;
        // A refiner's HMMER search aborts on the whole run at one target this long.
        if searchable(&orf.protein) {
            writeln!(
                self.sink,
                ">{}_{}\n{}",
                names[orf.contig], self.gene, orf.protein
            )?;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        self.sink.flush()?;
        Ok(())
    }
}

pub(super) fn inflate(compressed: &[u8], target: &Path) -> Result<()> {
    let mut decoder = flate2::read::GzDecoder::new(compressed);
    let mut sink = BufWriter::new(std::fs::File::create(target)?);
    std::io::copy(&mut decoder, &mut sink)?;
    sink.flush()?;
    Ok(())
}
