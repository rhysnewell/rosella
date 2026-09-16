use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use log::{debug, info, warn};

use crate::external::hmmer_engine::HmmerEngine;
use crate::quality::{Quality, orfs};

pub mod cache;
pub mod fragments;
pub mod hmm_table;
pub mod sets;
mod table;

pub use table::MarkerSet;

pub(crate) const HMM_GZ: &[u8] = include_bytes!("../../data/gtdb_markers.hmm.gz");

pub const DEFAULT_BAR_OFFSET: f64 = 10.0;

// hmmsearch refuses a target past this, and an ORF this long is an uncovered N span in a gold
// standard assembly or a scaffold gap, never a marker gene.
const MAX_SEARCH_RESIDUES: usize = 100_000;

/// A marker gene cut by a contig end is still that marker gene, but two halves of one gene on
/// two contigs are not two copies, so presence and duplication read different columns.
#[derive(Clone, Copy, Debug)]
pub struct MarkerRules {
    pub fragment_span: f64,
    pub bar_offset: f64,
}

impl Default for MarkerRules {
    fn default() -> Self {
        Self {
            fragment_span: fragments::DEFAULT_SPAN,
            bar_offset: DEFAULT_BAR_OFFSET,
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Tally {
    complete: u32,
    any: u32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Hit {
    pub marker: u16,
    pub partial: bool,
}


pub struct MarkerAnnotation {
    names: Vec<String>,
    per_contig: Vec<Vec<Hit>>,
    set: MarkerSet,
    rules: MarkerRules,
}

impl MarkerAnnotation {
    pub fn build(
        assembly: &str,
        min_contig_size: usize,
        threads: usize,
        shards: Option<usize>,
        rules: MarkerRules,
        cache: Option<&Path>,
    ) -> Result<Self> {
        let set = MarkerSet::embedded();
        let cached = cache
            .map(|directory| {
                cache::key(assembly, min_contig_size, rules.fragment_span)
                    .map(|key| (directory, key))
            })
            .transpose()?;
        if let Some(path) = cached
            .as_ref()
            .and_then(|(directory, key)| cache::find(directory, key))
        {
            match cache::read(&path, &set) {
                Ok((names, per_contig)) => {
                    info!("Read the marker annotation from {}", path.display());
                    return Ok(Self {
                        names,
                        per_contig,
                        set,
                        rules,
                    });
                }
                Err(error) => warn!("Ignoring {}: {error}", path.display()),
            }
        }
        // Annotating is most of a run and the temp directory is dropped on any failure, so the
        // search has to be known to work before the gene calling is paid for.
        HmmerEngine::check_installed()?;
        let engine = HmmerEngine::new(threads, shards);

        let directory = tempfile::tempdir()?;
        let hmm = directory.path().join("markers.hmm");
        inflate(HMM_GZ, &hmm)?;

        let (names, called, pieces) = {
            let _timer = crate::timing::scope("genes");
            let mut sink = engine.protein_shards(directory.path())?;
            let mut called: Vec<orfs::Orf> = Vec::new();
            let names = orfs::call_over(assembly, min_contig_size, |batch| {
                for mut orf in batch {
                    if searchable(&orf.protein) {
                        sink.write(called.len(), &orf.protein)?;
                    }
                    if !orf.partial {
                        orf.protein = String::new();
                    }
                    called.push(orf);
                }
                Ok(())
            })?;
            let pieces = sink.finish()?;
            info!("Called {} genes over {} contigs", called.len(), names.len());
            (names, called, pieces)
        };

        let bars = fragments::gathering(&hmm)?;
        let table = {
            let _timer = crate::timing::scope("search");
            let floor = fragments::floor(&bars, rules.fragment_span);
            engine.search(&hmm, &pieces, directory.path(), &floor)?
        };
        let mut hits = fragments::complete(&table, &bars);
        {
            let _timer = crate::timing::scope("fragments");
            let cut = |protein: usize| {
                called.get(protein).is_some_and(|orf| orf.partial) && !hits.contains_key(&protein)
            };
            let rescued = fragments::accepted(&table, &bars, rules.fragment_span, cut);
            debug!("{} markers rescued from cut genes", rescued.len());
            hits.extend(rescued);
        }

        let mut per_contig = vec![Vec::new(); names.len()];
        for (protein, (model, _)) in hits {
            let Some(marker) = set.id(&model) else {
                continue;
            };
            let Some(orf) = called.get(protein) else {
                continue;
            };
            per_contig[orf.contig].push(Hit {
                marker,
                partial: orf.partial,
            });
        }
        in_marker_order(&mut per_contig);
        let carriers = per_contig.iter().filter(|hits| !hits.is_empty()).count();
        debug!(
            "{carriers} of {} contigs carry a single copy marker",
            names.len()
        );
        if let Some((directory, key)) = cached.as_ref() {
            let path = cache::write_path(directory, key);
            match cache::write(&path, key, &set, &names, &per_contig) {
                Ok(()) => info!("Wrote the marker annotation to {}", path.display()),
                Err(error) => warn!("Could not write {}: {error}", path.display()),
            }
        }
        Ok(Self {
            names,
            per_contig,
            set,
            rules,
        })
    }

    pub fn report(&self, path: &Path) -> Result<()> {
        let mut sink = BufWriter::new(std::fs::File::create(path)?);
        writeln!(sink, "contig\tmodel\tpartial")?;
        for (contig, hits) in self.names.iter().zip(&self.per_contig) {
            for hit in hits {
                writeln!(
                    sink,
                    "{contig}\t{}\t{}",
                    self.set.name(hit.marker),
                    u8::from(hit.partial)
                )?;
            }
        }
        sink.flush()?;
        Ok(())
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Foreign bins may hold contigs this assembly never annotated, and a missing contig is a
    /// bin with no features rather than a reason to refuse the whole table.
    pub fn select_present(self, names: &[String]) -> Result<ContigMarkers> {
        self.selected(names, true)
    }

    pub fn select(self, names: &[String]) -> Result<ContigMarkers> {
        self.selected(names, false)
    }

    fn selected(self, names: &[String], tolerate_missing: bool) -> Result<ContigMarkers> {
        let index = self
            .names
            .iter()
            .enumerate()
            .map(|(position, name)| (name.as_str(), position))
            .collect::<HashMap<_, _>>();
        let mut per_contig = Vec::with_capacity(names.len());
        for name in names {
            match index.get(name.as_str()) {
                Some(position) => per_contig.push(self.per_contig[*position].clone()),
                None if tolerate_missing => per_contig.push(Default::default()),
                None => anyhow::bail!("the marker table does not hold {name}"),
            }
        }
        Ok(ContigMarkers::new(per_contig, self.set, self.rules))
    }
}

/// The search returns its hits in hash order, and the cache and the report are written from
/// these, so two annotations of one assembly would not diff against each other.
fn in_marker_order(per_contig: &mut [Vec<Hit>]) {
    for hits in per_contig {
        hits.sort_unstable_by_key(|hit| (hit.marker, hit.partial));
    }
}

pub struct ContigMarkers {
    per_contig: Vec<Vec<Hit>>,
    lengths: Vec<usize>,
    set: MarkerSet,
    rules: MarkerRules,
}

impl crate::quality::Scorer for ContigMarkers {
    fn features(&self, contigs: &[usize]) -> std::collections::HashSet<u32> {
        self.counts(contigs)
            .iter()
            .enumerate()
            .filter(|(_, tally)| tally.any > 0)
            .map(|(marker, _)| marker as u32)
            .collect()
    }

    /// Measured on bins that are whole genomes: the same number is a harsher test here than
    /// for a model that predicts the share of a genome, and this offset matches the two.
    fn completeness_bar(&self, requested: f64) -> f64 {
        (requested - self.rules.bar_offset).max(0.0)
    }

    fn set_name(&self, set: u16) -> &str {
        self.set.sets.name(set as usize)
    }

    fn smallest_scale(&self) -> f64 {
        self.set.sets.smallest_scale()
    }

    /// Read against whichever lineage the bin's pattern of absences fits, since a reduced
    /// genome is missing markers a whole one of another lineage would carry.
    fn score(&self, contigs: &[usize]) -> Quality {
        let counts = self.counts(contigs);
        let Some(chosen) = self
            .set
            .sets
            .choose(&observed(&counts), self.bin_bp(contigs))
        else {
            return Quality::default();
        };
        let (mut present, mut extra, mut total) = (0usize, 0usize, 0usize);
        for (marker, tally) in counts.iter().enumerate() {
            if !self.set.sets.holds(chosen, marker) {
                continue;
            }
            total += 1;
            present += usize::from(tally.any >= 1);
            extra += tally.complete.saturating_sub(1) as usize;
        }
        if total == 0 {
            return Quality::default();
        }
        Quality {
            completeness: 100.0 * present as f64 / total as f64,
            contamination: 100.0 * extra as f64 / total as f64,
            scale: self.set.sets.scale(chosen),
            set: chosen as u16,
        }
    }
}

fn observed(counts: &[Tally]) -> Vec<u16> {
    counts
        .iter()
        .enumerate()
        .filter(|(_, tally)| tally.any > 0)
        .map(|(marker, _)| marker as u16)
        .collect()
}

impl ContigMarkers {
    pub fn new(mut per_contig: Vec<Vec<Hit>>, set: MarkerSet, rules: MarkerRules) -> Self {
        in_marker_order(&mut per_contig);
        Self {
            per_contig,
            lengths: Vec::new(),
            set,
            rules,
        }
    }

    pub fn with_lengths(mut self, lengths: Vec<usize>) -> Self {
        self.lengths = lengths;
        self
    }

    fn bin_bp(&self, contigs: &[usize]) -> usize {
        contigs
            .iter()
            .filter_map(|contig| self.lengths.get(*contig))
            .sum()
    }

    fn counts(&self, contigs: &[usize]) -> Vec<Tally> {
        let mut counts = vec![Tally::default(); self.set.len()];
        for contig in contigs {
            for hit in &self.per_contig[*contig] {
                let tally = &mut counts[hit.marker as usize];
                tally.any += 1;
                if !hit.partial {
                    tally.complete += 1;
                }
            }
        }
        counts
    }
}

fn searchable(protein: &str) -> bool {
    !protein.is_empty() && protein.len() <= MAX_SEARCH_RESIDUES
}

fn inflate(compressed: &[u8], target: &Path) -> Result<()> {
    let mut decoder = flate2::read::GzDecoder::new(compressed);
    let mut sink = BufWriter::new(std::fs::File::create(target)?);
    std::io::copy(&mut decoder, &mut sink)?;
    sink.flush()?;
    Ok(())
}

