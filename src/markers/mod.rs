use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use log::{debug, info, warn};

use crate::external::hmmer_engine::HmmerEngine;
use crate::quality::{Quality, orfs};

pub mod cache;
pub mod fragments;
pub mod sets;

pub(crate) const HMM_GZ: &[u8] = include_bytes!("../../data/gtdb_markers.hmm.gz");
const TABLE: &str = include_str!("../../data/gtdb_markers.tsv");
const SET_TABLE: &str = include_str!("../../data/marker_sets.tsv");

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

pub struct MarkerSet {
    ids: HashMap<String, u16>,
    names: Vec<String>,
    sets: sets::Sets,
}

const FALLBACK_SETS: [(&str, &str); 2] = [("bac", "bac120"), ("ar", "ar53")];

impl MarkerSet {
    pub fn embedded() -> Self {
        Self::parse(TABLE).with_scales(SET_TABLE)
    }

    pub fn parse(table: &str) -> Self {
        let mut lines = table.lines();
        let header = lines
            .next()
            .unwrap_or_default()
            .split('\t')
            .collect::<Vec<_>>();
        let column = |name: &str| header.iter().position(|field| *field == name);
        let Some(name_at) = column("model_name") else {
            return Self {
                ids: HashMap::new(),
                names: Vec::new(),
                sets: sets::Sets::default(),
            };
        };
        let set_at = column("sets");
        let domain_at = column("domain");
        let group_names = match set_at {
            Some(_) => set_names(table),
            None => FALLBACK_SETS
                .iter()
                .map(|(set, _)| (*set).to_string())
                .collect(),
        };
        let rate_at = group_names
            .iter()
            .map(|group| column(&format!("ubiquity_{group}")))
            .collect::<Vec<_>>();

        let mut ids = HashMap::new();
        let mut names = Vec::new();
        let mut member_of = vec![Vec::new(); group_names.len()];
        let mut rates = vec![Vec::new(); group_names.len()];
        for line in lines {
            let fields = line.split('\t').collect::<Vec<_>>();
            let Some(name) = fields.get(name_at) else {
                continue;
            };
            if ids.contains_key(*name) {
                continue;
            }
            ids.insert((*name).to_string(), names.len() as u16);
            names.push((*name).to_string());
            for (group, held) in group_names.iter().enumerate() {
                let member = match set_at.and_then(|at| fields.get(at)) {
                    Some(listed) => listed.split(',').any(|entry| entry == held),
                    None => domain_at
                        .and_then(|at| fields.get(at))
                        .is_some_and(|domain| domain.contains(FALLBACK_SETS[group].1)),
                };
                member_of[group].push(member);
                let rate = rate_at[group]
                    .and_then(|at| fields.get(at))
                    .and_then(|field| field.parse::<f64>().ok())
                    .unwrap_or(f64::from(u8::from(member)));
                rates[group].push(rate);
            }
        }
        Self {
            ids,
            names,
            sets: sets::Sets::new(group_names, member_of, rates),
        }
    }

    pub fn with_scales(mut self, table: &str) -> Self {
        let mut expected = vec![0.0; self.sets.len()];
        let mut bounds = vec![0.0; self.sets.len()];
        let mut lines = table.lines();
        let header = lines
            .next()
            .unwrap_or_default()
            .split('\t')
            .collect::<Vec<_>>();
        let column = |name: &str| header.iter().position(|field| *field == name);
        let (Some(set_at), Some(bp_at)) = (column("set"), column("median_genome_bp")) else {
            return self;
        };
        let max_at = column("max_genome_bp");
        for line in lines {
            let fields = line.split('\t').collect::<Vec<_>>();
            let Some(found) = fields
                .get(set_at)
                .and_then(|name| (0..self.sets.len()).find(|set| self.sets.name(*set) == *name))
            else {
                continue;
            };
            if let Some(bp) = fields.get(bp_at).and_then(|field| field.parse().ok()) {
                expected[found] = bp;
            }
            if let Some(bp) = max_at
                .and_then(|at| fields.get(at))
                .and_then(|field| field.parse().ok())
            {
                bounds[found] = bp;
            }
        }
        self.sets = self.sets.with_scales(expected).with_bounds(bounds);
        self
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    pub fn id(&self, model: &str) -> Option<u16> {
        self.ids.get(model).copied()
    }

    pub fn name(&self, marker: u16) -> &str {
        self.names
            .get(marker as usize)
            .map(String::as_str)
            .unwrap_or_default()
    }
}

fn set_names(table: &str) -> Vec<String> {
    let mut lines = table.lines();
    let header = lines
        .next()
        .unwrap_or_default()
        .split('\t')
        .collect::<Vec<_>>();
    let Some(at) = header.iter().position(|field| *field == "sets") else {
        return Vec::new();
    };
    let mut found = Vec::new();
    for line in lines {
        let Some(listed) = line.split('\t').nth(at) else {
            continue;
        };
        for name in listed.split(',').filter(|name| !name.is_empty()) {
            if !found.iter().any(|held| held == name) {
                found.push(name.to_string());
            }
        }
    }
    found
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
        let directory = tempfile::tempdir()?;
        let hmm = directory.path().join("markers.hmm");
        inflate(HMM_GZ, &hmm)?;
        let proteins = directory.path().join("proteins.faa");

        let (names, called) = {
            let _timer = crate::timing::scope("genes");
            let mut sink = BufWriter::new(std::fs::File::create(&proteins)?);
            let mut called: Vec<orfs::Orf> = Vec::new();
            let names = orfs::call_over(assembly, min_contig_size, |batch| {
                for mut orf in batch {
                    if searchable(&orf.protein) {
                        writeln!(sink, ">{}\n{}", called.len(), orf.protein)?;
                    }
                    if !orf.partial {
                        orf.protein = String::new();
                    }
                    called.push(orf);
                }
                Ok(())
            })?;
            sink.flush()?;
            info!("Called {} genes over {} contigs", called.len(), names.len());
            (names, called)
        };

        HmmerEngine::check_installed()?;
        let engine = HmmerEngine::new(threads, shards);
        let mut hits = {
            let _timer = crate::timing::scope("search");
            engine.search(&hmm, &proteins, directory.path())?
        };
        {
            let _timer = crate::timing::scope("fragments");
            let cut = directory.path().join("fragments.faa");
            let found = write_fragments(&called, &hits, &cut)?;
            if found > 0 {
                let table = engine.search_domains(&hmm, &cut, directory.path())?;
                let bars = fragments::gathering(&hmm)?;
                let rescued = fragments::accepted(&table, &bars, rules.fragment_span);
                debug!("{} markers rescued from {found} cut genes", rescued.len());
                hits.extend(rescued);
            }
        }

        let mut per_contig = vec![Vec::new(); names.len()];
        for (protein, (model, _)) in hits {
            let Some(marker) = set.id(&model) else {
                continue;
            };
            let orf = &called[protein.parse::<usize>()?];
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

fn write_fragments(
    orfs: &[orfs::Orf],
    hits: &HashMap<String, (String, f64)>,
    target: &Path,
) -> Result<usize> {
    let mut sink = BufWriter::new(std::fs::File::create(target)?);
    let mut written = 0;
    for (position, orf) in orfs.iter().enumerate() {
        if !orf.partial || !searchable(&orf.protein) || hits.contains_key(&position.to_string()) {
            continue;
        }
        writeln!(sink, ">{position}\n{}", orf.protein)?;
        written += 1;
    }
    sink.flush()?;
    Ok(written)
}
