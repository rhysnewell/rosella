use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use log::info;

use crate::external::hmmer_engine::HmmerEngine;
use crate::quality::{Quality, orfs};

pub mod fragments;

const HMM_GZ: &[u8] = include_bytes!("../../data/gtdb_markers.hmm.gz");
const TABLE: &str = include_str!("../../data/gtdb_markers.tsv");

pub const DEFAULT_BAR_OFFSET: f64 = 10.0;

/// A marker gene cut by a contig end is still that marker gene, but two halves of one gene on
/// two contigs are not two copies, so presence and duplication read different columns.
#[derive(Clone, Copy, Debug)]
pub struct MarkerRules {
    pub fragment_span: f64,
    pub no_scale_floor: bool,
    pub bar_offset: f64,
}

impl Default for MarkerRules {
    fn default() -> Self {
        Self {
            fragment_span: fragments::DEFAULT_SPAN,
            no_scale_floor: false,
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

#[derive(Clone, Copy, Debug, Default)]
struct Domains {
    bacterial: bool,
    archaeal: bool,
}

pub struct MarkerSet {
    ids: HashMap<String, u16>,
    names: Vec<String>,
    domains: Vec<Domains>,
}

impl MarkerSet {
    pub fn embedded() -> Self {
        Self::parse(TABLE)
    }

    pub fn parse(table: &str) -> Self {
        let mut lines = table.lines();
        let header = lines
            .next()
            .unwrap_or_default()
            .split('\t')
            .collect::<Vec<_>>();
        let column = |name: &str| header.iter().position(|field| *field == name);
        let (Some(name_at), Some(domain_at)) = (column("model_name"), column("domain")) else {
            return Self {
                ids: HashMap::new(),
                names: Vec::new(),
                domains: Vec::new(),
            };
        };
        let mut ids = HashMap::new();
        let mut names = Vec::new();
        let mut domains = Vec::new();
        for line in lines {
            let fields = line.split('\t').collect::<Vec<_>>();
            let (Some(name), Some(domain)) = (fields.get(name_at), fields.get(domain_at)) else {
                continue;
            };
            if ids.contains_key(*name) {
                continue;
            }
            ids.insert((*name).to_string(), domains.len() as u16);
            names.push((*name).to_string());
            domains.push(Domains {
                bacterial: domain.contains("bac120"),
                archaeal: domain.contains("ar53"),
            });
        }
        Self {
            ids,
            names,
            domains,
        }
    }

    pub fn len(&self) -> usize {
        self.domains.len()
    }

    pub fn is_empty(&self) -> bool {
        self.domains.is_empty()
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
        genes: orfs::GeneRules,
        threads: usize,
        shards: Option<usize>,
        rules: MarkerRules,
    ) -> Result<Self> {
        let directory = tempfile::tempdir()?;
        let hmm = directory.path().join("markers.hmm");
        inflate(HMM_GZ, &hmm)?;
        let proteins = directory.path().join("proteins.faa");

        let (names, called) = {
            let _timer = crate::timing::scope("genes");
            let mut sink = BufWriter::new(std::fs::File::create(&proteins)?);
            let mut called: Vec<orfs::Orf> = Vec::new();
            let names = orfs::call_over(assembly, min_contig_size, genes, |batch| {
                for mut orf in batch {
                    if !orf.protein.is_empty() {
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
                info!("{} markers rescued from {found} cut genes", rescued.len());
                hits.extend(rescued);
            }
        }

        let set = MarkerSet::embedded();
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
        let carriers = per_contig.iter().filter(|hits| !hits.is_empty()).count();
        info!(
            "{carriers} of {} contigs carry a single copy marker",
            names.len()
        );
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

    pub fn select(self, names: &[String]) -> Result<ContigMarkers> {
        let index = self
            .names
            .iter()
            .enumerate()
            .map(|(position, name)| (name.as_str(), position))
            .collect::<HashMap<_, _>>();
        let mut per_contig = Vec::with_capacity(names.len());
        for name in names {
            let Some(position) = index.get(name.as_str()) else {
                anyhow::bail!("the marker table does not hold {name}");
            };
            per_contig.push(self.per_contig[*position].clone());
        }
        Ok(ContigMarkers::new(per_contig, self.set, self.rules))
    }
}

pub struct ContigMarkers {
    per_contig: Vec<Vec<Hit>>,
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

    /// A marker bin that reads whole is whole, so the genome floor is only needed while the
    /// scorer under-reports a fragmented genome.
    fn sees_scale(&self) -> bool {
        self.rules.no_scale_floor
    }

    /// Read against whichever domain the bin fills better, since a bin cannot be both.
    fn score(&self, contigs: &[usize]) -> Quality {
        let counts = self.counts(contigs);
        let tally = |in_set: fn(&Domains) -> bool| {
            let (mut present, mut extra, mut total) = (0usize, 0usize, 0usize);
            for (tally, domains) in counts.iter().zip(&self.set.domains) {
                if !in_set(domains) {
                    continue;
                }
                total += 1;
                present += usize::from(tally.any >= 1);
                extra += tally.complete.saturating_sub(1) as usize;
            }
            (present, extra, total)
        };
        let bacterial = tally(|domains| domains.bacterial);
        let archaeal = tally(|domains| domains.archaeal);
        let (present, extra, total) = match archaeal.0 > bacterial.0 {
            true => archaeal,
            false => bacterial,
        };
        if total == 0 {
            return Quality::default();
        }
        Quality {
            completeness: 100.0 * present as f64 / total as f64,
            contamination: 100.0 * extra as f64 / total as f64,
        }
    }
}

impl ContigMarkers {
    pub fn new(mut per_contig: Vec<Vec<Hit>>, set: MarkerSet, rules: MarkerRules) -> Self {
        for hits in &mut per_contig {
            hits.sort_unstable_by_key(|hit| (hit.marker, hit.partial));
        }
        Self {
            per_contig,
            set,
            rules,
        }
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
        if !orf.partial || orf.protein.is_empty() || hits.contains_key(&position.to_string()) {
            continue;
        }
        writeln!(sink, ">{position}\n{}", orf.protein)?;
        written += 1;
    }
    sink.flush()?;
    Ok(written)
}

