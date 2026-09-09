use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use log::info;

use crate::external::hmmer_engine::HmmerEngine;
use crate::quality::{Quality, orfs};

const HMM_GZ: &[u8] = include_bytes!("../../data/gtdb_markers.hmm.gz");
const TABLE: &str = include_str!("../../data/gtdb_markers.tsv");

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
                domains: Vec::new(),
            };
        };
        let mut ids = HashMap::new();
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
            domains.push(Domains {
                bacterial: domain.contains("bac120"),
                archaeal: domain.contains("ar53"),
            });
        }
        Self { ids, domains }
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
}

pub struct MarkerAnnotation {
    names: Vec<String>,
    per_contig: Vec<Vec<Hit>>,
    set: MarkerSet,
}

impl MarkerAnnotation {
    pub fn build(assembly: &str, min_contig_size: usize, threads: usize) -> Result<Self> {
        let (names, contigs) = orfs::read_over(assembly, min_contig_size)?;
        let called = orfs::call(&contigs, threads)?;
        info!(
            "Called {} genes over {} contigs",
            called.len(),
            contigs.len()
        );

        let directory = tempfile::tempdir()?;
        let hmm = directory.path().join("markers.hmm");
        inflate(HMM_GZ, &hmm)?;
        let proteins = directory.path().join("proteins.faa");
        write_proteins(&called, &proteins)?;
        let hits = HmmerEngine::new(threads).search(&hmm, &proteins, directory.path())?;

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
        })
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
        Ok(ContigMarkers::new(per_contig, self.set))
    }
}

pub struct ContigMarkers {
    per_contig: Vec<Vec<Hit>>,
    set: MarkerSet,
}

impl crate::quality::Scorer for ContigMarkers {
    /// A marker set is a fixed denominator, so presence is completeness and every second copy
    /// is contamination, read against whichever domain the bin fills better.
    fn score(&self, contigs: &[usize]) -> Quality {
        let counts = self.counts(contigs);
        let tally = |in_set: fn(&Domains) -> bool| {
            let (mut present, mut extra, mut total) = (0usize, 0usize, 0usize);
            for (count, domains) in counts.iter().zip(&self.set.domains) {
                if !in_set(domains) {
                    continue;
                }
                total += 1;
                present += usize::from(*count >= 1);
                extra += count.saturating_sub(1) as usize;
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
    pub fn new(mut per_contig: Vec<Vec<Hit>>, set: MarkerSet) -> Self {
        for hits in &mut per_contig {
            hits.sort_unstable_by_key(|hit| (hit.marker, hit.partial));
        }
        Self { per_contig, set }
    }

    fn counts(&self, contigs: &[usize]) -> Vec<u32> {
        let mut counts = vec![0u32; self.set.len()];
        for contig in contigs {
            for hit in &self.per_contig[*contig] {
                if !hit.partial {
                    counts[hit.marker as usize] += 1;
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

fn write_proteins(orfs: &[orfs::Orf], target: &Path) -> Result<()> {
    let mut sink = BufWriter::new(std::fs::File::create(target)?);
    for (position, orf) in orfs.iter().enumerate() {
        if orf.protein.is_empty() {
            continue;
        }
        writeln!(sink, ">{position}\n{}", orf.protein)?;
    }
    sink.flush()?;
    Ok(())
}
