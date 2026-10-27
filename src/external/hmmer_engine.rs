use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    sync::Once,
};

use anyhow::{Result, anyhow, bail};
use log::debug;
use rayon::prelude::*;
use tongs::{
    alphabet::Alphabet,
    hmm::Hmm,
    hmmfile,
    pipeline::ZSetBy,
    search::{Query, Scan},
    sequence::Sequence,
    subcmd::hmmsearch::write_domtblout,
    tophits::TopHits,
    util::simd_env::FlushToZero,
};

// The proteins are searched a slice at a time, so a large assembly never holds them all at once.
const HELD_RESIDUES: usize = 1 << 28;

pub struct ProteinFile {
    sink: BufWriter<File>,
    path: PathBuf,
    filled: bool,
}

impl ProteinFile {
    pub fn create(path: PathBuf) -> Result<Self> {
        Ok(Self {
            sink: BufWriter::new(File::create(&path)?),
            path,
            filled: false,
        })
    }

    pub fn write(&mut self, id: usize, protein: &str) -> Result<()> {
        writeln!(self.sink, ">{id}\n{protein}")?;
        self.filled = true;
        Ok(())
    }

    pub fn finish(mut self) -> Result<Option<PathBuf>> {
        self.sink.flush()?;
        Ok(self.filled.then_some(self.path))
    }
}

pub fn each_protein(path: &Path, mut visit: impl FnMut(usize, &str) -> Result<()>) -> Result<()> {
    let mut lines = BufReader::new(File::open(path)?).lines();
    while let (Some(header), Some(protein)) = (lines.next(), lines.next()) {
        let header = header?;
        let id = header
            .strip_prefix('>')
            .and_then(|id| id.parse::<usize>().ok())
            .ok_or_else(|| anyhow!("{} holds a bad header {header}", path.display()))?;
        visit(id, &protein?)?;
    }
    Ok(())
}

// The log-sum table fills once. Filling it sets flush to zero on the calling thread. The guard
// gives the thread its mode back.
fn init() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let _mode = FlushToZero::enter();
        tongs::logsum::p7_flogsuminit();
    });
}

fn sequence(abc: &Alphabet, id: usize, protein: &str) -> Result<Sequence> {
    let dsq = abc
        .digitize_checked(protein.as_bytes())
        .map_err(|error| anyhow!("HMMER refused protein {id}: {error}"))?;
    let n = dsq.len() - 2;
    Ok(Sequence {
        name: id.to_string(),
        acc: String::new(),
        desc: String::new(),
        dsq,
        n,
        l: n,
        taxid: -1,
    })
}

fn read(hmm: &Path) -> Result<Vec<Hmm>> {
    let hmms = hmmfile::read_hmm_file_auto(hmm)
        .map_err(|error| anyhow!("HMMER could not read {}: {error}", hmm.display()))?;
    if hmms.is_empty() {
        bail!("{} holds no models", hmm.display());
    }
    Ok(hmms)
}

pub fn align(profile: &str, proteins: &[&str]) -> Result<String> {
    init();
    let hmms = hmmfile::read_hmms_auto(BufReader::new(profile.as_bytes()))
        .map_err(|error| anyhow!("HMMER could not read the profile: {error}"))?;
    let [hmm] = &hmms[..] else {
        bail!("an alignment takes one profile, not {}", hmms.len());
    };
    let abc = Alphabet::new(hmm.abc_type);
    let sequences = proteins
        .iter()
        .enumerate()
        .map(|(at, protein)| sequence(&abc, at, protein))
        .collect::<Result<Vec<_>>>()?;
    Ok(tongs::align::a2m(hmm, &sequences))
}

// Set up as hmmsearch -T floor --domT floor -Z 1 --domZ 1 --noali sets up each query.
fn query(hmm: &Hmm, floor: f64) -> Result<Query<'_>> {
    let mut query = Query::new(hmm);
    let pli = &mut query.pli;
    pli.t = Some(floor);
    pli.by_e = false;
    pli.dom_t = Some(floor);
    pli.dom_by_e = false;
    pli.z = 1.0;
    pli.z_setby = ZSetBy::Option;
    pli.domz = 1.0;
    pli.domz_setby = ZSetBy::Option;
    pli.do_alignment_display = false;
    pli.new_model_thresholds(&hmm.cutoff)
        .map_err(|error| anyhow!("{error} for model {}", hmm.name))?;
    Ok(query)
}

// HMMER writes the table with flush to zero on, so an E-value under the smallest normal double
// prints as 0.
fn rows(query: &Query, scan: Scan) -> Result<String> {
    let _mode = FlushToZero::enter();
    let mut found = TopHits::new();
    found.hits = scan.hits;
    found.sort_by_sortkey();
    found.threshold(&query.pli, 1.0, 1.0);
    let hmm = query.hmm;
    let mut text = Vec::new();
    write_domtblout(&mut text, &hmm.name, hmm.acc.as_deref(), hmm.m, &found, 1.0, 1.0, false);
    Ok(String::from_utf8(text)?)
}

pub fn search(hmm: &Path, proteins: &Path, floor: f64) -> Result<String> {
    search_held(hmm, proteins, floor, HELD_RESIDUES)
}

// One search serves both readings. The reporting floor sits under the lowest score either can
// accept, so the table is a superset of what a gathering-cutoff run would report.
pub fn search_held(hmm: &Path, proteins: &Path, floor: f64, held: usize) -> Result<String> {
    let started = std::time::Instant::now();
    init();
    let hmms = read(hmm)?;
    let abc = Alphabet::new(hmms[0].abc_type);
    let queries = crate::pool::install(|| {
        hmms.par_iter()
            .map(|hmm| query(hmm, floor))
            .collect::<Result<Vec<_>>>()
    })?;
    let progress =
        crate::progress::counted(crate::progress::Stage::SearchingModels, queries.len() as u64);
    let mut table = String::new();
    let mut slices = 0;
    // A slice waits on its slowest model, so the longest start first and the short ones fill in.
    let mut order = (0..queries.len()).collect::<Vec<_>>();
    order.sort_by_key(|model| std::cmp::Reverse(queries[*model].hmm.m));
    let mut targets = Vec::new();
    let mut residues = 0;
    // Z is fixed, so a target's place in the file never reaches a score and every slice starts at 0.
    let mut slice = |targets: &mut Vec<Sequence>| -> Result<()> {
        slices += 1;
        progress.set_length((slices * queries.len()) as u64);
        let found = crate::pool::install(|| {
            order
                .par_iter()
                .with_max_len(1)
                .map(|model| {
                    let query = &queries[*model];
                    let text = rows(query, query.scan(targets, 0));
                    progress.inc(1);
                    text.map(|text| (*model, text))
                })
                .collect::<Result<Vec<_>>>()
        })?;
        let mut tables = vec![String::new(); queries.len()];
        for (model, text) in found {
            tables[model] = text;
        }
        table.extend(tables);
        targets.clear();
        Ok(())
    };
    each_protein(proteins, |id, protein| {
        let held_one = sequence(&abc, id, protein)?;
        residues += held_one.n;
        targets.push(held_one);
        if residues < held {
            return Ok(());
        }
        slice(&mut targets)?;
        residues = 0;
        Ok(())
    })?;
    if !targets.is_empty() {
        slice(&mut targets)?;
    }
    crate::progress::done(&progress);
    debug!(
        "HMMER: {} models, {slices} slices, {:.1}s",
        queries.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(table)
}
