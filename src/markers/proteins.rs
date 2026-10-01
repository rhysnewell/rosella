use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use anyhow::{Result, anyhow};

use crate::external::hmmer_engine::HmmerEngine;
use crate::markers::checkm::{self, Counted};
use crate::markers::{ContigMarkers, inflate};
use crate::quality::orfs::Orf;

const UNPLACED: u32 = u32::MAX;
const BINNED_STEM: &str = "binned";
const SHARED_TABLE: &str = "shared.tbl";

// Calling genes again costs more than the search it would save, so a called band keeps its
// proteins on disk until the bins are known and CheckM searches the binned ones alone.
pub(super) struct Proteins {
    directory: tempfile::TempDir,
    pieces: Vec<PathBuf>,
    contig_of: Vec<u32>,
    covered: Vec<u32>,
    engine: HmmerEngine,
}

impl Proteins {
    pub(super) fn keep(
        directory: tempfile::TempDir,
        pieces: Vec<PathBuf>,
        table: &str,
        called: &[Orf],
        covered: impl Iterator<Item = bool>,
        engine: HmmerEngine,
    ) -> Result<Self> {
        std::fs::write(directory.path().join(SHARED_TABLE), table)?;
        Ok(Self {
            directory,
            pieces,
            contig_of: called.iter().map(|orf| orf.contig as u32).collect(),
            covered: covered
                .enumerate()
                .filter_map(|(row, called)| called.then_some(row as u32))
                .collect(),
            engine,
        })
    }

    pub(super) fn renumbered(mut self, placed: &[Option<usize>]) -> Self {
        let to = |row: u32| {
            placed
                .get(row as usize)
                .copied()
                .flatten()
                .map_or(UNPLACED, |contig| contig as u32)
        };
        for contig in &mut self.contig_of {
            *contig = to(*contig);
        }
        self.covered = self
            .covered
            .iter()
            .map(|row| to(*row))
            .filter(|contig| *contig != UNPLACED)
            .collect();
        self
    }

    fn wanted(&self, protein: usize, wanted: &[bool]) -> Option<usize> {
        let contig = *self.contig_of.get(protein)? as usize;
        wanted.get(contig).copied()?.then_some(contig)
    }

    fn search(&self, panel: &checkm::Panel, wanted: &[bool]) -> Result<Option<String>> {
        let directory = self.directory.path();
        let mut sink = self.engine.shards(directory, BINNED_STEM)?;
        for (protein_id, protein) in self.read()? {
            if self.wanted(protein_id, wanted).is_some() {
                sink.write(protein_id, &protein)?;
            }
        }
        let pieces = sink.finish()?;
        if pieces.is_empty() {
            return Ok(None);
        }
        let hmm = directory.join("checkm.hmm");
        inflate(checkm::HMM_GZ, &hmm)?;
        let mut table = std::fs::read_to_string(directory.join(SHARED_TABLE))?;
        table += &self.engine.search(
            &hmm,
            &pieces,
            directory,
            &format!("{:.2}", panel.floor()),
            "checkm",
        )?;
        Ok(Some(table))
    }

    fn read(&self) -> Result<impl Iterator<Item = (usize, String)> + '_> {
        let mut proteins = Vec::new();
        for piece in &self.pieces {
            let mut lines = BufReader::new(File::open(piece)?).lines();
            while let (Some(header), Some(protein)) = (lines.next(), lines.next()) {
                let header = header?;
                let protein_id = header
                    .strip_prefix('>')
                    .and_then(|id| id.parse::<usize>().ok())
                    .ok_or_else(|| anyhow!("{} holds a bad header {header}", piece.display()))?;
                proteins.push((protein_id, protein?));
            }
        }
        Ok(proteins.into_iter())
    }
}

impl ContigMarkers {
    pub(super) fn search_checkm(&mut self, contigs: &[usize]) -> Result<Vec<usize>> {
        let _timer = crate::timing::scope("checkm");
        let mut wanted = vec![false; self.checkm.len()];
        for contig in contigs {
            wanted[*contig] = self.checkm[*contig].is_none();
        }
        let mut searched = Vec::new();
        for kept in std::mem::take(&mut self.proteins) {
            let covered = kept
                .covered
                .iter()
                .map(|contig| *contig as usize)
                .filter(|contig| wanted[*contig])
                .collect::<Vec<_>>();
            if covered.is_empty() {
                continue;
            }
            let panel = &self.set.checkm;
            let mut copies = match kept.search(panel, &wanted)? {
                Some(table) => panel.tally(
                    &table,
                    |protein| kept.wanted(protein, &wanted),
                    wanted.len(),
                ),
                None => vec![Vec::new(); wanted.len()],
            };
            for contig in covered {
                self.checkm[contig] = Some(std::mem::take(&mut copies[contig]));
                wanted[contig] = false;
                searched.push(contig);
            }
        }
        Ok(searched)
    }

    pub(super) fn checkm_copies(&mut self) -> Result<Vec<(Counted, String)>> {
        let wanted = vec![true; self.checkm.len()];
        let mut found = Vec::new();
        for kept in std::mem::take(&mut self.proteins) {
            let Some(table) = kept.search(&self.set.checkm, &wanted)? else {
                continue;
            };
            let counted = self
                .set
                .checkm
                .counted(&table, |protein| kept.wanted(protein, &wanted));
            let ids = counted
                .iter()
                .flat_map(|copy| [Some(copy.protein), copy.partner])
                .flatten()
                .collect::<HashSet<_>>();
            let sequences = kept
                .read()?
                .filter(|(id, _)| ids.contains(id))
                .collect::<HashMap<_, _>>();
            for copy in counted {
                let joined = [Some(copy.protein), copy.partner]
                    .into_iter()
                    .flatten()
                    .filter_map(|id| sequences.get(&id).map(String::as_str))
                    .collect::<String>();
                found.push((copy, joined));
            }
        }
        Ok(found)
    }
}
