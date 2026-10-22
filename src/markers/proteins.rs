use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Result;

use crate::external::hmmer_engine::{self, ProteinFile};
use crate::markers::checkm::{self, Counted};
use crate::markers::{ContigMarkers, inflate};
use crate::quality::orfs::Orf;

const UNPLACED: u32 = u32::MAX;
const BINNED: &str = "binned.faa";
const SHARED_TABLE: &str = "shared.tbl";

// Calling genes again costs more than the search it would save, so a called band keeps its
// proteins on disk until the bins are known and CheckM searches the binned ones alone.
pub(super) struct Proteins {
    directory: tempfile::TempDir,
    proteins: Option<PathBuf>,
    contig_of: Vec<u32>,
    covered: Vec<u32>,
}

impl Proteins {
    pub(super) fn keep(
        directory: tempfile::TempDir,
        proteins: Option<PathBuf>,
        table: &str,
        called: &[Orf],
        covered: impl Iterator<Item = bool>,
    ) -> Result<Self> {
        std::fs::write(directory.path().join(SHARED_TABLE), table)?;
        Ok(Self {
            directory,
            proteins,
            contig_of: called.iter().map(|orf| orf.contig as u32).collect(),
            covered: covered
                .enumerate()
                .filter_map(|(row, called)| called.then_some(row as u32))
                .collect(),
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

    // `keep` sees each protein as it streams into the search, so a caller that needs sequences
    // back holds them without a second read and a caller that does not holds nothing.
    fn search(
        &self,
        panel: &checkm::Panel,
        wanted: &[bool],
        mut keep: impl FnMut(usize, &str),
    ) -> Result<Option<String>> {
        let directory = self.directory.path();
        let mut sink = ProteinFile::create(directory.join(BINNED))?;
        self.each(|protein_id, protein| {
            if self.wanted(protein_id, wanted).is_some() {
                sink.write(protein_id, protein)?;
                keep(protein_id, protein);
            }
            Ok(())
        })?;
        let Some(binned) = sink.finish()? else {
            return Ok(None);
        };
        let hmm = directory.join("checkm.hmm");
        inflate(checkm::HMM_GZ, &hmm)?;
        let mut table = std::fs::read_to_string(directory.join(SHARED_TABLE))?;
        table += &hmmer_engine::search(&hmm, &binned, panel.floor())?;
        Ok(Some(table))
    }

    fn each(&self, visit: impl FnMut(usize, &str) -> Result<()>) -> Result<()> {
        match &self.proteins {
            Some(proteins) => hmmer_engine::each_protein(proteins, visit),
            None => Ok(()),
        }
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
            let mut copies = match kept.search(panel, &wanted, |_, _| {})? {
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
            let mut sequences = HashMap::new();
            let searched = kept.search(&self.set.checkm, &wanted, |id, protein| {
                sequences.insert(id, protein.to_string());
            })?;
            let Some(table) = searched else {
                continue;
            };
            let counted = self
                .set
                .checkm
                .counted(&table, |protein| kept.wanted(protein, &wanted));
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
