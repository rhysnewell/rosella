use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};

use anyhow::Result;

use crate::external::hmmer_engine::ProteinFile;
use crate::markers::annotation::searchable;
use crate::markers::annotator::write_contigs;
use crate::markers::checkm::Counted;
use crate::quality::orfs;

// A cached annotation holds no proteins, so the chosen contigs are called again.
pub(crate) struct Called {
    pub(crate) proteins: Option<PathBuf>,
    homes: Vec<(usize, bool)>,
    offsets: Vec<u64>,
    spans: HashMap<usize, Range<usize>>,
}

impl Called {
    pub(crate) fn over(
        assembly: &str,
        names: &[String],
        contigs: &BTreeSet<usize>,
        folder: &Path,
    ) -> Result<Self> {
        let subset = folder.join("binned.fna");
        let wanted = contigs
            .iter()
            .map(|contig| names[*contig].clone())
            .collect::<Vec<_>>();
        write_contigs(assembly, &wanted, &subset)?;
        let index = contigs
            .iter()
            .map(|contig| (names[*contig].as_str(), *contig))
            .collect::<HashMap<_, _>>();
        let mut sink = ProteinFile::create(folder.join("proteins.faa"))?;
        let mut called = Self {
            proteins: None,
            homes: Vec::new(),
            offsets: Vec::new(),
            spans: HashMap::new(),
        };
        orfs::call_over(&subset.to_string_lossy(), 0..usize::MAX, |batch, walked| {
            for orf in batch {
                let Some(contig) = index.get(walked[orf.contig].as_str()).copied() else {
                    continue;
                };
                let id = called.homes.len();
                called.offsets.push(sink.offset());
                if searchable(&orf.protein) {
                    sink.write(id, &orf.protein)?;
                }
                called.homes.push((contig, orf.partial()));
                called.spans.entry(contig).or_insert(id..id).end = id + 1;
            }
            Ok(())
        })?;
        called.offsets.push(sink.offset());
        called.proteins = sink.finish()?;
        Ok(called)
    }

    pub(crate) fn home(&self, protein: usize) -> Option<(usize, bool)> {
        self.homes.get(protein).copied()
    }

    pub(crate) fn proteins_in(&self, contigs: &[usize]) -> usize {
        contigs
            .iter()
            .filter_map(|contig| self.spans.get(contig))
            .map(ExactSizeIterator::len)
            .sum()
    }

    fn sequence(&self, file: &mut File, protein: usize) -> Result<String> {
        let (Some(from), Some(to)) = (self.offsets.get(protein), self.offsets.get(protein + 1))
        else {
            return Ok(String::new());
        };
        let mut record = String::new();
        file.seek(SeekFrom::Start(*from))?;
        Read::take(&mut *file, to - from).read_to_string(&mut record)?;
        Ok(record.lines().nth(1).unwrap_or_default().to_string())
    }

    pub(crate) fn joined(&self, file: &mut File, copy: &Counted) -> Result<String> {
        let mut protein = self.sequence(file, copy.protein)?;
        if let Some(partner) = copy.partner {
            protein.push_str(&self.sequence(file, partner)?);
        }
        Ok(protein)
    }

    pub(crate) fn copy_contigs(&self, contigs: &BTreeSet<usize>, target: &Path) -> Result<bool> {
        let Some(proteins) = &self.proteins else {
            return Ok(false);
        };
        let mut source = File::open(proteins)?;
        let mut sink = BufWriter::new(File::create(target)?);
        let mut filled = false;
        for span in contigs.iter().filter_map(|contig| self.spans.get(contig)) {
            let (from, to) = (self.offsets[span.start], self.offsets[span.end]);
            source.seek(SeekFrom::Start(from))?;
            filled |= std::io::copy(&mut Read::take(&mut source, to - from), &mut sink)? > 0;
        }
        sink.flush()?;
        Ok(filled)
    }
}
