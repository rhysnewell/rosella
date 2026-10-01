pub mod clr;
pub mod kmer_counting;
pub mod sketch;

use anyhow::Result;
use needletail::Sequence;
use needletail::parser::SequenceRecord;
use rayon::prelude::*;

// Contigs are read in chunks rather than whole so the parallel pass does not hold the
// assembly in memory beside what it builds.
const CHUNK: usize = 512;

pub(crate) fn each_named(
    assembly: &str,
    wanted: impl Fn(&str) -> bool,
    mut each: impl FnMut(&str, &SequenceRecord) -> Result<()>,
) -> Result<()> {
    let mut reader = needletail::parse_fastx_file(assembly)?;
    while let Some(record) = reader.next() {
        let record = record?;
        let name = crate::contig_id(record.id())?;
        if wanted(name) {
            each(name, &record)?;
        }
    }
    Ok(())
}

// A contig under `floor` still comes back by name, so a caller keeping every name keeps the
// assembly's order without measuring what it will never read.
pub(crate) fn measured<T: Send>(
    assembly: &str,
    floor: usize,
    measure: impl Fn(&[u8]) -> T + Sync,
    mut absorb: impl FnMut(Vec<(String, Option<T>)>),
) -> Result<()> {
    let mut reader = needletail::parse_fastx_file(assembly)?;
    let mut chunk = Vec::with_capacity(CHUNK);
    loop {
        let mut bases = 0;
        while chunk.len() < CHUNK && bases < crate::defaults::CHUNK_BASES {
            let Some(record) = reader.next() else { break };
            let record = record?;
            let sequence =
                (record.num_bases() >= floor).then(|| record.normalize(false).into_owned());
            bases += sequence.as_ref().map_or(0, Vec::len);
            chunk.push((crate::contig_id(record.id())?.to_string(), sequence));
        }
        if chunk.is_empty() {
            return Ok(());
        }
        let measures = chunk
            .par_iter()
            .map(|(_, sequence)| sequence.as_deref().map(&measure))
            .collect::<Vec<_>>();
        absorb(
            chunk
                .drain(..)
                .map(|(name, _)| name)
                .zip(measures)
                .collect(),
        );
    }
}
