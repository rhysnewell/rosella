pub mod clr;
pub mod halves;
pub mod kmer_counting;
pub mod scan;
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
    measure: impl Fn(&[u8], usize) -> T + Sync,
    mut absorb: impl FnMut(Vec<(String, Option<T>)>),
) -> Result<()> {
    pipelined(
        assembly,
        |record| {
            let bases = record.num_bases();
            let sequence = (bases >= floor).then(|| record.raw_seq().to_vec());
            Ok((crate::contig_id(record.id())?.to_string(), bases, sequence))
        },
        |chunk| {
            let measures = chunk
                .par_iter()
                .map(|(_, bases, raw)| {
                    raw.as_deref()
                        .map(|raw| measure(&raw.normalize(false), *bases))
                })
                .collect::<Vec<_>>();
            absorb(
                chunk
                    .into_iter()
                    .map(|(name, _, _)| name)
                    .zip(measures)
                    .collect(),
            );
            Ok(())
        },
    )
}

// Inflating the assembly is one core's work. Run apart, it reads the next chunk while the
// caller works through this one, rather than leaving the pool idle through every read.
pub(crate) fn pipelined<T: Send>(
    assembly: &str,
    take: impl Fn(&SequenceRecord) -> Result<T> + Send,
    mut each: impl FnMut(Vec<T>) -> Result<()>,
) -> Result<()> {
    if rayon::current_num_threads() == 1 {
        return chunks(assembly, take, |chunk| each(chunk).map(|()| true));
    }
    let (send, receive) = std::sync::mpsc::sync_channel(1);
    std::thread::scope(|scope| {
        let reader =
            scope.spawn(move || chunks(assembly, take, |chunk| Ok(send.send(chunk).is_ok())));
        let consumed = receive.iter().try_for_each(&mut each);
        drop(receive);
        let read = reader
            .join()
            .unwrap_or_else(|panic| std::panic::resume_unwind(panic));
        consumed.and(read)
    })
}

fn chunks<T>(
    assembly: &str,
    take: impl Fn(&SequenceRecord) -> Result<T>,
    mut each: impl FnMut(Vec<T>) -> Result<bool>,
) -> Result<()> {
    let mut reader = needletail::parse_fastx_file(assembly)?;
    loop {
        let mut chunk = Vec::with_capacity(CHUNK);
        let mut bases = 0;
        while chunk.len() < CHUNK && bases < crate::defaults::CHUNK_BASES {
            let Some(record) = reader.next() else { break };
            let record = record?;
            bases += record.num_bases();
            chunk.push(take(&record)?);
        }
        if chunk.is_empty() || !each(chunk)? {
            return Ok(());
        }
    }
}
