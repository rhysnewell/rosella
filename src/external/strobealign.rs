use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};
use std::fmt::Write as _;
use std::io::{BufReader, Write};
use std::path::Path;
use std::sync::Mutex;
use std::sync::mpsc::{channel, sync_channel};
use std::time::Instant;

use anyhow::{Result, anyhow};
use fastrand::Rng;
use log::{debug, info};
use strobealign::aligner::{Aligner, Scores};
use strobealign::chain::get_sorted_chains;
use strobealign::chainer::{Chainer, ChainingParameters};
use strobealign::index::StrobemerIndex;
use strobealign::indexer::make_index;
use strobealign::insertsize::InsertSizeDistribution;
use strobealign::io::SequenceIOError;
use strobealign::io::reads::{
    PeekableSequenceReader, interleaved_record_iterator, open_reads, record_iterator,
};
use strobealign::io::record::RecordPair;
use strobealign::io::xopen::xopen;
use strobealign::modes::align::{
    MappingParameters, SamOutput, align_paired_end_read, align_single_end_read,
};
use strobealign::packed_seq::PackedSeq;
use strobealign::refseq::RefSequence;
use strobealign::seeding::SeedingParameters;

use crate::coverage::coverage_calculator::Reads;
use crate::external::provenance::write_bam;
use crate::external::split_merge::merge;
use crate::report_sink::pending_beside;

pub const STROBEALIGN: &str = "strobealign";
pub const VERSION: &str = env!("ROSELLA_STROBEALIGN_VERSION");

const CHUNK_BASES: usize = 1_000_000;
const FILTER_FRACTION: f64 = 0.0002;
const BANDWIDTH: usize = 1024;

type Pairs = Box<dyn Iterator<Item = Result<RecordPair, SequenceIOError>> + Send>;

struct Part<'a> {
    refseq: RefSequence,
    index: StrobemerIndex,
    seeding: &'a SeedingParameters,
    mapping: MappingParameters,
    chainer: Chainer,
    output: SamOutput,
}

pub fn map(
    assembly: &str,
    set: &Reads,
    threads: usize,
    secondaries: usize,
    bam: &Path,
    budget: u64,
) -> Result<()> {
    let (names, lengths) = contigs(assembly)?;
    let seeding = SeedingParameters::new(read_length(set)?);
    let bytes = index_bytes(&seeding);
    let parts = parts(&lengths, &bytes, budget);
    if parts.len() > 1 {
        info!(
            "Index parts: {} (strobealign index {:.1} GB, {:.1} GB free for it)",
            parts.len(),
            bytes(lengths.iter().sum()) as f64 / 1e9,
            budget as f64 / 1e9
        );
    }
    let header = header(&names, &lengths, &seeding, secondaries);
    if let [part] = parts.as_slice() {
        return map_part(
            assembly,
            part,
            (&seeding, secondaries),
            set,
            threads,
            &header,
            bam,
        );
    }
    let pieces = parts
        .iter()
        .map(|part| {
            let piece = pending_beside(bam)?;
            map_part(
                assembly,
                part,
                (&seeding, secondaries),
                set,
                threads,
                &header,
                piece.path(),
            )?;
            Ok(piece)
        })
        .collect::<Result<Vec<_>>>()?;
    let paths = pieces.iter().map(|piece| piece.path()).collect::<Vec<_>>();
    merge(&paths, bam, threads)
}

// The seed vector is reserved a fifth over its estimate, and each bucket start is a usize.
pub fn index_bytes(seeding: &SeedingParameters) -> impl Fn(usize) -> u64 + use<> {
    let spacing = seeding.syncmer.k - seeding.syncmer.s + 1;
    move |bases| {
        let syncmers = bases / spacing + 1;
        let bits = ((syncmers as f64).log2() as u32).clamp(9, 32) - 1;
        (bases / 4 + syncmers * 16 * 6 / 5 + (8 << bits)) as u64
    }
}

// Each part costs a pass over the reads, so the contigs are dealt longest first to the lightest
// part and no part is left smaller than it has to be.
pub fn parts(lengths: &[usize], bytes: impl Fn(usize) -> u64, budget: u64) -> Vec<Vec<usize>> {
    let total = lengths.iter().sum::<usize>();
    let count = bytes(total)
        .div_ceil(budget.max(1))
        .clamp(1, lengths.len().max(1) as u64);
    let mut lightest = (0..count as usize)
        .map(|part| Reverse((0, part)))
        .collect::<BinaryHeap<_>>();
    let mut order = (0..lengths.len()).collect::<Vec<_>>();
    order.sort_by_key(|at| Reverse(lengths[*at]));
    let mut parts = vec![Vec::new(); count as usize];
    for at in order {
        let Some(Reverse((held, part))) = lightest.pop() else {
            break;
        };
        parts[part].push(at);
        lightest.push(Reverse((held + lengths[at], part)));
    }
    parts.retain(|part| !part.is_empty());
    parts.iter_mut().for_each(|part| part.sort_unstable());
    parts
}

fn contigs(assembly: &str) -> Result<(Vec<String>, Vec<usize>)> {
    let mut reader = needletail::parse_fastx_file(assembly)?;
    let (mut names, mut lengths) = (Vec::new(), Vec::new());
    while let Some(record) = reader.next() {
        let record = record?;
        names.push(crate::contig_id(record.id())?.to_string());
        lengths.push(record.num_bases());
    }
    Ok((names, lengths))
}

fn reference(assembly: &str, part: &[usize]) -> Result<RefSequence> {
    let mut reader = needletail::parse_fastx_file(assembly)?;
    let (mut sequence, mut starts, mut names) = (PackedSeq::new(), Vec::new(), Vec::new());
    let (mut at, mut wanted) = (0, part.iter().peekable());
    while let Some(record) = reader.next() {
        let record = record?;
        if wanted.next_if_eq(&&at).is_some() {
            starts.push(sequence.len());
            names.push(crate::contig_id(record.id())?.to_string());
            sequence.extend(record.seq().iter().copied());
        }
        at += 1;
    }
    Ok(RefSequence::new(sequence, starts, names)?)
}

fn opened(path: &str) -> Result<PeekableSequenceReader> {
    let file = xopen(path).map_err(|error| anyhow!("cannot read {path}: {error}"))?;
    Ok(PeekableSequenceReader::new(open_reads(file)))
}

fn read_length(set: &Reads) -> Result<usize> {
    let head = opened(set.first())?.peek(500)?;
    let bases = head.iter().map(|record| record.len()).sum::<usize>();
    Ok(bases.checked_div(head.len()).unwrap_or(150))
}

fn pairs(set: &Reads) -> Result<Pairs> {
    let first = opened(set.first())?;
    Ok(match set {
        Reads::Interleaved(_) => interleaved_record_iterator(first),
        Reads::Paired(_, second) => record_iterator(first, Some(second))?,
        Reads::Single(_) | Reads::Long(_) => record_iterator(first, None)?,
    })
}

// Every part writes the whole assembly's header, so its hits name contigs the final BAM knows.
fn header(
    names: &[String],
    lengths: &[usize],
    seeding: &SeedingParameters,
    secondaries: usize,
) -> String {
    let mut text = String::from("@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
    for (name, length) in names.iter().zip(lengths) {
        let _ = writeln!(text, "@SQ\tSN:{name}\tLN:{length}");
    }
    let _ = writeln!(
        text,
        "@PG\tID:{STROBEALIGN}\tPN:{STROBEALIGN}\tVN:{VERSION}\tDS:{}",
        described(seeding, secondaries)
    );
    text
}

pub fn settings(set: &Reads, secondaries: usize) -> Result<String> {
    Ok(described(
        &SeedingParameters::new(read_length(set)?),
        secondaries,
    ))
}

fn described(seeding: &SeedingParameters, secondaries: usize) -> String {
    let (syncmer, strobe) = (&seeding.syncmer, &seeding.randstrobe);
    format!(
        "k={} s={} w_min={} w_max={} max_dist={} secondaries={secondaries} \
         filter={FILTER_FRACTION} bandwidth={BANDWIDTH}",
        syncmer.k, syncmer.s, strobe.w_min, strobe.w_max, strobe.max_dist
    )
}

fn map_part(
    assembly: &str,
    part: &[usize],
    (seeding, secondaries): (&SeedingParameters, usize),
    set: &Reads,
    threads: usize,
    header: &str,
    bam: &Path,
) -> Result<()> {
    let started = Instant::now();
    let contigs = part.len();
    let refseq = reference(assembly, part)?;
    let bits = seeding.syncmer.pick_bits(&refseq);
    let (index, _) = make_index(&refseq, seeding.clone(), bits, FILTER_FRACTION, threads);
    let indexed = started.elapsed();
    let part = Part {
        chainer: Chainer::new(index.k(), ChainingParameters::default()),
        mapping: MappingParameters {
            max_secondary: secondaries,
            ..MappingParameters::default()
        },
        output: SamOutput::new(false, false, None, false),
        refseq,
        index,
        seeding,
    };
    let aligner = Aligner::new(Scores::default(), part.index.k(), BANDWIDTH);
    let reads = pairs(set)?;
    let (sam, mut writer) = std::io::pipe()?;
    std::thread::scope(|scope| {
        let written =
            scope.spawn(|| write_bam(BufReader::with_capacity(1 << 20, sam), bam, threads));
        let mapped = writer
            .write_all(header.as_bytes())
            .map_err(anyhow::Error::from)
            .and_then(|()| part.stream(reads, &aligner, threads, &mut writer));
        drop(writer);
        let written = written.join().expect("the BAM writer panicked");
        debug!(
            "strobealign part of {contigs} contigs: index {:.1} s, mapping {:.1} s",
            indexed.as_secs_f64(),
            (started.elapsed() - indexed).as_secs_f64()
        );
        mapped.and(written)
    })
}

impl Part<'_> {
    // Chunks are mapped in any order and written in the order read, so a part's output follows
    // the read file and the parts can be merged read by read.
    fn stream(
        &self,
        reads: Pairs,
        aligner: &Aligner,
        threads: usize,
        out: &mut impl Write,
    ) -> Result<()> {
        let (chunks, taken) = sync_channel(threads.max(1));
        let taken = Mutex::new(Some(taken));
        let (done, finished) = channel();
        std::thread::scope(|scope| {
            scope.spawn(move || {
                for chunk in chunked(reads).enumerate() {
                    if chunks.send(chunk).is_err() {
                        break;
                    }
                }
            });
            for _ in 0..threads.max(1) {
                let (done, taken) = (done.clone(), &taken);
                let mut aligner = aligner.clone();
                scope.spawn(move || {
                    loop {
                        let next = taken.lock().expect("a mapping thread panicked");
                        let Some(Ok((at, chunk))) = next.as_ref().map(|taken| taken.recv()) else {
                            break;
                        };
                        drop(next);
                        if done.send((at, self.chunk(chunk, &mut aligner))).is_err() {
                            break;
                        }
                    }
                });
            }
            drop(done);
            let (mut waiting, mut next) = (BTreeMap::new(), 0);
            let written = finished.iter().try_for_each(|(at, sam)| {
                waiting.insert(at, sam);
                while let Some(sam) = waiting.remove(&next) {
                    out.write_all(&sam?)?;
                    next += 1;
                }
                Ok(())
            });
            // Dropping the queue stops the reader and every mapping thread once its chunk is done.
            if written.is_err() {
                taken.lock().expect("a mapping thread panicked").take();
            }
            written
        })
    }

    fn chunk(
        &self,
        chunk: Vec<Result<RecordPair, SequenceIOError>>,
        aligner: &mut Aligner,
    ) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        let mut rng = Rng::with_seed(0);
        let mut sizes = InsertSizeDistribution::new();
        for record in chunk {
            let (first, second) = record?;
            let (mut chains, chain, hits) = self.chains(&first.sequence, &mut rng);
            let records = match second {
                Some(second) => {
                    let (other, other_chain, other_hits) = self.chains(&second.sequence, &mut rng);
                    let mut both = [chains, other];
                    align_paired_end_read(
                        &first,
                        &second,
                        &self.refseq,
                        &self.mapping,
                        &self.output,
                        self.seeding,
                        &mut sizes,
                        [hits, other_hits],
                        [chain, other_chain],
                        &mut both,
                        aligner,
                        &mut rng,
                    )
                    .0
                }
                None => {
                    align_single_end_read(
                        &first,
                        &self.index,
                        &self.refseq,
                        &self.mapping,
                        &self.output,
                        hits,
                        chain,
                        &mut chains,
                        aligner,
                        &mut rng,
                    )
                    .0
                }
            };
            for record in records {
                writeln!(out, "{record}")?;
            }
        }
        Ok(out)
    }

    fn chains(
        &self,
        sequence: &[u8],
        rng: &mut Rng,
    ) -> (
        Vec<strobealign::chain::Chain>,
        strobealign::details::ChainDetails,
        strobealign::hit::HitsDetails,
    ) {
        get_sorted_chains(
            sequence,
            &self.index,
            &self.chainer,
            self.mapping.rescue_distance,
            self.mapping.mcs_strategy,
            rng,
        )
    }
}

fn chunked(mut reads: Pairs) -> impl Iterator<Item = Vec<Result<RecordPair, SequenceIOError>>> {
    std::iter::from_fn(move || {
        let (mut chunk, mut bases) = (Vec::new(), 0);
        for record in reads.by_ref() {
            bases += record.as_ref().map_or(0, |(first, second)| {
                first.len() + second.as_ref().map_or(0, |second| second.len())
            });
            chunk.push(record);
            if bases > CHUNK_BASES {
                break;
            }
        }
        (!chunk.is_empty()).then_some(chunk)
    })
}
