use std::fs::File;
use std::num::NonZero;
use std::path::Path;

use anyhow::{Context, Result, bail};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::{self as sam};

// strobealign keeps a secondary that scores at least this share of the best, in per cent.
const SECONDARY_SHARE: i64 = 95;

struct Piece {
    reader: bam::io::Reader<bgzf::io::Reader<File>>,
    next: Option<RecordBuf>,
}

impl Piece {
    fn advance(&mut self, header: &sam::Header) -> Result<()> {
        let mut record = RecordBuf::default();
        self.next = (self.reader.read_record_buf(header, &mut record)? != 0).then_some(record);
        Ok(())
    }

    fn take(&mut self, header: &sam::Header, name: Option<&[u8]>) -> Result<Vec<RecordBuf>> {
        let mut held = Vec::new();
        while self.next.as_ref().is_some_and(|next| named(next) == name) {
            held.extend(self.next.take());
            self.advance(header)?;
        }
        Ok(held)
    }
}

fn named(record: &RecordBuf) -> Option<&[u8]> {
    record.name().map(|name| name.as_ref())
}

// Every part saw every read in the same order, so the parts are read side by side one read at
// a time. MAPQ and mate fields stay as the part that wrote each hit set them.
pub fn merge(paths: &[&Path], out: &Path, threads: usize) -> Result<()> {
    let mut pieces = Vec::with_capacity(paths.len());
    let mut header = None;
    for path in paths {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let mut reader = bam::io::Reader::new(file);
        let found = reader.read_header()?;
        header.get_or_insert(found);
        pieces.push(Piece { reader, next: None });
    }
    let header = header.context("no mapped parts to merge")?;
    for piece in &mut pieces {
        piece.advance(&header)?;
    }
    let workers = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
    let mut writer = bam::io::Writer::from(bgzf::io::MultithreadedWriter::with_worker_count(
        workers,
        File::create(out)?,
    ));
    writer.write_header(&header)?;
    while let Some(first) = pieces[0].next.as_ref() {
        let name = named(first).map(<[u8]>::to_vec);
        let mut group = Vec::new();
        for piece in &mut pieces {
            let held = piece.take(&header, name.as_deref())?;
            if held.is_empty() {
                bail!("the index parts' reads fell out of step");
            }
            group.extend(held);
        }
        for record in resolve(group) {
            writer.write_alignment_record(&header, &record)?;
        }
    }
    if pieces.iter().any(|piece| piece.next.is_some()) {
        bail!("an index part holds reads the others lack");
    }
    writer.into_inner().finish()?;
    Ok(())
}

pub fn resolve(group: Vec<RecordBuf>) -> Vec<RecordBuf> {
    let (last, first): (Vec<_>, Vec<_>) = group
        .into_iter()
        .partition(|record| record.flags().is_last_segment());
    let mut kept = settle(first);
    kept.extend(settle(last));
    kept
}

fn settle(mate: Vec<RecordBuf>) -> Vec<RecordBuf> {
    let (mapped, unmapped): (Vec<_>, Vec<_>) = mate
        .into_iter()
        .partition(|record| !record.flags().is_unmapped());
    let Some(best) = mapped.iter().map(score).max() else {
        return unmapped.into_iter().take(1).collect();
    };
    let primary = mapped
        .iter()
        .position(|record| !record.flags().is_secondary() && score(record) == best)
        .or_else(|| mapped.iter().position(|record| score(record) == best));
    mapped
        .into_iter()
        .enumerate()
        .filter(|(_, record)| score(record) * 100 >= best * SECONDARY_SHARE)
        .map(|(at, mut record)| {
            record
                .flags_mut()
                .set(Flags::SECONDARY, Some(at) != primary);
            record
        })
        .collect()
}

fn score(record: &RecordBuf) -> i64 {
    record
        .data()
        .get(&Tag::ALIGNMENT_SCORE)
        .and_then(|value| value.as_int())
        .unwrap_or(0)
}
