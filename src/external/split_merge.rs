use std::fs::File;
use std::path::Path;

use anyhow::{Context, Result, bail};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
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
    record.name().map(AsRef::as_ref)
}

// Every part saw every read in the same order, so the parts are read side by side one read at
// a time. MAPQ and mate fields stay as the part that wrote each hit set them.
pub fn merge(paths: &[&Path], out: &Path, threads: usize, secondaries: usize) -> Result<()> {
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
    let mut writer = crate::external::provenance::bam_writer(out, threads)?;
    writer.write_header(&header)?;
    while let Some(first) = pieces[0].next.as_ref() {
        let name = named(first).map(<[u8]>::to_vec);
        let mut group = Vec::new();
        for (part, piece) in pieces.iter_mut().enumerate() {
            let held = piece.take(&header, name.as_deref())?;
            if held.is_empty() {
                bail!("the index parts' reads fell out of step");
            }
            group.extend(held.into_iter().map(|record| (part, record)));
        }
        for record in resolve(group, secondaries) {
            writer.write_alignment_record(&header, &record)?;
        }
    }
    if pieces.iter().any(|piece| piece.next.is_some()) {
        bail!("an index part holds reads the others lack");
    }
    writer.into_inner().finish()?;
    Ok(())
}

fn resolve(group: Vec<(usize, RecordBuf)>, secondaries: usize) -> Vec<RecordBuf> {
    let (last, first): (Vec<_>, Vec<_>) = group
        .into_iter()
        .partition(|(_, record)| record.flags().is_last_segment());
    let mut kept = settle(first, secondaries);
    kept.extend(settle(last, secondaries));
    kept
}

// A supplementary record is a piece of its own part's primary, so it is kept only beside that
// primary and never stands for the read.
fn settle(mate: Vec<(usize, RecordBuf)>, secondaries: usize) -> Vec<RecordBuf> {
    let (unmapped, placed): (Vec<_>, Vec<_>) = mate
        .into_iter()
        .partition(|(_, record)| record.flags().is_unmapped());
    let (pieces, mut mapped): (Vec<_>, Vec<_>) = placed
        .into_iter()
        .partition(|(_, record)| record.flags().is_supplementary());
    let Some(best) = mapped.iter().map(|(_, record)| score(record)).max() else {
        return unmapped
            .into_iter()
            .take(1)
            .map(|(_, record)| record)
            .collect();
    };
    let at = mapped
        .iter()
        .position(|(_, record)| !record.flags().is_secondary() && score(record) == best)
        .or_else(|| mapped.iter().position(|(_, record)| score(record) == best))
        .expect("the best score belongs to a mapped record");
    let (home, mut primary) = mapped.remove(at);
    primary.flags_mut().remove(Flags::SECONDARY);
    mapped.retain(|(_, record)| score(record) * 100 >= best * SECONDARY_SHARE);
    mapped.sort_by_key(|(_, record)| std::cmp::Reverse(score(record)));
    std::iter::once(primary)
        .chain(
            pieces
                .into_iter()
                .filter(|(part, _)| *part == home)
                .map(|(_, record)| record),
        )
        .chain(mapped.into_iter().take(secondaries).map(|(_, mut record)| {
            record.flags_mut().insert(Flags::SECONDARY);
            record.data_mut().remove(&Tag::OTHER_ALIGNMENTS);
            record
        }))
        .collect()
}

fn score(record: &RecordBuf) -> i64 {
    record
        .data()
        .get(&Tag::ALIGNMENT_SCORE)
        .and_then(Value::as_int)
        .unwrap_or(0)
}
