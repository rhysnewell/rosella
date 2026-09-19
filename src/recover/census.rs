use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use anyhow::Result;
use log::debug;

pub const STAGES_FILE: &str = "stages.tsv";

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

struct Row {
    stage: String,
    bins: usize,
    binned: usize,
    binned_bp: usize,
    unbinned: usize,
    unbinned_bp: usize,
    digest: u64,
    unbinned_digest: u64,
}

/// Fixed keys rather than a hasher seeded per process, so two runs of the same binary can be
/// compared, and folded with a commutative add so bin order and bin numbering do not enter.
fn digest_of(contigs: &mut [usize]) -> u64 {
    contigs.sort_unstable();
    let mut hash = FNV_OFFSET;
    for contig in contigs.iter() {
        for byte in (*contig as u64).to_le_bytes() {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(FNV_PRIME);
        }
    }
    hash
}

/// Every stage logs its own delta in its own shape, so nothing says where the contigs are
/// after each one. This is the same numbers for every stage, and it needs no truth to read.
#[derive(Default)]
pub struct Census {
    rows: Vec<Row>,
}

impl Census {
    pub fn record<B, C, U>(&mut self, stage: &str, bins: B, unbinned: U, lengths: &[usize])
    where
        B: IntoIterator<Item = C>,
        C: IntoIterator<Item = usize>,
        U: IntoIterator<Item = usize>,
    {
        let mut row = Row {
            stage: stage.to_string(),
            bins: 0,
            binned: 0,
            binned_bp: 0,
            unbinned: 0,
            unbinned_bp: 0,
            digest: 0,
            unbinned_digest: 0,
        };
        let mut members = Vec::new();
        for contigs in bins {
            row.bins += 1;
            members.clear();
            members.extend(contigs);
            row.binned += members.len();
            row.binned_bp += members.iter().map(|contig| lengths[*contig]).sum::<usize>();
            row.digest = row.digest.wrapping_add(digest_of(&mut members));
        }
        members.clear();
        members.extend(unbinned);
        row.unbinned = members.len();
        row.unbinned_bp += members.iter().map(|contig| lengths[*contig]).sum::<usize>();
        row.unbinned_digest = digest_of(&mut members);
        self.rows.push(row);
    }

    pub fn write(&self, path: impl AsRef<Path>) -> Result<()> {
        let mut out = BufWriter::new(File::create(path)?);
        writeln!(
            out,
            "stage\tbins\tbinned_contigs\tbinned_bp\tunbinned_contigs\tunbinned_bp\tdigest\tunbinned_digest"
        )?;
        for row in &self.rows {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:016x}\t{:016x}",
                row.stage,
                row.bins,
                row.binned,
                row.binned_bp,
                row.unbinned,
                row.unbinned_bp,
                row.digest,
                row.unbinned_digest
            )?;
            debug!(
                "{} {} bins, {} contigs {} bp binned, {} contigs {} bp unbinned, digest {:016x}",
                row.stage,
                row.bins,
                row.binned,
                row.binned_bp,
                row.unbinned,
                row.unbinned_bp,
                row.digest
            );
        }
        out.flush()?;
        Ok(())
    }
}
