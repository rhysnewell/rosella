pub mod assembly_graph;
pub mod bins;
pub mod cli;
pub mod clustering;
pub mod coverage;
pub mod defaults;
pub mod digest;
pub mod embedding;
pub mod external;
pub mod kmers;
pub mod markers;
pub mod palette;
pub mod pool;
pub mod progress;
pub mod quality;
pub mod recover;
pub mod refine;
pub mod report_sink;
pub mod rows;
pub mod seeds;
pub mod tables;
pub mod timing;
pub mod tuning;

#[macro_use]
extern crate anyhow;

use anyhow::Result;
use flate2::read::MultiGzDecoder;
use std::{
    io::{BufRead, BufReader},
    path::Path,
};

pub const AUTHOR_AND_EMAIL: &str = "Rhys J. P. Newell, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology <rhys.newell94 near gmail.com>";

// A FASTA id ends at the first whitespace, so a header carrying assembler annotation still
// matches the bare name a depth table holds.
pub fn contig_id(id: &[u8]) -> Result<&str> {
    Ok(std::str::from_utf8(id)?
        .split_whitespace()
        .next()
        .unwrap_or_default())
}

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

// The magic decides rather than the extension, so a gzipped table under any name reads, and
// a plain one named .gz does not turn into garbage.
pub fn get_file_reader<P: AsRef<Path>>(file_path: P) -> Result<Box<dyn BufRead>> {
    let mut reader = BufReader::new(std::fs::File::open(file_path)?);
    if reader.fill_buf()?.starts_with(&GZIP_MAGIC) {
        return Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))));
    }
    Ok(Box::new(reader))
}
