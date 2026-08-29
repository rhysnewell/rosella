pub mod cli;
pub mod clustering;
pub mod coverage;
pub mod embedding;
pub mod external;
pub mod kmers;
pub mod recover;
pub mod refine;
pub mod seeds;
pub mod timing;

#[macro_use]
extern crate anyhow;

use anyhow::Result;
use std::{
    io::{BufRead, BufReader},
    path::Path,
};

pub const AUTHOR: &str = "Rhys J. P. Newell, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology";
pub const AUTHOR_AND_EMAIL: &str = "Rhys J. P. Newell, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology <rhys.newell94 near gmail.com>";
pub const EMAIL: &str = "rhys.newell94 near gmail.com";

// Enum for exclusion out here so long read can find it
pub enum GenomeExclusionTypes {
    SeparatorType,
    NoneType,
    GenomesAndContigsType,
}

/// read any file into a buffered reader, optionally unzipping it
pub fn get_file_reader<P: AsRef<Path>>(file_path: P) -> Result<Box<dyn BufRead>> {
    let reader = BufReader::new(Box::new(std::fs::File::open(file_path)?));
    Ok(Box::new(reader))
}
