use std::fs::File;
use std::io::BufRead;
use std::num::NonZero;
use std::path::Path;

use anyhow::Result;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::header::record::value::Map;
use noodles_sam::header::record::value::map::Program;
use noodles_sam::header::record::value::map::program::tag::{
    COMMAND_LINE, DESCRIPTION, NAME, VERSION,
};
use noodles_sam::{self as sam, alignment::io::Write as _};

use crate::external::mapper::{preset_secondaries, settings};
use crate::external::strobealign::{self, STROBEALIGN};

// strobealign runs in process with no command line, so its version and seeding decide.
pub fn seeded_alike(header: &sam::Header, settings: &str) -> bool {
    let programs = header.programs().as_ref();
    programs.get(STROBEALIGN.as_bytes()).is_some_and(|program| {
        let field = |tag| {
            program
                .other_fields()
                .get(&tag)
                .map(|value| value.as_slice())
        };
        field(VERSION) == Some(strobealign::VERSION.as_bytes())
            && field(DESCRIPTION) == Some(settings.as_bytes())
    })
}

pub fn mapped_alike(
    header: &sam::Header,
    mapper: &str,
    extra: Option<&str>,
    secondaries: usize,
) -> bool {
    let Some(line) = header
        .programs()
        .as_ref()
        .get(b"minimap2".as_slice())
        .and_then(|program| program.other_fields().get(&COMMAND_LINE))
    else {
        return false;
    };
    let line = String::from_utf8_lossy(line);
    let expected = settings(mapper, extra, secondaries);
    let (found, cap) = capped(options(line.split_whitespace().skip(1)), mapper);
    let (wanted, _) = capped(options(expected.iter().map(String::as_str)), mapper);
    found == wanted && cap == secondaries
}

// minimap2 takes the last -N it is given.
fn capped<'a>(options: Vec<&'a str>, mapper: &str) -> (Vec<&'a str>, usize) {
    let mut cap = preset_secondaries(mapper);
    let mut kept = Vec::new();
    let mut tokens = options.into_iter();
    while let Some(token) = tokens.next() {
        match token {
            "-N" => cap = tokens.next().and_then(|n| n.parse().ok()).unwrap_or(cap),
            option => kept.push(option),
        }
    }
    (kept, cap)
}

// Verbosity, threads and the split prefix leave the hits alone, and the index path ends the
// options.
fn options<'a>(tokens: impl IntoIterator<Item = &'a str>) -> Vec<&'a str> {
    let mut tokens = tokens.into_iter();
    let mut kept = Vec::new();
    while let Some(token) = tokens.next() {
        match token {
            "-v" | "-t" | "--split-prefix" => {
                tokens.next();
            }
            index if index.ends_with(".mmi") => break,
            option => kept.push(option),
        }
    }
    kept
}

pub fn write_bam(sam: impl BufRead, path: &Path, threads: usize) -> Result<()> {
    let mut reader = sam::io::Reader::new(sam);
    let mut header = reader.read_header()?;
    let line = std::env::args().collect::<Vec<_>>().join(" ");
    let rosella = Map::<Program>::builder()
        .insert(NAME, "rosella")
        .insert(VERSION, crate::VERSION)
        .insert(COMMAND_LINE, line)
        .build()?;
    header.programs_mut().add("rosella", rosella)?;
    let workers = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
    let mut writer = bam::io::Writer::from(bgzf::io::MultithreadedWriter::with_worker_count(
        workers,
        File::create(path)?,
    ));
    writer.write_header(&header)?;
    // htslib tools expect each integer tag in its smallest type. Eager parsing stores it so.
    let mut record = sam::alignment::RecordBuf::default();
    while reader.read_record_buf(&header, &mut record)? != 0 {
        writer.write_alignment_record(&header, &record)?;
    }
    writer.into_inner().finish()?;
    Ok(())
}
