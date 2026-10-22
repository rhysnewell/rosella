use std::ffi::OsString;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, PipeReader, Read, Seek, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::thread::JoinHandle;

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use log::info;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::{self as sam, alignment::io::Write as _};

use crate::cli::MappingParams;
use crate::coverage::bam_depth::{Sample, complete};
use crate::coverage::coverage_calculator::{CoverageInputs, Reads};
use crate::external::memory::free_memory;
use crate::external::minimap2::minimap2;

pub struct Mapper<'a> {
    assembly: &'a str,
    directory: PathBuf,
    threads: usize,
    params: &'a MappingParams,
}

impl<'a> Mapper<'a> {
    pub fn new(inputs: &'a CoverageInputs<'a>) -> Self {
        Self {
            assembly: inputs.assembly,
            directory: Path::new(inputs.output_directory).join("bams"),
            threads: inputs.threads,
            params: &inputs.reads.mapping,
        }
    }

    // BAMs are kept, so a rerun into the same directory maps only what it lacks.
    pub fn bams(&self, reads: &[Reads]) -> Result<Vec<Sample>> {
        let (mut short, mut long) = (Vec::new(), Vec::new());
        for set in reads {
            let kept = self.kept(set);
            if kept.is_file() && complete(&kept)? {
                info!("Mapped reads reused: {}", kept.display());
            } else if set.is_long() {
                long.push(set);
            } else {
                short.push(set);
            }
        }
        if self.assembly.is_empty() && !(short.is_empty() && long.is_empty()) {
            bail!("mapping reads needs the assembly passed with -r");
        }
        self.make(&short, &self.params.mapper)?;
        self.make(&long, &self.params.longread_mapper)?;
        reads
            .iter()
            .map(|set| {
                let path = self.kept(set);
                if !path.is_file() {
                    bail!("mapping left no BAM at {}", path.display());
                }
                Ok(Sample {
                    name: set.name().to_string(),
                    path,
                })
            })
            .collect()
    }

    // Named for the sample, so a later run given these BAMs with -b finds the same kept depth.
    fn kept(&self, set: &Reads) -> PathBuf {
        self.directory.join(format!("{}.bam", set.name()))
    }

    // The index is kept beside the BAMs, so a directory filled one sample at a time indexes the
    // assembly once.
    fn index(&self, mapper: &str) -> Result<Index> {
        let preset = preset(mapper);
        let built = self.built(mapper)?;
        let free = free_memory();
        let size = match built.first() {
            Some(index) => std::fs::metadata(&index.path)?.len(),
            None => self.estimated(mapper, preset, free)?,
        };
        // A part takes a quarter of the free memory, since the mapping threads need about as much
        // again and other programs grow while it runs. Each extra part costs one more pass.
        let parts = size.div_ceil((free / 4).max(1)).max(1);
        if parts > 1 {
            info!(
                "The {mapper} index takes {:.1} GB against {:.1} GB free, so it needs {parts} parts",
                size as f64 / 1e9,
                free as f64 / 1e9
            );
        }
        // Free memory moves between samples, so any index split at least as finely is kept
        // rather than built again for a slightly different count.
        if let Some(index) = built
            .into_iter()
            .filter(|index| index.parts >= parts)
            .min_by_key(|index| index.parts)
        {
            return Ok(index);
        }
        let path = self.indexed(mapper, preset, parts)?;
        Ok(Index { path, parts })
    }

    fn path(&self, mapper: &str, parts: u64) -> PathBuf {
        let stem = file_name(self.assembly);
        match parts {
            1 => self.directory.join(format!("{stem}.{mapper}.mmi")),
            _ => self.directory.join(format!("{stem}.{mapper}.{parts}.mmi")),
        }
    }

    fn built(&self, mapper: &str) -> Result<Vec<Index>> {
        let assembly = std::fs::metadata(self.assembly)?.modified()?;
        let prefix = format!("{}.{mapper}.", file_name(self.assembly));
        let mut built = Vec::new();
        for entry in std::fs::read_dir(&self.directory)? {
            let name = entry?.file_name();
            let Some(parts) = name.to_str().and_then(|name| index_parts(name, &prefix)) else {
                continue;
            };
            let path = self.path(mapper, parts);
            let fresh = std::fs::metadata(&path).and_then(|index| index.modified());
            if fresh.is_ok_and(|fresh| fresh >= assembly) {
                built.push(Index { path, parts });
            }
        }
        Ok(built)
    }

    // Building the whole index only to learn its size is what swaps a small machine, so a slice
    // of the assembly small enough to index in a sliver of free memory gives bytes per base.
    fn estimated(&self, mapper: &str, preset: Option<&str>, free: u64) -> Result<u64> {
        let total = bases(self.assembly)?;
        let slice = self
            .directory
            .join(format!("{}.{mapper}.slice.fa", file_name(self.assembly)));
        let taken = sliced(self.assembly, &slice, (free / 64).clamp(1, total.max(1)))?;
        let index = slice.with_extension("mmi");
        self.build(mapper, preset, &slice, &index, None)?;
        let size = std::fs::metadata(&index)?.len();
        std::fs::remove_file(&slice)?;
        std::fs::remove_file(&index)?;
        Ok((size as u128 * total as u128 / taken.max(1) as u128) as u64)
    }

    fn indexed(&self, mapper: &str, preset: Option<&str>, parts: u64) -> Result<PathBuf> {
        let index = self.path(mapper, parts);
        let split = match parts {
            1 => None,
            _ => Some(bases(self.assembly)?.div_ceil(parts)),
        };
        self.build(mapper, preset, Path::new(self.assembly), &index, split)?;
        Ok(index)
    }

    fn build(
        &self,
        mapper: &str,
        preset: Option<&str>,
        reference: &Path,
        index: &Path,
        split: Option<u64>,
    ) -> Result<()> {
        let partial = index.with_extension("mmi.partial");
        let mut args = self.tuned(preset);
        if let Some(split) = split {
            args.extend(["-I".into(), split.to_string().into()]);
        }
        args.extend([
            "-t".into(),
            self.threads.to_string().into(),
            "-d".into(),
            partial.clone().into(),
            reference.into(),
        ]);
        mapping(format!("indexing the assembly for {mapper}"), &args, |_| {
            Ok(())
        })?;
        std::fs::rename(&partial, index)?;
        Ok(())
    }

    fn tuned(&self, preset: Option<&str>) -> Vec<OsString> {
        let mut args = Vec::new();
        if let Some(preset) = preset {
            args.extend(["-x".into(), preset.into()]);
        }
        if let Some(extra) = &self.params.minimap2_params {
            args.extend(extra.split_whitespace().map(OsString::from));
        }
        args
    }

    fn make(&self, sets: &[&Reads], mapper: &str) -> Result<()> {
        if sets.is_empty() {
            return Ok(());
        }
        std::fs::create_dir_all(&self.directory)?;
        let index = self.index(mapper)?;
        sets.iter()
            .try_for_each(|set| self.map(set, mapper, &index))
    }

    // The BAM stays in mapper order with every equally good hit, so the depth pass can split a
    // read that maps as well to several contigs.
    fn map(&self, set: &Reads, mapper: &str, index: &Index) -> Result<()> {
        let files = set.files();
        // minimap2 exits the whole process on a read file it cannot open.
        for file in &files {
            File::open(file).with_context(|| format!("cannot read {file}"))?;
        }
        let mut args: Vec<OsString> = vec!["-a".into(), "--secondary=yes".into()];
        args.extend(self.tuned(preset(mapper)));
        if index.parts > 1 {
            args.extend([
                "--split-prefix".into(),
                self.directory.join(format!("{}.split", set.name())).into(),
            ]);
        }
        args.extend([
            "-t".into(),
            self.threads.to_string().into(),
            index.path.clone().into(),
        ]);
        args.extend(files.into_iter().map(OsString::from));
        let kept = self.kept(set);
        let partial = kept.with_extension("bam.partial");
        let check = whole(set);
        mapping(format!("{} with {mapper}", set.name()), &args, |sam| {
            write_bam(sam, &partial, self.threads)
        })?;
        if let Err(damaged) = check.join().expect("the read file check panicked") {
            std::fs::remove_file(&partial)?;
            return Err(damaged);
        }
        std::fs::rename(&partial, &kept)?;
        Ok(())
    }
}

struct Index {
    path: PathBuf,
    parts: u64,
}

pub fn index_parts(name: &str, prefix: &str) -> Option<u64> {
    match name.strip_prefix(prefix)?.strip_suffix("mmi")? {
        "" => Some(1),
        parts => parts.strip_suffix('.')?.parse().ok(),
    }
}

fn sliced(assembly: &str, slice: &Path, wanted: u64) -> Result<u64> {
    let mut reader = needletail::parse_fastx_file(assembly)?;
    let mut writer = BufWriter::new(File::create(slice)?);
    let mut taken = 0;
    while taken < wanted {
        let Some(record) = reader.next() else {
            break;
        };
        let record = record?;
        writer.write_all(b">")?;
        writer.write_all(record.id())?;
        writer.write_all(b"\n")?;
        writer.write_all(&record.seq())?;
        writer.write_all(b"\n")?;
        taken += record.num_bases() as u64;
    }
    writer.flush()?;
    Ok(taken)
}

fn bases(assembly: &str) -> Result<u64> {
    let mut reader = needletail::parse_fastx_file(assembly)?;
    let mut total = 0;
    while let Some(record) = reader.next() {
        total += record?.num_bases() as u64;
    }
    Ok(total)
}

// minimap2 stops at a damaged gzip block and still exits cleanly, so every read file is read to
// its end beside the mapping and its BAM is kept only if all of them were whole.
fn whole(set: &Reads) -> JoinHandle<Result<()>> {
    let files = set.files();
    std::thread::spawn(move || files.iter().try_for_each(|path| intact(Path::new(path))))
}

pub fn intact(path: &Path) -> Result<()> {
    let mut file = File::open(path)?;
    let mut magic = [0; 2];
    if file.read_exact(&mut magic).is_err() || magic != [0x1f, 0x8b] {
        return Ok(());
    }
    file.rewind()?;
    let mut reads = MultiGzDecoder::new(BufReader::with_capacity(1 << 20, file));
    std::io::copy(&mut reads, &mut std::io::sink()).with_context(|| {
        format!(
            "{} is damaged, so the mapper would stop part way through it",
            path.display()
        )
    })?;
    Ok(())
}

pub const MAPPERS: [(&str, Option<&str>); 6] = [
    ("minimap2-sr", Some("sr")),
    ("minimap2-ont", Some("map-ont")),
    ("minimap2-pb", Some("map-pb")),
    ("minimap2-hifi", Some("map-hifi")),
    ("minimap2-lr-hq", Some("lr:hq")),
    ("minimap2-no-preset", None),
];

fn preset(mapper: &str) -> Option<&'static str> {
    MAPPERS
        .iter()
        .find(|(name, _)| *name == mapper)
        .and_then(|(_, preset)| *preset)
}

fn mapping<T: Send>(
    message: String,
    args: &[OsString],
    read: impl FnOnce(&mut BufReader<PipeReader>) -> Result<T> + Send,
) -> Result<T> {
    let progress = crate::progress::spinning(crate::progress::Stage::MappingReads);
    progress.set_message(message);
    let outcome = minimap2(args, read);
    crate::progress::done(&progress);
    outcome
}

pub fn write_bam(sam: impl BufRead, path: &Path, threads: usize) -> Result<()> {
    let mut reader = sam::io::Reader::new(sam);
    let header = reader.read_header()?;
    let workers = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
    let mut writer = bam::io::Writer::from(bgzf::io::MultithreadedWriter::with_worker_count(
        workers,
        File::create(path)?,
    ));
    writer.write_header(&header)?;
    let mut record = sam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        writer.write_alignment_record(&header, &record)?;
    }
    writer.into_inner().finish()?;
    Ok(())
}

fn file_name(path: &str) -> &str {
    Path::new(path)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or(path)
}
