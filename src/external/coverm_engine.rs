use std::fs::File;
use std::io::{BufReader, Read, Seek};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::thread::JoinHandle;

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use log::info;

use crate::cli::MappingParams;
use crate::coverage::bam_depth::{Sample, complete};
use crate::coverage::coverage_calculator::{CoverageInputs, Reads};

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
                    bail!("coverm make left no BAM at {}", path.display());
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

    // CoverM names a BAM for its reference and the sample's first read file.
    fn made(&self, reference: &Path, set: &Reads) -> PathBuf {
        self.directory.join(format!(
            "{}.{}.bam",
            file_name(&reference.to_string_lossy()),
            file_name(set.first())
        ))
    }

    // CoverM maps a lone sample against the assembly itself, which indexes it again on every
    // call, so a directory filled one sample at a time keeps the index beside its BAMs.
    fn index(&self, mapper: &str) -> Result<Option<PathBuf>> {
        let Some(preset) = preset(mapper) else {
            return Ok(None);
        };
        let index = self
            .directory
            .join(format!("{}.{mapper}.mmi", file_name(self.assembly)));
        let built = std::fs::metadata(&index).and_then(|index| index.modified());
        let assembly = std::fs::metadata(self.assembly)?.modified()?;
        if built.is_ok_and(|built| built >= assembly) {
            return Ok(Some(index));
        }
        let partial = index.with_extension("mmi.partial");
        let mut command = Command::new("minimap2");
        if let Some(preset) = preset {
            command.arg("-x").arg(preset);
        }
        if let Some(extra) = &self.params.minimap2_params {
            command.args(extra.split_whitespace());
        }
        command
            .arg("-t")
            .arg(self.threads.to_string())
            .arg("-d")
            .arg(&partial)
            .arg(self.assembly);
        run(command, format!("indexing the assembly for {mapper}"))?;
        std::fs::rename(&partial, &index)?;
        Ok(Some(index))
    }

    fn make(&self, sets: &[&Reads], mapper: &str) -> Result<()> {
        if sets.is_empty() {
            return Ok(());
        }
        std::fs::create_dir_all(&self.directory)?;
        let index = self.index(mapper)?;
        if let Some(index) = &index {
            return sets.iter().try_for_each(|set| self.map(set, mapper, index));
        }
        let reference = Path::new(self.assembly);
        let mut command = Command::new("coverm");
        command
            .arg("make")
            .arg("--reference")
            .arg(reference)
            .arg("--output-directory")
            .arg(&self.directory)
            .arg("--threads")
            .arg(self.threads.to_string())
            .arg("--mapper")
            .arg(mapper)
            .arg("--discard-unmapped");
        if let Some(bwa) = &self.params.bwa_params {
            command.arg("--bwa-params").arg(bwa);
        }
        let mut lists: [(&str, Vec<&str>); 4] = [
            ("-1", Vec::new()),
            ("-2", Vec::new()),
            ("--interleaved", Vec::new()),
            ("--single", Vec::new()),
        ];
        for set in sets {
            match set {
                Reads::Paired(forward, reverse) => {
                    lists[0].1.push(forward);
                    lists[1].1.push(reverse);
                }
                Reads::Interleaved(path) => lists[2].1.push(path),
                Reads::Single(path) | Reads::Long(path) => lists[3].1.push(path),
            }
        }
        for (flag, paths) in lists.iter().filter(|(_, paths)| !paths.is_empty()) {
            command.arg(flag).args(paths);
        }

        let check = whole(sets);
        run(command, format!("{} samples with {mapper}", sets.len()))?;
        check.join().expect("the read file check panicked")?;
        for set in sets {
            let made = self.made(reference, set);
            std::fs::rename(&made, self.kept(set))
                .with_context(|| format!("coverm make left no BAM at {}", made.display()))?;
        }
        Ok(())
    }

    // The BAM stays in mapper order with every equally good hit, so the depth pass can split a
    // read that maps as well to several contigs.
    fn map(&self, set: &Reads, mapper: &str, index: &Path) -> Result<()> {
        let mut minimap2 = Command::new("minimap2");
        minimap2.arg("-a").arg("--secondary=yes");
        if let Some(Some(preset)) = preset(mapper) {
            minimap2.arg("-x").arg(preset);
        }
        if let Some(extra) = &self.params.minimap2_params {
            minimap2.args(extra.split_whitespace());
        }
        minimap2.arg("-t").arg(self.threads.to_string()).arg(index);
        match set {
            Reads::Paired(forward, reverse) => minimap2.arg(forward).arg(reverse),
            Reads::Interleaved(path) | Reads::Single(path) | Reads::Long(path) => {
                minimap2.arg(path)
            }
        };
        let kept = self.kept(set);
        let partial = kept.with_extension("bam.partial");
        let mut samtools = Command::new("samtools");
        samtools
            .arg("view")
            .arg("-b")
            .arg("-@")
            .arg(self.threads.to_string())
            .arg("-o")
            .arg(&partial)
            .arg("-");
        let check = whole(&[set]);
        piped(minimap2, samtools, format!("{} with {mapper}", set.name()))?;
        if let Err(damaged) = check.join().expect("the read file check panicked") {
            std::fs::remove_file(&partial)?;
            return Err(damaged);
        }
        std::fs::rename(&partial, &kept)?;
        Ok(())
    }
}

// minimap2 and bwa stop at a damaged gzip block and still exit cleanly, so every read file is
// read to its end beside the mapping and its BAM is kept only if all of them were whole.
fn whole(sets: &[&Reads]) -> JoinHandle<Result<()>> {
    let files: Vec<String> = sets.iter().flat_map(|set| set.files()).collect();
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

fn preset(mapper: &str) -> Option<Option<&'static str>> {
    match mapper {
        "minimap2-sr" => Some(Some("sr")),
        "minimap2-ont" => Some(Some("map-ont")),
        "minimap2-pb" => Some(Some("map-pb")),
        "minimap2-hifi" => Some(Some("map-hifi")),
        "minimap2-lr-hq" => Some(Some("lr:hq")),
        "minimap2-no-preset" => Some(None),
        _ => None,
    }
}

// The mapper's log is drained as it runs, since a full pipe would stall it.
fn piped(mut mapper: Command, mut writer: Command, message: String) -> Result<()> {
    let progress = crate::progress::spinning(crate::progress::Stage::MappingReads);
    progress.set_message(message);
    let mut first = spawn(mapper.stdout(Stdio::piped()).stderr(Stdio::piped()))?;
    let stdout = first.stdout.take().expect("piped stdout");
    let second = spawn(writer.stdin(stdout).stderr(Stdio::piped()));
    let log = first.stderr.take().map(|mut stderr| {
        std::thread::spawn(move || {
            let mut log = String::new();
            stderr.read_to_string(&mut log).map(|_| log)
        })
    });
    let written = second.and_then(|second| Ok(second.wait_with_output()?));
    let mapped = first.wait()?;
    crate::progress::done(&progress);
    let log = log
        .and_then(|log| log.join().ok())
        .transpose()?
        .unwrap_or_default();
    for (program, ok, log) in [
        (&mapper, mapped.success(), log.as_str()),
        (
            &writer,
            written.as_ref().is_ok_and(|output| output.status.success()),
            written.as_ref().map_or("", |output| {
                std::str::from_utf8(&output.stderr).unwrap_or("")
            }),
        ),
    ] {
        if !ok {
            bail!(
                "{} failed: {}",
                program.get_program().to_string_lossy(),
                log.trim()
            );
        }
    }
    written.map(|_| ())
}

fn spawn(command: &mut Command) -> Result<Child> {
    let program = command.get_program().to_string_lossy().into_owned();
    command.spawn().map_err(|error| match error.kind() {
        std::io::ErrorKind::NotFound => {
            anyhow::anyhow!("{program} is not on PATH. Install it or pass sorted BAMs with -b")
        }
        _ => anyhow::Error::new(error).context(format!("running {program}")),
    })
}

fn run(mut command: Command, message: String) -> Result<()> {
    let program = command.get_program().to_string_lossy().into_owned();
    let progress = crate::progress::spinning(crate::progress::Stage::MappingReads);
    progress.set_message(message);
    let outcome = command.output();
    crate::progress::done(&progress);
    let output = match outcome {
        Ok(output) => output,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => bail!(
            "{program} is not on PATH. Install CoverM and minimap2 or pass sorted BAMs with -b"
        ),
        Err(error) => return Err(error).with_context(|| format!("running {program}")),
    };
    if !output.status.success() {
        bail!(
            "{program} failed with {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        );
    }
    Ok(())
}

fn file_name(path: &str) -> &str {
    Path::new(path)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or(path)
}
