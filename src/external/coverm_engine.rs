use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{Context, Result, bail};
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

    // Sorted BAMs are kept, so a rerun into the same directory maps only what it lacks.
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
        let preset = match mapper {
            "minimap2-sr" => Some("sr"),
            "minimap2-ont" => Some("map-ont"),
            "minimap2-pb" => Some("map-pb"),
            "minimap2-hifi" => Some("map-hifi"),
            "minimap2-lr-hq" => Some("lr:hq"),
            "minimap2-no-preset" => None,
            _ => return Ok(None),
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
        let reference = index.as_deref().unwrap_or(Path::new(self.assembly));
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
        if index.is_some() {
            command.arg("--minimap2-reference-is-index");
        }
        // Secondary alignments never count towards depth, so they are not written at all.
        if mapper.starts_with("minimap2") {
            let extra = self.params.minimap2_params.as_deref().unwrap_or_default();
            command
                .arg("--minimap2-params")
                .arg(format!("--secondary=no {extra}").trim_end());
        }
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

        run(command, format!("{} samples with {mapper}", sets.len()))?;
        for set in sets {
            let made = self.made(reference, set);
            std::fs::rename(&made, self.kept(set))
                .with_context(|| format!("coverm make left no BAM at {}", made.display()))?;
        }
        Ok(())
    }
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
