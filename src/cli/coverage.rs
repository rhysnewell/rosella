use clap::builder::PossibleValuesParser;
use clap::{ArgAction, ArgGroup, Args};

use crate::external::mapper::MAPPERS;

#[derive(Args, Debug, Clone)]
pub struct ReadArgs {
    #[command(flatten)]
    pub coverage: CoverageSource,

    #[command(flatten)]
    pub mapping: MappingParams,
}

// Any one source is enough, so clap requires the group rather than any single member.
#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Coverage input")]
#[command(group(ArgGroup::new("coverage-source").required(true).multiple(true).args([
    "read1", "coupled", "interleaved", "single", "longreads", "bam_files",
    "longread_bam_files", "coverage_file", "depth_directory",
])))]
pub struct CoverageSource {
    /// Coverage table made before the run. It gives each contig one depth and no window
    /// depth. Samples from reads or BAMs given beside it are added to it
    #[arg(
        short = 'C',
        long = "coverage-file",
        short_alias = 'i',
        alias = "coverage-values"
    )]
    pub coverage_file: Option<String>,

    #[command(flatten)]
    pub samples: Samples,

    /// Where each sample's depth read from a BAM is kept: contigs.tsv with each contig's
    /// length, and `<sample>.depth.zst` with its mean, variance and read count per contig and
    /// its 100 bp window means. A later run reuses a sample's file while it is newer than its
    /// BAM. Given alone, every sample kept in it makes the table, in name order. Defaults to
    /// depth/ under --output-directory
    #[arg(long = "depth-directory", alias = "window-depth")]
    pub depth_directory: Option<String>,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Coverage input")]
pub struct Samples {
    /// Forward read files, paired with --read2
    #[arg(short = '1', long, num_args = 1.., action = ArgAction::Append, requires = "read2")]
    pub read1: Vec<String>,

    /// Reverse read files, paired with --read1
    #[arg(short = '2', long, num_args = 1.., action = ArgAction::Append, requires = "read1")]
    pub read2: Vec<String>,

    /// Paired read files, given as forward and reverse in turn
    #[arg(short = 'c', long, num_args = 1.., action = ArgAction::Append)]
    pub coupled: Vec<String>,

    /// Interleaved paired read files
    #[arg(long, num_args = 1.., action = ArgAction::Append)]
    pub interleaved: Vec<String>,

    /// Unpaired read files
    #[arg(long, num_args = 1.., action = ArgAction::Append)]
    pub single: Vec<String>,

    /// Long read files
    #[arg(long, num_args = 1.., action = ArgAction::Append)]
    pub longreads: Vec<String>,

    /// BAM files, one per sample, short or long reads. Coordinate sorted, or in the mapper's
    /// order as rosella keeps them. In the mapper's order a read that maps equally well to
    /// several contigs is shared between them by their depth
    #[arg(short = 'b', long = "bam-files", num_args = 1.., action = ArgAction::Append)]
    pub bam_files: Vec<String>,

    /// Long read BAM files, read the same way as --bam-files
    #[arg(short = 'l', long = "longread-bam-files", num_args = 1.., action = ArgAction::Append)]
    pub longread_bam_files: Vec<String>,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Read mapping")]
pub struct MappingParams {
    /// minimap2 preset for short reads
    #[arg(short = 'p', long, default_value = crate::defaults::SHORT_READ_MAPPER, value_parser = mappers())]
    pub mapper: String,

    /// minimap2 preset for long reads
    #[arg(long = "longread-mapper", default_value = crate::defaults::LONG_READ_MAPPER, value_parser = mappers())]
    pub longread_mapper: String,

    /// Extra minimap2 options for indexing and mapping. '-a' and '--secondary=yes' are always
    /// passed
    #[arg(
        long = "minimap2-parameters",
        alias = "minimap2-params",
        allow_hyphen_values = true
    )]
    pub minimap2_params: Option<String>,
}

fn mappers() -> PossibleValuesParser {
    PossibleValuesParser::new(MAPPERS.map(|(name, _)| name))
}
