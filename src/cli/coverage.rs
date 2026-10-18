use clap::{ArgAction, ArgGroup, Args};

#[derive(Args, Debug, Clone)]
pub struct ReadArgs {
    #[command(flatten)]
    pub coverage: CoverageSource,

    #[command(flatten)]
    pub mapping: MappingParams,
}

/// Where coverage comes from. Any one of these is enough, so clap requires the group
/// rather than any single member.
#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Coverage input")]
#[command(group(ArgGroup::new("coverage-source").required(true).multiple(true).args([
    "read1", "coupled", "interleaved", "single", "longreads", "bam_files",
    "longread_bam_files", "coverage_file", "depth_directory",
])))]
pub struct CoverageSource {
    /// Precomputed coverage table, the lowest tier: whole-contig depth only, no window depth.
    /// Samples from reads or BAMs given beside it are added to it
    #[arg(
        short = 'C',
        long = "coverage-file",
        short_alias = 'i',
        alias = "coverage-values"
    )]
    pub coverage_file: Option<String>,

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

    /// BAM files, one per sample, short or long reads, coordinate sorted or in the mapper's order
    /// as rosella keeps them. In the mapper's order, a read that maps as well to several contigs
    /// is split between them by their depth. Read once for the coverage table and the window depth
    #[arg(short = 'b', long = "bam-files", num_args = 1.., action = ArgAction::Append)]
    pub bam_files: Vec<String>,

    /// Long read BAM files, read the same way as --bam-files
    #[arg(short = 'l', long = "longread-bam-files", num_args = 1.., action = ArgAction::Append)]
    pub longread_bam_files: Vec<String>,

    /// Where each sample's depth read from a BAM is kept: contigs.tsv with each contig's
    /// length, and `<sample>.depth.zst` with its mean, variance and read count per contig and
    /// its 100 bp window means. A later run reuses a sample's file while it is newer than its
    /// BAM. Given alone, every sample kept in it makes the table, in name order. Defaults to
    /// depth/ under --output-directory
    #[arg(long = "depth-directory", alias = "window-depth")]
    pub depth_directory: Option<String>,
}

// CoverM owns the list of mapper names and validates them.
#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Read mapping")]
pub struct MappingParams {
    /// Mapping software for short reads. Any name CoverM accepts
    #[arg(short = 'p', long, default_value = crate::defaults::SHORT_READ_MAPPER)]
    pub mapper: String,

    /// Mapping software for long reads, such as minimap2-hifi or minimap2-lr-hq. Any name
    /// CoverM accepts
    #[arg(long = "longread-mapper", default_value = crate::defaults::LONG_READ_MAPPER)]
    pub longread_mapper: String,

    /// Extra parameters for minimap2, for indexing and for mapping. '-a' and '--secondary=yes'
    /// are always passed
    #[arg(
        long = "minimap2-parameters",
        alias = "minimap2-params",
        allow_hyphen_values = true
    )]
    pub minimap2_params: Option<String>,

    /// Extra parameters for BWA or BWA-MEM2
    #[arg(
        long = "bwa-parameters",
        alias = "bwa-params",
        allow_hyphen_values = true
    )]
    pub bwa_params: Option<String>,
}
