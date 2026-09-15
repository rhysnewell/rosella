use clap::{ArgAction, ArgGroup, Args};

/// Where coverage comes from. Any one of these is enough, so clap requires the group
/// rather than any single member.
#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Coverage input")]
#[command(group(ArgGroup::new("coverage-source").required(true).multiple(true).args([
    "read1", "coupled", "interleaved", "single", "longreads", "bam_files",
    "longread_bam_files", "coverage_file",
])))]
pub struct CoverageSource {
    /// Precomputed CoverM coverage table, in place of mapping anything. One left in
    /// --output-directory by an earlier run is picked up without this
    #[arg(short = 'C', long = "coverage-file")]
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

    /// Reference sorted BAM files. No read mapping is undertaken
    #[arg(short = 'b', long = "bam-files", num_args = 1.., action = ArgAction::Append)]
    pub bam_files: Vec<String>,

    /// Reference sorted long read BAM files. No read mapping is undertaken
    #[arg(short = 'l', long = "longread-bam-files", num_args = 1.., action = ArgAction::Append)]
    pub longread_bam_files: Vec<String>,
}

/// Passed straight to CoverM, which owns the list of mapper names and validates them.
#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Read mapping")]
pub struct MappingParams {
    /// Mapping software for short reads. Any name CoverM accepts
    #[arg(short = 'p', long)]
    pub mapper: Option<String>,

    /// Mapping software for long reads. Any name CoverM accepts
    #[arg(long = "longread-mapper")]
    pub longread_mapper: Option<String>,

    /// Extra parameters for minimap2, for indexing and for mapping. '-a' is always passed
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

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Read filtering")]
pub struct ReadFiltering {
    /// Exclude reads aligned over fewer bases than this
    #[arg(long = "min-read-aligned-length")]
    pub min_read_aligned_length: Option<u32>,

    /// Exclude reads by percent identity to the reference, between 0 and 100
    #[arg(long = "min-read-percent-identity")]
    pub min_read_percent_identity: Option<f32>,

    /// Exclude reads by aligned percent of their length, between 0 and 100
    #[arg(long = "min-read-aligned-percent", default_value = "0.0")]
    pub min_read_aligned_percent: f32,

    /// As --min-read-aligned-length, but the pair is excluded together
    #[arg(long = "min-read-aligned-length-pair")]
    pub min_read_aligned_length_pair: Option<u32>,

    /// As --min-read-percent-identity, but the pair is excluded together
    #[arg(long = "min-read-percent-identity-pair")]
    pub min_read_percent_identity_pair: Option<f32>,

    /// As --min-read-aligned-percent, but the pair is excluded together
    #[arg(long = "min-read-aligned-percent-pair")]
    pub min_read_aligned_percent_pair: Option<f32>,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Read filtering")]
pub struct AlignmentFlags {
    /// Count only reads mapped in proper pairs
    #[arg(long = "proper-pairs-only", action = ArgAction::SetTrue)]
    pub proper_pairs_only: bool,

    /// Count secondary alignments
    #[arg(long = "include-secondary", action = ArgAction::SetTrue)]
    pub include_secondary: bool,

    /// Skip supplementary alignments
    #[arg(long = "exclude-supplementary", action = ArgAction::SetTrue)]
    pub exclude_supplementary: bool,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Coverage trimming")]
pub struct CoverageTrimming {
    /// Bases to ignore at each contig end
    #[arg(long = "contig-end-exclusion", default_value = "75")]
    pub contig_end_exclusion: usize,

    /// Discard this percent of the lowest coverage positions
    #[arg(long = "trim-min", default_value = "5.0")]
    pub trim_min: f32,

    /// Keep positions up to this percent of the coverage distribution
    #[arg(long = "trim-max", default_value = "95.0")]
    pub trim_max: f32,

    /// Report zero for contigs covered across less than this fraction
    #[arg(long = "min-covered-fraction", default_value = "0.0")]
    pub min_covered_fraction: f32,
}
