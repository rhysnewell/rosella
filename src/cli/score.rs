use clap::{ArgAction, ArgGroup, Args};

use super::common::Logging;

#[derive(Args, Debug, Clone)]
#[command(group(ArgGroup::new("genomes").required(true).multiple(true)
    .args(["genome_fasta_files", "genome_fasta_directory"])))]
pub struct ScoreArgs {
    /// Assembly the bins were built from
    #[arg(short = 'r', long, alias = "reference")]
    pub assembly: String,

    /// Bins to score
    #[arg(short = 'f', long = "genome-fasta-files", num_args = 1.., action = ArgAction::Append)]
    pub genome_fasta_files: Vec<String>,

    /// Directory holding the bins to score
    #[arg(short = 'd', long = "genome-fasta-directory")]
    pub genome_fasta_directory: Option<String>,

    /// Extension of the bins inside --genome-fasta-directory
    #[arg(short = 'x', long = "genome-fasta-extension", default_value = "fna")]
    pub genome_fasta_extension: String,

    /// Where the quality table is written
    #[arg(short = 'o', long = "output-file", default_value = "quality.tsv")]
    pub output_file: String,

    /// Contigs shorter than this are left out, as they are in a recover run
    #[arg(long = "min-contig-size", default_value = "1500")]
    pub min_contig_size: usize,

    #[arg(short = 't', long, default_value = "10")]
    pub threads: usize,

    #[arg(long = "hmm-shards")]
    pub hmm_shards: Option<u8>,

    #[arg(long = "gene-min-length", default_value = "0")]
    pub gene_min_length: usize,

    #[arg(long = "gene-model-depth", default_value = "0")]
    pub gene_model_depth: usize,

    #[arg(long = "marker-fragment-span", default_value_t = crate::markers::fragments::DEFAULT_SPAN)]
    pub marker_fragment_span: f64,

    #[arg(long = "marker-bar-offset", default_value_t = crate::markers::DEFAULT_BAR_OFFSET,
          value_parser = crate::cli::common::percentage)]
    pub marker_bar_offset: f64,

    #[command(flatten)]
    pub logging: Logging,
}
