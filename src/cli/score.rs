use clap::{ArgAction, ArgGroup, Args};

use super::common::{Logging, MarkerParams};

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
    #[arg(short = 'x', long = "genome-fasta-extension", default_value = crate::defaults::FASTA_EXTENSION)]
    pub genome_fasta_extension: String,

    /// Where the quality table is written
    #[arg(short = 'o', long = "output-file", default_value = crate::defaults::QUALITY_FILE)]
    pub output_file: String,

    /// Contigs shorter than this are left out, as they are in a recover run
    #[arg(long = "min-contig-size", default_value_t = crate::defaults::MIN_CONTIG_SIZE)]
    pub min_contig_size: usize,

    /// Threads for the gene search
    #[arg(short = 't', long, default_value_t = crate::defaults::THREADS)]
    pub threads: usize,

    #[command(flatten)]
    pub markers: MarkerParams,

    #[command(flatten)]
    pub logging: Logging,
}
