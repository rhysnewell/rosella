use clap::{ArgAction, Args};

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Input and output")]
pub struct Common {
    /// Where bins and tables are written. A composition table already here is reused, so give
    /// a run with other inputs its own directory
    #[arg(short, long = "output-directory")]
    pub output_directory: String,

    /// k-mer table made before the run, read in place of counting. Its width sets the k-mer
    /// size and --kmer-size is ignored. One an earlier run left in --output-directory is read
    /// without this
    #[arg(short = 'K', long = "kmer-frequency-file", alias = "kmer-frequencies")]
    pub kmer_frequency_file: Option<String>,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Input and output")]
pub struct Genomes {
    /// Bins to read
    #[arg(short = 'f', long = "genome-fasta-files", num_args = 1.., action = ArgAction::Append)]
    pub genome_fasta_files: Vec<String>,

    /// Directory holding the bins to read
    #[arg(short = 'd', long = "genome-fasta-directory")]
    pub genome_fasta_directory: Option<String>,

    /// Extension of the bins inside --genome-fasta-directory
    #[arg(short = 'x', long = "genome-fasta-extension",
          default_value = crate::defaults::FASTA_EXTENSION)]
    pub genome_fasta_extension: String,
}

impl Genomes {
    pub fn discover(&self) -> anyhow::Result<Vec<std::path::PathBuf>> {
        crate::bins::discover(
            &self.genome_fasta_files,
            self.genome_fasta_directory.as_deref(),
            &self.genome_fasta_extension,
        )
    }
}

// Each stochastic stage draws from its own stream, so a probe can move one and hold the rest.
// The partition and neighbour seeds are fixed, so a user never picks a seed to trust the bins.
#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Reproducibility")]
pub struct SeedParams {
    /// Seed for the samples drawn while binning. The partition and the neighbour search have
    /// fixed seeds
    #[arg(long, default_value_t = 42)]
    pub seed: u64,

    /// Seed for the neighbour graph. A probe: the shipped value is fixed
    #[arg(long = "knn-seed", hide_short_help = true)]
    pub knn: Option<u64>,

    /// Base seed for the partition ensemble. A probe: the shipped value is fixed
    #[arg(id = "partition-seed", long = "partition-seed", hide_short_help = true)]
    pub partition: Option<u64>,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Runtime")]
pub struct Runtime {
    /// Threads for rosella and everything it calls
    #[arg(short, long, default_value_t = crate::defaults::THREADS)]
    pub threads: usize,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Runtime")]
pub struct Logging {
    /// Log at debug level
    #[arg(short, long, action = ArgAction::SetTrue)]
    pub verbose: bool,

    /// Log errors only
    #[arg(short, long, action = ArgAction::SetTrue)]
    pub quiet: bool,
}

// Declared here so the help flags sit with the run controls. Nothing reads the fields.
#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Runtime")]
pub struct HelpFlags {
    /// Print the flags in everyday use
    #[arg(short = 'h', long = "help", action = ArgAction::HelpShort)]
    help: Option<bool>,

    /// Print every flag, including the ones the short help leaves out
    #[arg(short = 'H', long = "full-help", action = ArgAction::HelpLong)]
    full_help: Option<bool>,
}

pub(crate) fn percentage(value: &str) -> Result<f64, String> {
    bounded(value, 0.0, 100.0)
}

pub(crate) fn non_negative(value: &str) -> Result<f64, String> {
    let parsed: f64 = value
        .parse()
        .map_err(|_| format!("`{value}` is not a number"))?;
    if parsed >= 0.0 {
        Ok(parsed)
    } else {
        Err(format!("`{parsed}` is below 0.0"))
    }
}

pub(crate) fn unit_interval(value: &str) -> Result<f64, String> {
    bounded(value, 0.0, 1.0)
}

fn bounded<T>(value: &str, low: T, high: T) -> Result<T, String>
where
    T: std::str::FromStr + PartialOrd + std::fmt::Display + Copy,
{
    let parsed = value
        .parse::<T>()
        .map_err(|_| format!("`{value}` is not a number"))?;
    if parsed < low || parsed > high {
        return Err(format!("`{value}` is outside {low} to {high}"));
    }
    Ok(parsed)
}
