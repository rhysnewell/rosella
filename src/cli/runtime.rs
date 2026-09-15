use clap::{ArgAction, Args};

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Input and output")]
pub struct Common {
    /// Where bins and the run's tables are written
    #[arg(short, long = "output-directory")]
    pub output_directory: String,

    /// Precomputed tetranucleotide frequency table, in place of counting them
    #[arg(short = 'K', long = "kmer-frequency-file")]
    pub kmer_frequency_file: Option<String>,
}

/// Each stochastic stage draws from its own stream, so a run can hold three still and move
/// the fourth. Unset means the master seed, which is what keeps the default path unchanged.
#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Reproducibility")]
pub struct SeedParams {
    /// Seeds the embedding and every sample taken during clustering
    #[arg(long, default_value_t = 42)]
    pub seed: u64,

    /// Seed for the nearest neighbour graph. Defaults to --seed
    #[arg(long = "knn-seed", hide_short_help = true)]
    pub knn: Option<u64>,

    /// Seed for the samples the codelength score and the refiner take. Defaults to --seed
    #[arg(long = "sample-seed", hide_short_help = true)]
    pub sample: Option<u64>,

    /// Seed for the node order a graph partition visits. Defaults to --seed
    #[arg(id = "partition-seed", long = "partition-seed", hide_short_help = true)]
    pub partition: Option<u64>,
}

#[derive(Args, Debug, Clone)]
#[command(next_help_heading = "Runtime")]
pub struct Runtime {
    /// Threads for the rayon pool and for every tool rosella calls
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

/// Declared here rather than left to clap so the two help flags sit in a section with the
/// rest of the run controls instead of alone above them. Nothing reads the fields; clap
/// acts during parse and needs them only to hang the actions on.
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

pub(crate) fn knn_candidates_in_range(value: &str) -> Result<usize, String> {
    bounded(value, 2, 256)
}

pub(crate) fn above_zero(value: &str) -> Result<f64, String> {
    let parsed: f64 = value
        .parse()
        .map_err(|_| format!("`{value}` is not a number"))?;
    if parsed > 0.0 && parsed.is_finite() {
        Ok(parsed)
    } else {
        Err(format!("`{value}` is not above 0"))
    }
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
