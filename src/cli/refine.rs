use clap::*;

use super::common::*;

pub(crate) fn command() -> Command {
    Command::new("refine")
        .about("Refine MAGs using UMAP and HDBSCAN clustering.")
        .arg_required_else_help(true)
        .args(full_help_args())
        .arg(
            Arg::new("assembly")
                .short('r')
                .long("assembly")
                .alias("reference")
                .required_unless_present_any(["full-help", "full-help-roff"])
                .required_unless_present_all(["coverage-file", "kmer-frequency-file"]),
        )
        .arg(output_directory())
        .arg(
            Arg::new("genome-fasta-files")
                .short('f')
                .long("genome-fasta-files")
                .num_args(1..)
                .required_unless_present_any([
                    "full-help",
                    "full-help-roff",
                    "genome-fasta-directory",
                ]),
        )
        .arg(
            Arg::new("genome-fasta-directory")
                .short('d')
                .long("genome-fasta-directory")
                .required_unless_present_any(["full-help", "full-help-roff", "genome-fasta-files"]),
        )
        .arg(
            Arg::new("genome-fasta-extension")
                .short('x')
                .long("genome-fasta-extension")
                .default_value("fna"),
        )
        .arg(
            Arg::new("checkm-results")
                .long("checkm-results")
                .required(false),
        )
        .arg(threads())
        .args(read_inputs())
        .args(mapping_params())
        .args(read_filtering())
        .args(alignment_flags())
        .args(coverage_trimming())
        .arg(min_covered_fraction())
        .arg(coverage_file())
        .arg(seed())
        .arg(kmer_frequency_file())
        .args(binning_params())
        .arg(n_neighbours())
        .args(embedding_overrides())
        .arg(
            Arg::new("max-contamination")
                .long("max-contamination")
                .value_parser(value_parser!(f64))
                .default_value("15.0"),
        )
        .arg(max_retries())
        .arg(
            Arg::new("bin-tag")
                .long("bin-tag")
                .value_parser(value_parser!(String))
                .default_value("refined_1"),
        )
        .args(logging_args())
}
