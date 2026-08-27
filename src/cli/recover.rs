use clap::*;

use super::common::*;

pub(crate) fn command() -> Command {
    Command::new("recover")
        .about("Recover MAGs from contigs using UMAP and HDBSCAN clustering.")
        .arg_required_else_help(true)
        .args(full_help_args())
        .arg(
            Arg::new("assembly")
                .short('r')
                .long("assembly")
                .alias("reference")
                .required_unless_present_any(["full-help", "full-help-roff"]),
        )
        .arg(output_directory())
        .args(read_inputs())
        .arg(threads())
        .args(mapping_params())
        .args(read_filtering())
        .arg(min_covered_fraction())
        .args(alignment_flags())
        .args(coverage_trimming())
        .arg(coverage_file())
        .arg(seed())
        .arg(kmer_frequency_file())
        .args(binning_params())
        .arg(n_neighbours())
        .args(embedding_overrides())
        .arg(max_retries())
        .arg(Arg::new("refine").long("refine").action(ArgAction::SetTrue))
        .args(logging_args())
}
