use clap::*;

pub(crate) const MAPPING_SOFTWARE_LIST: &[&str] = &[
    "bwa-mem",
    "bwa-mem2",
    "minimap2-sr",
    "minimap2-ont",
    "minimap2-pb",
    "minimap2-hifi",
    "minimap2-no-preset",
];
pub(crate) const DEFAULT_MAPPING_SOFTWARE: &str = "minimap2-sr";

pub(crate) const LONGREAD_MAPPING_SOFTWARE_LIST: &[&str] =
    &["minimap2-ont", "minimap2-pb", "minimap2-hifi"];
pub(crate) const DEFAULT_LONGREAD_MAPPING_SOFTWARE: &str = "minimap2-ont";

const HELP_FLAGS: [&str; 2] = ["full-help", "full-help-roff"];

/// Every read or alignment input is mutually sufficient, so each one is only required when
/// none of the others is present.
const READ_SOURCES: [&str; 8] = [
    "coverage-file",
    "bam-files",
    "read1",
    "coupled",
    "interleaved",
    "single",
    "longreads",
    "longread-bam-files",
];

/// `read1` and `read2` satisfy each other through `requires`, so neither counts as the
/// other's alternative source.
fn read_input(name: &'static str, satisfied_by_others_except: &[&str]) -> Arg {
    let alternatives = READ_SOURCES
        .iter()
        .filter(|source| !satisfied_by_others_except.contains(*source))
        .chain(HELP_FLAGS.iter())
        .copied()
        .collect::<Vec<&'static str>>();
    Arg::new(name)
        .long(name)
        .action(ArgAction::Append)
        .num_args(1..)
        .required_unless_present_any(alternatives)
}

pub(crate) fn full_help_args() -> [Arg; 2] {
    [
        Arg::new("full-help")
            .short('H')
            .long("full-help")
            .required(false)
            .action(ArgAction::SetTrue),
        Arg::new("full-help-roff")
            .long("full-help-roff")
            .required(false)
            .action(ArgAction::SetTrue),
    ]
}

pub(crate) fn output_directory() -> Arg {
    Arg::new("output-directory")
        .short('o')
        .long("output-directory")
        .required_unless_present_any(HELP_FLAGS)
}

pub(crate) fn threads() -> Arg {
    Arg::new("threads")
        .short('t')
        .long("threads")
        .value_parser(value_parser!(usize))
        .default_value("10")
}

pub(crate) fn read_inputs() -> [Arg; 8] {
    [
        read_input("read1", &["read1"]).short('1').requires("read2"),
        read_input("read2", &["read1"]).short('2').requires("read1"),
        read_input("coupled", &["coupled"]).short('c'),
        read_input("interleaved", &["interleaved"]),
        read_input("single", &["single"]),
        read_input("longreads", &["longreads"]),
        read_input("bam-files", &["bam-files"]).short('b'),
        read_input("longread-bam-files", &["longread-bam-files"]).short('l'),
    ]
}

pub(crate) fn mapping_params() -> [Arg; 5] {
    [
        Arg::new("mapper")
            .short('p')
            .long("mapper")
            .value_parser(MAPPING_SOFTWARE_LIST.iter().collect::<Vec<_>>())
            .default_value(DEFAULT_MAPPING_SOFTWARE),
        Arg::new("longread-mapper")
            .long("longread-mapper")
            .value_parser(LONGREAD_MAPPING_SOFTWARE_LIST.iter().collect::<Vec<_>>())
            .default_value(DEFAULT_LONGREAD_MAPPING_SOFTWARE),
        Arg::new("minimap2-params")
            .long("minimap2-parameters")
            .alias("minimap2-params")
            .allow_hyphen_values(true),
        Arg::new("minimap2-reference-is-index").long("minimap2-reference-is-index"),
        Arg::new("bwa-params")
            .long("bwa-parameters")
            .alias("bwa-params")
            .allow_hyphen_values(true),
    ]
}

pub(crate) fn read_filtering() -> [Arg; 6] {
    [
        Arg::new("min-read-aligned-length")
            .long("min-read-aligned-length")
            .value_parser(value_parser!(u32)),
        Arg::new("min-read-percent-identity")
            .long("min-read-percent-identity")
            .value_parser(value_parser!(f32)),
        Arg::new("min-read-aligned-percent")
            .long("min-read-aligned-percent")
            .value_parser(value_parser!(f32))
            .default_value("0.0"),
        Arg::new("min-read-aligned-length-pair")
            .long("min-read-aligned-length-pair")
            .value_parser(value_parser!(u32)),
        Arg::new("min-read-percent-identity-pair")
            .long("min-read-percent-identity-pair")
            .value_parser(value_parser!(f32)),
        Arg::new("min-read-aligned-percent-pair")
            .long("min-read-aligned-percent-pair")
            .value_parser(value_parser!(f32)),
    ]
}

pub(crate) fn min_covered_fraction() -> Arg {
    Arg::new("min-covered-fraction")
        .long("min-covered-fraction")
        .value_parser(value_parser!(f32))
        .default_value("0.0")
}

pub(crate) fn alignment_flags() -> [Arg; 3] {
    [
        Arg::new("proper-pairs-only")
            .long("proper-pairs-only")
            .action(ArgAction::SetTrue),
        Arg::new("include-secondary")
            .long("include-secondary")
            .action(ArgAction::SetTrue),
        Arg::new("exclude-supplementary")
            .long("exclude-supplementary")
            .action(ArgAction::SetTrue),
    ]
}

pub(crate) fn coverage_trimming() -> [Arg; 3] {
    [
        Arg::new("contig-end-exclusion")
            .long("contig-end-exclusion")
            .value_parser(value_parser!(usize))
            .default_value("75"),
        Arg::new("trim-min")
            .long("trim-min")
            .value_parser(value_parser!(f32))
            .default_value("5.0"),
        Arg::new("trim-max")
            .long("trim-max")
            .value_parser(value_parser!(f32))
            .default_value("95.0"),
    ]
}

pub(crate) fn coverage_file() -> Arg {
    Arg::new("coverage-file")
        .long("coverage-file")
        .short('C')
        .value_parser(value_parser!(String))
        .required_unless_present_any([
            "read1",
            "read2",
            "coupled",
            "interleaved",
            "single",
            "longreads",
            "longread-bam-files",
            "bam-files",
            "full-help",
            "full-help-roff",
        ])
}

pub(crate) fn seed() -> Arg {
    Arg::new("seed")
        .long("seed")
        .value_parser(value_parser!(u64))
        .default_value("42")
}

pub(crate) fn kmer_frequency_file() -> Arg {
    Arg::new("kmer-frequency-file")
        .long("kmer-frequency-file")
        .short('K')
        .value_parser(value_parser!(String))
}

pub(crate) fn binning_params() -> [Arg; 4] {
    [
        Arg::new("min-contig-size")
            .long("min-contig-size")
            .value_parser(value_parser!(usize))
            .default_value("1500"),
        Arg::new("min-bin-size")
            .long("min-bin-size")
            .value_parser(value_parser!(usize))
            .default_value("200000"),
        Arg::new("max-bin-size")
            .long("max-bin-size")
            .value_parser(value_parser!(usize))
            .default_value("15000000"),
        Arg::new("min-contig-count")
            .long("min-contig-count")
            .value_parser(value_parser!(usize))
            .default_value("10"),
    ]
}

pub(crate) fn n_neighbours() -> Arg {
    Arg::new("n-neighbours")
        .long("n-neighbours")
        .alias("n-neighbors")
        .value_parser(value_parser!(usize))
        .default_value("100")
}

pub(crate) fn embedding_overrides() -> [Arg; 3] {
    [
        Arg::new("n-components")
            .long("n-components")
            .value_parser(n_components_in_range),
        Arg::new("umap-a")
            .long("umap-a")
            .value_parser(umap_a_in_range),
        Arg::new("umap-b")
            .long("umap-b")
            .value_parser(umap_b_in_range),
    ]
}

/// The embedding overrides bypass the bounds the derived values are clamped to, so without
/// a range a typo reaches the optimiser and comes back as a silently useless embedding.
fn n_components_in_range(value: &str) -> Result<usize, String> {
    bounded(value, 2, 100)
}

fn umap_a_in_range(value: &str) -> Result<f32, String> {
    bounded(value, 0.01, 10.0)
}

fn umap_b_in_range(value: &str) -> Result<f32, String> {
    bounded(value, 0.01, 5.0)
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

pub(crate) fn max_retries() -> Arg {
    Arg::new("max-retries")
        .long("max-retries")
        .value_parser(value_parser!(usize))
        .default_value("5")
}

pub(crate) fn logging_args() -> [Arg; 2] {
    [
        Arg::new("verbose")
            .short('v')
            .long("verbose")
            .action(ArgAction::SetTrue),
        Arg::new("quiet")
            .short('q')
            .long("quiet")
            .action(ArgAction::SetTrue),
    ]
}
