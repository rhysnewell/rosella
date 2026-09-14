use std::io::IsTerminal;
use std::sync::OnceLock;
use std::time::Duration;

use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use log::{Log, Metadata, Record};

static BARS: OnceLock<MultiProgress> = OnceLock::new();

/// Rosella plumage, crimson through to rose, laid out in the order a run reaches the stages.
#[derive(Debug, Clone, Copy)]
pub enum Stage {
    MappingReads,
    CountingKmers,
    NearestNeighbours,
    Partitioning,
    CallingGenes,
    SearchingModels,
    RefiningBins,
    RescuingUnbinned,
    WritingBins,
}

impl Stage {
    pub const ALL: [Self; 9] = [
        Self::MappingReads,
        Self::CountingKmers,
        Self::NearestNeighbours,
        Self::Partitioning,
        Self::CallingGenes,
        Self::SearchingModels,
        Self::RefiningBins,
        Self::RescuingUnbinned,
        Self::WritingBins,
    ];

    fn name(self) -> &'static str {
        match self {
            Self::MappingReads => "Mapping reads",
            Self::CountingKmers => "Counting k-mers",
            Self::NearestNeighbours => "Nearest neighbours",
            Self::Partitioning => "Partitioning",
            Self::CallingGenes => "Calling genes",
            Self::SearchingModels => "Searching models",
            Self::RefiningBins => "Refining bins",
            Self::RescuingUnbinned => "Rescuing unbinned",
            Self::WritingBins => "Writing bins",
        }
    }

    fn colour(self) -> u8 {
        match self {
            Self::MappingReads => 160,
            Self::CountingKmers => 166,
            Self::NearestNeighbours => 178,
            Self::Partitioning => 107,
            Self::CallingGenes => 72,
            Self::SearchingModels => 68,
            Self::RefiningBins => 62,
            Self::RescuingUnbinned => 97,
            Self::WritingBins => 168,
        }
    }
}

const TRACK: u8 = 238;

fn bars() -> &'static MultiProgress {
    BARS.get_or_init(|| MultiProgress::with_draw_target(ProgressDrawTarget::hidden()))
}

fn styled(template: &str) -> ProgressStyle {
    ProgressStyle::with_template(template)
        .unwrap_or_else(|_| ProgressStyle::default_bar())
        .progress_chars("=> ")
}

/// The bars and the log share stderr, so every record is written while they are cleared away.
pub fn install(logger: env_logger::Logger, quiet: bool) -> Result<(), log::SetLoggerError> {
    let target = match !quiet && std::io::stderr().is_terminal() {
        true => ProgressDrawTarget::stderr(),
        false => ProgressDrawTarget::hidden(),
    };
    let _ = BARS.set(MultiProgress::with_draw_target(target));
    let level = logger.filter();
    log::set_boxed_logger(Box::new(Through(logger)))?;
    log::set_max_level(level);
    Ok(())
}

pub fn counted_template(stage: Stage) -> String {
    let colour = stage.colour();
    format!("{{prefix:<20.{colour}}} [{{bar:28.{colour}/{TRACK}}}] {{pos}}/{{len}} {{msg}}")
}

pub fn spinning_template(stage: Stage) -> String {
    let colour = stage.colour();
    format!("{{prefix:<20.{colour}}} {{spinner:.{colour}}} {{msg}}")
}

pub fn counted(stage: Stage, total: u64) -> ProgressBar {
    let bar = bars().add(ProgressBar::new(total));
    bar.set_style(styled(&counted_template(stage)));
    bar.set_prefix(stage.name());
    bar
}

pub fn spinning(stage: Stage) -> ProgressBar {
    let bar = bars().add(ProgressBar::new_spinner());
    bar.set_style(styled(&spinning_template(stage)));
    bar.set_prefix(stage.name());
    bar.enable_steady_tick(Duration::from_millis(120));
    bar
}

struct Through(env_logger::Logger);

impl Log for Through {
    fn enabled(&self, metadata: &Metadata) -> bool {
        self.0.enabled(metadata)
    }

    fn log(&self, record: &Record) {
        bars().suspend(|| self.0.log(record));
    }

    fn flush(&self) {
        self.0.flush();
    }
}
