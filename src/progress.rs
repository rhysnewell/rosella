use std::io::IsTerminal;
use std::sync::OnceLock;
use std::time::Duration;

use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use log::{Log, Metadata, Record};

static BARS: OnceLock<MultiProgress> = OnceLock::new();

const COUNTED: &str = "{prefix:<20} [{bar:28}] {pos}/{len} {msg}";
const SPINNING: &str = "{prefix:<20} {spinner} {msg}";

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

pub fn counted(stage: &'static str, total: u64) -> ProgressBar {
    let bar = bars().add(ProgressBar::new(total));
    bar.set_style(styled(COUNTED));
    bar.set_prefix(stage);
    bar
}

pub fn spinning(stage: &'static str) -> ProgressBar {
    let bar = bars().add(ProgressBar::new_spinner());
    bar.set_style(styled(SPINNING));
    bar.set_prefix(stage);
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
