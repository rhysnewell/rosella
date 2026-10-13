use std::borrow::Cow;
use std::io::{IsTerminal, Write};
use std::sync::{Mutex, MutexGuard, OnceLock, PoisonError};
use std::time::Duration;

use env_logger::Builder;
use env_logger::fmt::style::{Ansi256Color, Color, Style};
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressState, ProgressStyle};
use log::{Level, LevelFilter, Log, Metadata, Record};

use crate::palette;

static BARS: OnceLock<MultiProgress> = OnceLock::new();
// A pass that searches neighbours and partitions more than once holds one bar, so the run prints
// a line per pass rather than one per inner call.
static HOLD: Mutex<Option<Hold>> = Mutex::new(None);

// One width for the stage names and the log levels alike, so a bar and an info line put their
// bodies in the same column.
const GUTTER: usize = 20;
const BAR: usize = 22;
const TICK: Duration = Duration::from_millis(80);
// Bars inside a stage, one per bin, finish in well under this and would flood the screen.
const KEPT: Duration = Duration::from_secs(2);
const FILL: &str = "█▉▊▋▌▍▎▏ ";
const TICKS: &str = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏ ";

// Rosella plumage across the genus, rose through to lilac, in the order a run reaches the
// stages. Pale headed yellows and blues sit in the middle where most of a run is spent.
#[derive(Debug, Clone, Copy)]
pub enum Stage {
    MappingReads,
    ReadingDepth,
    CountingKmers,
    LearningSpace,
    NearestNeighbours,
    Partitioning,
    CallingGenes,
    SearchingModels,
    RefiningBins,
    RescuingUnbinned,
    AttachingContigs,
    WritingBins,
}

impl Stage {
    pub const ALL: [Self; 12] = [
        Self::MappingReads,
        Self::ReadingDepth,
        Self::CountingKmers,
        Self::LearningSpace,
        Self::NearestNeighbours,
        Self::Partitioning,
        Self::CallingGenes,
        Self::SearchingModels,
        Self::RefiningBins,
        Self::RescuingUnbinned,
        Self::AttachingContigs,
        Self::WritingBins,
    ];

    fn name(self) -> &'static str {
        match self {
            Self::MappingReads => "Mapping reads",
            Self::ReadingDepth => "Reading depth",
            Self::CountingKmers => "Counting k-mers",
            Self::LearningSpace => "Learning space",
            Self::NearestNeighbours => "Nearest neighbours",
            Self::Partitioning => "Partitioning",
            Self::CallingGenes => "Calling genes",
            Self::SearchingModels => "Searching models",
            Self::RefiningBins => "Refining bins",
            Self::RescuingUnbinned => "Rescuing unbinned",
            Self::AttachingContigs => "Attaching contigs",
            Self::WritingBins => "Writing bins",
        }
    }

    fn colour(self) -> u8 {
        match self {
            Self::MappingReads => palette::ROSE,
            Self::ReadingDepth => palette::BLUSH,
            Self::CountingKmers => palette::APRICOT,
            Self::LearningSpace => palette::PEACH,
            Self::NearestNeighbours => palette::GOLD,
            Self::Partitioning => palette::LIME,
            Self::CallingGenes => palette::MINT,
            Self::SearchingModels => palette::TEAL,
            Self::RefiningBins => palette::SKY,
            Self::RescuingUnbinned => palette::LILAC,
            Self::AttachingContigs => palette::LAVENDER,
            Self::WritingBins => palette::ORCHID,
        }
    }
}

fn bars() -> &'static MultiProgress {
    BARS.get_or_init(|| MultiProgress::with_draw_target(ProgressDrawTarget::hidden()))
}

fn styled(template: &str) -> ProgressStyle {
    ProgressStyle::with_template(template)
        .unwrap_or_else(|_| ProgressStyle::default_bar())
        .with_key(
            "stamp",
            |_: &ProgressState, out: &mut dyn std::fmt::Write| {
                let _ = out.write_str(&stamp());
            },
        )
        .progress_chars(FILL)
        .tick_chars(TICKS)
}

fn stamp() -> String {
    jiff::Timestamp::now()
        .round(jiff::Unit::Second)
        .map(|now| now.to_string())
        .unwrap_or_default()
}

fn tinted(colour: u8) -> Style {
    Style::new().fg_color(Some(Color::Ansi256(Ansi256Color(colour))))
}

fn level_tint(level: Level) -> Style {
    match level {
        Level::Error => tinted(palette::CORAL).bold(),
        Level::Warn => tinted(palette::APRICOT),
        Level::Info => tinted(palette::GREY),
        _ => tinted(palette::TRACK),
    }
}

// Debug lines say which module spoke, since at that level the level itself says nothing.
fn gutter<'a>(record: &Record<'a>) -> &'a str {
    match record.level() {
        Level::Error => "error",
        Level::Warn => "warn",
        Level::Info => "info",
        _ => record.target().rsplit("::").next().unwrap_or_default(),
    }
}

// The bars and the log share stderr, so every record is written while they are cleared away.
// Off a terminal there are no bars to line up with, so the record carries a timestamp instead.
pub fn install(level: LevelFilter) -> Result<(), log::SetLoggerError> {
    let attached = std::io::stderr().is_terminal();
    let target = match attached && level > LevelFilter::Error {
        true => ProgressDrawTarget::stderr(),
        false => ProgressDrawTarget::hidden(),
    };
    let _ = BARS.set(MultiProgress::with_draw_target(target));

    let mut builder = Builder::new();
    builder.filter_level(level);
    if let Ok(filters) = std::env::var("RUST_LOG") {
        builder.parse_filters(&filters);
    }
    builder.format(move |out, record| match attached {
        true => {
            let tint = level_tint(record.level());
            let track = tinted(palette::TRACK);
            writeln!(
                out,
                "{}{}{} {}{:<width$}{} {}",
                track.render(),
                stamp(),
                track.render_reset(),
                tint.render(),
                gutter(record),
                tint.render_reset(),
                record.args(),
                width = GUTTER
            )
        }
        false => {
            let stamp = out.timestamp();
            writeln!(
                out,
                "{stamp} {:<5} {} {}",
                record.level().as_str().to_ascii_lowercase(),
                record.target(),
                record.args()
            )
        }
    });

    let logger = builder.build();
    let filter = logger.filter();
    log::set_boxed_logger(Box::new(Through(logger)))?;
    log::set_max_level(filter);
    Ok(())
}

pub fn counted_template(stage: Stage) -> String {
    let colour = stage.colour();
    let track = palette::TRACK;
    format!(
        "{{stamp:.{track}}} {{prefix:<{GUTTER}.{colour}}} ▕{{bar:{BAR}.{colour}/{track}}}▏ {{human_pos:>7}}/{{human_len:<7}} {{eta:>4.{track}}}  {{msg}}"
    )
}

pub fn stepped_template(stage: Stage) -> String {
    let colour = stage.colour();
    let track = palette::TRACK;
    format!(
        "{{stamp:.{track}}} {{prefix:<{GUTTER}.{colour}}} ▕{{bar:{BAR}.{colour}/{track}}}▏ {{elapsed:>4.{track}}} left {{eta:<4.{track}}}  {{msg}}"
    )
}

pub fn spinning_template(stage: Stage) -> String {
    let colour = stage.colour();
    let track = palette::TRACK;
    format!(
        "{{stamp:.{track}}} {{prefix:<{GUTTER}.{colour}}} {{spinner:.{colour}}} {{elapsed:>4.{track}}}  {{msg}}"
    )
}

struct Hold {
    bar: ProgressBar,
    label: String,
}

enum Role {
    Alone,
    Holding,
    Relaying(Stage),
}

pub struct Bar {
    bar: ProgressBar,
    role: Role,
}

impl Bar {
    pub fn inc(&self, delta: u64) {
        self.bar.inc(delta);
        self.relay();
    }

    pub fn set_length(&self, length: u64) {
        self.bar.set_length(length);
        self.relay();
    }

    pub fn set_message(&self, message: impl Into<Cow<'static, str>>) {
        let message = message.into();
        if let Role::Holding = self.role
            && let Some(hold) = held().as_mut()
        {
            hold.label = message.to_string();
        }
        self.bar.set_message(message);
        self.relay();
    }

    pub fn message(&self) -> String {
        self.bar.message()
    }

    fn relay(&self) {
        let Role::Relaying(stage) = self.role else {
            return;
        };
        // A pane is often 100 columns wide, so a counted call shows its count and not its message.
        if let Some(hold) = held().as_ref() {
            let progress = match self.bar.length() {
                Some(length) => format!("{}/{length}", self.bar.position()),
                None => self.bar.message(),
            };
            let inner = format!("{} {progress}", stage.name().to_lowercase());
            hold.bar.set_message(match hold.label.is_empty() {
                true => inner.trim_end().to_string(),
                false => format!("{}, {}", hold.label, inner.trim_end()),
            });
        }
    }
}

fn held() -> MutexGuard<'static, Option<Hold>> {
    HOLD.lock().unwrap_or_else(PoisonError::into_inner)
}

fn made(stage: Stage, bar: ProgressBar, template: String) -> Bar {
    if held().is_some() {
        bar.set_draw_target(ProgressDrawTarget::hidden());
        return Bar {
            bar,
            role: Role::Relaying(stage),
        };
    }
    let bar = bars().add(bar);
    bar.set_style(styled(&template));
    bar.set_prefix(stage.name());
    bar.enable_steady_tick(TICK);
    Bar {
        bar,
        role: Role::Alone,
    }
}

pub fn counted(stage: Stage, total: u64) -> Bar {
    made(stage, ProgressBar::new(total), counted_template(stage))
}

// For work counted in units too fine to read, where the message names the coarse unit.
pub fn stepped(stage: Stage, total: u64) -> Bar {
    made(stage, ProgressBar::new(total), stepped_template(stage))
}

pub fn spinning(stage: Stage) -> Bar {
    made(stage, ProgressBar::new_spinner(), spinning_template(stage))
}

// A bar already reporting into another stays that way, since only one pass is held at a time.
pub fn hold(bar: Bar) -> Bar {
    let mut held = held();
    match (&bar.role, held.is_none()) {
        (Role::Alone, true) => {
            *held = Some(Hold {
                bar: bar.bar.clone(),
                label: bar.bar.message(),
            });
            Bar {
                role: Role::Holding,
                ..bar
            }
        }
        _ => bar,
    }
}

// A finished stage stays on screen as one line, so the next bar draws below it.
pub fn done(bar: &Bar) {
    match bar.role {
        Role::Relaying(_) => {
            if let Some(hold) = held().as_ref() {
                hold.bar.set_message(hold.label.clone());
            }
            bar.bar.finish_and_clear();
            return;
        }
        Role::Holding => {
            held().take();
        }
        Role::Alone => {}
    }
    let bar = &bar.bar;
    let elapsed = bar.elapsed();
    if elapsed >= KEPT && !bars().is_hidden() {
        let tint = tinted(palette::TRACK);
        let detail = match (bar.message(), bar.length()) {
            (message, _) if !message.is_empty() => message,
            (_, Some(length)) => format!("{}/{length}", bar.position()),
            _ => String::new(),
        };
        bars().suspend(|| {
            eprintln!(
                "{}{}{} {:<GUTTER$} ✓ {}{:>7}  {detail}{}",
                tint.render(),
                stamp(),
                tint.render_reset(),
                bar.prefix(),
                tint.render(),
                clock(elapsed),
                tint.render_reset(),
            )
        });
    }
    bar.finish_and_clear();
}

fn clock(elapsed: Duration) -> String {
    let seconds = elapsed.as_secs();
    match seconds {
        0..60 => format!("{seconds}s"),
        60..3600 => format!("{}m {:02}s", seconds / 60, seconds % 60),
        _ => format!("{}h {:02}m", seconds / 3600, seconds % 3600 / 60),
    }
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
