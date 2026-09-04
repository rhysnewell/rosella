use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
    sync::{Mutex, OnceLock},
    time::{Duration, Instant},
};

use anyhow::Result;
use log::info;

pub const TIMINGS_FILE: &str = "timings.tsv";

static STAGES: OnceLock<Mutex<BTreeMap<&'static str, (u32, Duration)>>> = OnceLock::new();
static START: OnceLock<Instant> = OnceLock::new();

type Stages = Mutex<BTreeMap<&'static str, (u32, Duration)>>;

fn stages() -> &'static Stages {
    STAGES.get_or_init(|| Mutex::new(BTreeMap::new()))
}

pub fn start() {
    let _ = START.set(Instant::now());
}

pub struct Scope {
    name: &'static str,
    started: Instant,
}

impl Drop for Scope {
    fn drop(&mut self) {
        let elapsed = self.started.elapsed();
        let mut stages = stages().lock().unwrap();
        let entry = stages.entry(self.name).or_insert((0, Duration::ZERO));
        entry.0 += 1;
        entry.1 += elapsed;
    }
}

/// Stages must not nest inside one another. They may run concurrently, in which case the
/// accumulated time is thread seconds and `report` divides by that rather than by the wall.
pub fn scope(name: &'static str) -> Scope {
    Scope {
        name,
        started: Instant::now(),
    }
}

pub fn report(path: impl AsRef<Path>) -> Result<()> {
    let total = START.get().map(Instant::elapsed).unwrap_or_default();
    let stages = stages().lock().unwrap();
    let measured = stages
        .values()
        .map(|(_, elapsed)| *elapsed)
        .sum::<Duration>();

    let mut rows = stages
        .iter()
        .map(|(name, (calls, elapsed))| (*name, calls.to_string(), *elapsed))
        .collect::<Vec<_>>();
    rows.sort_by(|left, right| right.2.cmp(&left.2).then(left.0.cmp(right.0)));
    rows.push(("unaccounted", String::new(), total.saturating_sub(measured)));
    rows.push(("total", String::new(), total));
    let share_of = measured.max(total);

    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "stage\tcalls\tseconds\tpercent")?;
    for (name, calls, elapsed) in &rows {
        let percent = if share_of.is_zero() {
            0.0
        } else {
            elapsed.as_secs_f64() / share_of.as_secs_f64() * 100.0
        };
        writeln!(
            out,
            "{}\t{}\t{:.3}\t{:.1}",
            name,
            calls,
            elapsed.as_secs_f64(),
            percent
        )?;
        info!(
            "{} {:.3}s {:.1}%{}",
            name,
            elapsed.as_secs_f64(),
            percent,
            if calls.is_empty() {
                String::new()
            } else {
                format!(" over {} calls", calls)
            }
        );
    }
    out.flush()?;
    Ok(())
}
