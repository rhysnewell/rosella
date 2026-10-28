use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
    sync::{Mutex, MutexGuard, OnceLock, PoisonError},
    time::{Duration, Instant},
};

use anyhow::Result;
use log::debug;

pub const TIMINGS_FILE: &str = "timings.tsv";

#[derive(Default)]
struct Stage {
    calls: u32,
    wall: Duration,
    cpu: Duration,
}

struct Ledger {
    stages: BTreeMap<&'static str, Stage>,
    open: BTreeMap<&'static str, u32>,
    cpu_at: Duration,
}

// Process CPU is read at every scope boundary and shared over the scopes open since the last
// one, so stages running side by side split the cores rather than each claiming all of them.
impl Ledger {
    fn settle(&mut self) {
        let now = process_cpu();
        let spent = now.saturating_sub(self.cpu_at);
        self.cpu_at = now;
        let open = self.open.values().sum::<u32>();
        if open == 0 {
            return;
        }
        for (name, scopes) in &self.open {
            self.stages.entry(name).or_default().cpu += spent * *scopes / open;
        }
    }
}

static LEDGER: Mutex<Ledger> = Mutex::new(Ledger {
    stages: BTreeMap::new(),
    open: BTreeMap::new(),
    cpu_at: Duration::ZERO,
});
static START: OnceLock<Instant> = OnceLock::new();

fn ledger() -> MutexGuard<'static, Ledger> {
    LEDGER.lock().unwrap_or_else(PoisonError::into_inner)
}

fn process_cpu() -> Duration {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::zeroed();
    if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } != 0 {
        return Duration::ZERO;
    }
    let usage = unsafe { usage.assume_init() };
    let time = |at: libc::timeval| Duration::new(at.tv_sec as u64, at.tv_usec as u32 * 1000);
    time(usage.ru_utime) + time(usage.ru_stime)
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
        let mut ledger = ledger();
        ledger.settle();
        if let Some(scopes) = ledger.open.get_mut(self.name) {
            *scopes -= 1;
            if *scopes == 0 {
                ledger.open.remove(self.name);
            }
        }
        let stage = ledger.stages.entry(self.name).or_default();
        stage.calls += 1;
        stage.wall += elapsed;
    }
}

// Stages must not nest inside one another. They may run concurrently, in which case the
// accumulated time is thread seconds and `report` divides by that rather than by the wall.
pub fn scope(name: &'static str) -> Scope {
    let mut ledger = ledger();
    ledger.settle();
    *ledger.open.entry(name).or_default() += 1;
    Scope {
        name,
        started: Instant::now(),
    }
}

pub fn report(path: impl AsRef<Path>) -> Result<()> {
    let total = START.get().map(Instant::elapsed).unwrap_or_default();
    let mut ledger = ledger();
    ledger.settle();
    let total_cpu = ledger.cpu_at;
    let measured = ledger
        .stages
        .values()
        .map(|stage| stage.wall)
        .sum::<Duration>();
    let measured_cpu = ledger
        .stages
        .values()
        .map(|stage| stage.cpu)
        .sum::<Duration>();

    let mut rows = ledger
        .stages
        .iter()
        .map(|(name, stage)| (*name, stage.calls.to_string(), stage.wall, stage.cpu))
        .collect::<Vec<_>>();
    rows.sort_by(|left, right| right.2.cmp(&left.2).then(left.0.cmp(right.0)));
    rows.push((
        "unaccounted",
        String::new(),
        total.saturating_sub(measured),
        total_cpu.saturating_sub(measured_cpu),
    ));
    rows.push(("total", String::new(), total, total_cpu));
    let share_of = measured.max(total);

    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "stage\tcalls\tseconds\tpercent\tcpu_seconds")?;
    for (name, calls, elapsed, cpu) in &rows {
        let percent = if share_of.is_zero() {
            0.0
        } else {
            elapsed.as_secs_f64() / share_of.as_secs_f64() * 100.0
        };
        writeln!(
            out,
            "{}\t{}\t{:.3}\t{:.1}\t{:.3}",
            name,
            calls,
            elapsed.as_secs_f64(),
            percent,
            cpu.as_secs_f64()
        )?;
        debug!(
            "{} {:.3}s {:.1}% {:.3}s cpu{}",
            name,
            elapsed.as_secs_f64(),
            percent,
            cpu.as_secs_f64(),
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
