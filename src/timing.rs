use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
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
    let (user, sys) = user_and_sys();
    user + sys
}

fn rusage() -> Option<libc::rusage> {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::zeroed();
    (unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } == 0)
        .then(|| unsafe { usage.assume_init() })
}

fn user_and_sys() -> (Duration, Duration) {
    let time = |at: libc::timeval| Duration::new(at.tv_sec as u64, at.tv_usec as u32 * 1000);
    rusage().map_or((Duration::ZERO, Duration::ZERO), |usage| {
        (time(usage.ru_utime), time(usage.ru_stime))
    })
}

// Linux resets the resident high-water mark on request, so a stage reports its own peak.
// Elsewhere the peak is the run's so far.
fn forget_peak() {
    let _ = std::fs::write("/proc/self/clear_refs", "5");
}

fn peak_bytes() -> u64 {
    let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    let high_water = status
        .lines()
        .find_map(|line| line.strip_prefix("VmHWM:"))
        .and_then(|kb| kb.trim().trim_end_matches("kB").trim().parse::<u64>().ok());
    match high_water {
        Some(kb) => kb * 1024,
        None => rusage().map_or(0, |usage| usage.ru_maxrss as u64),
    }
}

pub const STAGE_COLUMNS: &str = "stage\tseconds\tuser_seconds\tsys_seconds\tpeak_rss_mb";

pub struct Stopwatch {
    started: Instant,
    user: Duration,
    sys: Duration,
}

impl Stopwatch {
    pub fn start() -> Self {
        forget_peak();
        let (user, sys) = user_and_sys();
        Self {
            started: Instant::now(),
            user,
            sys,
        }
    }

    pub fn row(&self, stage: &str) -> String {
        let (user, sys) = user_and_sys();
        format!(
            "{stage}\t{:.3}\t{:.3}\t{:.3}\t{:.1}",
            self.started.elapsed().as_secs_f64(),
            user.saturating_sub(self.user).as_secs_f64(),
            sys.saturating_sub(self.sys).as_secs_f64(),
            peak_bytes() as f64 / (1024.0 * 1024.0)
        )
    }
}

pub fn stages_path(directory: &Path, sample: &str) -> PathBuf {
    directory.join(format!("{sample}.timings.tsv"))
}

// A rerun replaces only the stage it repeats, so depth remade from a kept BAM keeps the
// mapping's row.
pub fn record(path: &Path, row: &str) -> Result<()> {
    let stage = row.split('\t').next().unwrap_or_default();
    let kept = match std::fs::read_to_string(path) {
        Ok(text) => text,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => String::new(),
        Err(error) => return Err(error.into()),
    };
    let mut rows = kept
        .lines()
        .skip(1)
        .filter(|line| line.split('\t').next() != Some(stage))
        .chain([row])
        .collect::<Vec<_>>();
    rows.insert(0, STAGE_COLUMNS);
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    crate::report_sink::write_atomically(path, |mut file| {
        Ok(file.write_all(format!("{}\n", rows.join("\n")).as_bytes())?)
    })
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
