use std::{
    fs, thread,
    time::{Duration, Instant},
};

use rosella::timing::{self, TIMINGS_FILE};

struct Row {
    stage: String,
    calls: String,
    seconds: f64,
    percent: f64,
    cpu_seconds: f64,
    peak_mb: f64,
}

fn read_report(path: &std::path::Path) -> Vec<Row> {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    assert_eq!(
        lines.next().unwrap(),
        "stage\tcalls\tseconds\tpercent\tcpu_seconds\tpeak_rss_mb"
    );
    lines
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            assert_eq!(fields.len(), 6, "malformed row: {line}");
            Row {
                stage: fields[0].to_string(),
                calls: fields[1].to_string(),
                seconds: fields[2].parse().unwrap(),
                percent: fields[3].parse().unwrap(),
                cpu_seconds: fields[4].parse().unwrap(),
                peak_mb: fields[5].parse().unwrap(),
            }
        })
        .collect()
}

const HEAVY_MB: usize = 256;

// Mapped directly, since the allocator may hold a freed block and the next stage would see it.
fn touch_and_release(megabytes: usize) {
    let bytes = megabytes << 20;
    unsafe {
        let at = libc::mmap(
            std::ptr::null_mut(),
            bytes,
            libc::PROT_READ | libc::PROT_WRITE,
            libc::MAP_ANON | libc::MAP_PRIVATE,
            -1,
            0,
        );
        assert_ne!(at, libc::MAP_FAILED);
        std::ptr::write_bytes(at.cast::<u8>(), 1, bytes);
        libc::munmap(at, bytes);
    }
}

fn find<'a>(rows: &'a [Row], stage: &str) -> &'a Row {
    rows.iter()
        .find(|row| row.stage == stage)
        .unwrap_or_else(|| panic!("no {stage} row"))
}

/// One test, because the accumulator is process global and parallel tests would see each
/// other's scopes.
#[test]
fn repeated_scopes_accumulate_and_the_report_accounts_for_the_whole_run() {
    timing::start();

    for _ in 0..2 {
        let _timer = timing::scope("alpha");
        thread::sleep(Duration::from_millis(40));
    }
    {
        let _timer = timing::scope("beta");
        thread::sleep(Duration::from_millis(10));
    }
    {
        let _timer = timing::scope("spin");
        let until = Instant::now() + Duration::from_millis(60);
        while Instant::now() < until {
            std::hint::black_box(0);
        }
    }
    {
        let _timer = timing::scope("heavy");
        touch_and_release(HEAVY_MB);
    }
    {
        let _timer = timing::scope("light");
        thread::sleep(Duration::from_millis(5));
    }
    thread::sleep(Duration::from_millis(30));

    let directory = std::env::temp_dir().join("rosella_timing_test");
    fs::create_dir_all(&directory).unwrap();
    let path = directory.join(TIMINGS_FILE);
    timing::report(&path).unwrap();

    let rows = read_report(&path);

    let alpha = find(&rows, "alpha");
    assert_eq!(
        alpha.calls, "2",
        "two scopes of one name collapse to one row"
    );
    assert!(
        alpha.seconds >= 0.08,
        "alpha only counted {}s",
        alpha.seconds
    );

    let beta = find(&rows, "beta");
    assert_eq!(beta.calls, "1");

    let alpha_position = rows.iter().position(|row| row.stage == "alpha").unwrap();
    let beta_position = rows.iter().position(|row| row.stage == "beta").unwrap();
    assert!(alpha_position < beta_position, "rows are not slowest first");

    let spin = find(&rows, "spin");
    assert!(
        spin.cpu_seconds >= 0.04,
        "a busy scope was charged {}s of cpu",
        spin.cpu_seconds
    );
    assert!(
        alpha.cpu_seconds < 0.02,
        "a sleeping scope was charged {}s of cpu",
        alpha.cpu_seconds
    );

    let (heavy, light) = (find(&rows, "heavy"), find(&rows, "light"));
    assert!(
        heavy.peak_mb - light.peak_mb > HEAVY_MB as f64 / 2.0,
        "a stage after a {HEAVY_MB} MB one peaked at {} MB against its {} MB",
        light.peak_mb,
        heavy.peak_mb
    );

    let unaccounted = find(&rows, "unaccounted");
    assert!(unaccounted.calls.is_empty());
    assert!(
        unaccounted.seconds >= 0.03,
        "time outside any scope went missing: {}s",
        unaccounted.seconds
    );

    let total = find(&rows, "total");
    assert_eq!(total.percent, 100.0);
    assert!(total.peak_mb >= heavy.peak_mb);
    let partitioned = rows
        .iter()
        .filter(|row| row.stage != "total")
        .map(|row| row.seconds)
        .sum::<f64>();
    assert!(
        (partitioned - total.seconds).abs() < 0.01,
        "stages plus unaccounted was {partitioned}s against a total of {}s",
        total.seconds
    );
    let charged = rows
        .iter()
        .filter(|row| row.stage != "total")
        .map(|row| row.cpu_seconds)
        .sum::<f64>();
    assert!(
        (charged - total.cpu_seconds).abs() < 0.01,
        "stage cpu plus unaccounted was {charged}s against a total of {}s",
        total.cpu_seconds
    );

    fs::remove_file(&path).unwrap();
}

#[test]
fn a_sample_keeps_one_row_per_stage_and_a_rerun_replaces_only_its_own() {
    let directory = tempfile::tempdir().unwrap();
    let path = timing::stages_path(directory.path(), "s1");
    let watch = timing::Stopwatch::start();
    thread::sleep(Duration::from_millis(30));
    for stage in ["mapping", "depth", "mapping"] {
        timing::record(&path, &watch.row(stage)).unwrap();
    }
    let text = fs::read_to_string(&path).unwrap();
    let rows = text.lines().collect::<Vec<_>>();
    assert_eq!(rows[0], timing::STAGE_COLUMNS);
    let stages = rows[1..]
        .iter()
        .map(|row| row.split('\t').next().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(stages, ["depth", "mapping"]);
    let fields = rows[1]
        .split('\t')
        .skip(1)
        .map(|field| field.parse::<f64>().unwrap());
    let [seconds, user, sys, peak] = fields.collect::<Vec<_>>()[..] else {
        panic!("{}", rows[1]);
    };
    assert!(
        seconds >= 0.03 && user >= 0.0 && sys >= 0.0 && peak > 0.0,
        "{}",
        rows[1]
    );
}
