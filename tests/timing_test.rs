use std::{fs, thread, time::Duration};

use rosella::timing::{self, TIMINGS_FILE};

struct Row {
    stage: String,
    calls: String,
    seconds: f64,
    percent: f64,
}

fn read_report(path: &std::path::Path) -> Vec<Row> {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    assert_eq!(lines.next().unwrap(), "stage\tcalls\tseconds\tpercent");
    lines
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            assert_eq!(fields.len(), 4, "malformed row: {line}");
            Row {
                stage: fields[0].to_string(),
                calls: fields[1].to_string(),
                seconds: fields[2].parse().unwrap(),
                percent: fields[3].parse().unwrap(),
            }
        })
        .collect()
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

    let unaccounted = find(&rows, "unaccounted");
    assert!(unaccounted.calls.is_empty());
    assert!(
        unaccounted.seconds >= 0.03,
        "time outside any scope went missing: {}s",
        unaccounted.seconds
    );

    let total = find(&rows, "total");
    assert_eq!(total.percent, 100.0);
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

    fs::remove_file(&path).unwrap();
}
