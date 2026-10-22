use std::path::Path;

use rosella::external::memory::cgroup_room;

fn write(root: &Path, at: &str, files: &[(&str, &str)]) {
    let directory = root.join(at);
    std::fs::create_dir_all(&directory).unwrap();
    for (name, text) in files {
        std::fs::write(directory.join(name), text).unwrap();
    }
}

#[test]
fn the_tightest_cgroup_above_the_job_bounds_its_room_less_reclaimable_cache() {
    let unified = tempfile::tempdir().unwrap();
    write(unified.path(), "", &[("memory.current", "900000000000")]);
    write(
        unified.path(),
        "pbs.slice/job.1",
        &[
            ("memory.max", "64000000000"),
            ("memory.current", "10000000000"),
            ("memory.stat", "anon 6000000000\ninactive_file 4000000000\n"),
        ],
    );
    write(
        unified.path(),
        "pbs.slice/job.1/task",
        &[("memory.max", "max"), ("memory.current", "9000000000")],
    );
    let split = tempfile::tempdir().unwrap();
    write(
        split.path(),
        "memory",
        &[
            ("memory.limit_in_bytes", "9223372036854771712"),
            ("memory.usage_in_bytes", "500000000000"),
        ],
    );
    write(
        split.path(),
        "memory/pbs/job.2",
        &[
            ("memory.limit_in_bytes", "8000000000"),
            ("memory.usage_in_bytes", "7000000000"),
            (
                "memory.stat",
                "inactive_file 5\ntotal_inactive_file 1000000000\n",
            ),
        ],
    );
    for (membership, root, room) in [
        (
            "0::/pbs.slice/job.1/task\n",
            unified.path(),
            Some(58_000_000_000),
        ),
        ("0::/\n", unified.path(), None),
        (
            "12:cpu,cpuacct:/pbs/job.2\n4:memory:/pbs/job.2\n0::/user.slice\n",
            split.path(),
            Some(2_000_000_000),
        ),
    ] {
        assert_eq!(cgroup_room(membership, root), room, "{membership}");
    }
}
