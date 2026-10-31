use std::path::Path;

use sysinfo::System;

pub fn free_memory() -> u64 {
    let host = reclaimable();
    let job = std::fs::read_to_string("/proc/self/cgroup")
        .ok()
        .and_then(|membership| cgroup_room(&membership, Path::new("/sys/fs/cgroup")));
    job.map_or(host, |job| job.min(host))
}

// sysinfo takes the compressor's pages off free and inactive pages, though it holds neither.
#[cfg(target_os = "macos")]
fn reclaimable() -> u64 {
    let mut stat = unsafe { std::mem::zeroed::<libc::vm_statistics64>() };
    let mut count = libc::HOST_VM_INFO64_COUNT;
    #[allow(deprecated)]
    let host = unsafe { libc::mach_host_self() };
    let status = unsafe {
        libc::host_statistics64(
            host,
            libc::HOST_VM_INFO64,
            &mut stat as *mut libc::vm_statistics64 as libc::host_info64_t,
            &mut count,
        )
    };
    if status != libc::KERN_SUCCESS {
        return available();
    }
    let page = unsafe { libc::vm_page_size } as u64;
    let pages = stat.free_count as u64 + stat.inactive_count as u64 + stat.purgeable_count as u64;
    pages * page
}

#[cfg(not(target_os = "macos"))]
fn reclaimable() -> u64 {
    available()
}

fn available() -> u64 {
    let mut system = System::new();
    system.refresh_memory();
    system.available_memory()
}

// A batch job's limit can sit on any cgroup above the process, and the host knows none of them.
// A job whose page cache fills its limit reads as no room, which would split the index per contig.
pub fn cgroup_room(membership: &str, root: &Path) -> Option<u64> {
    let (base, limit, usage, cache) = membership.lines().find_map(|line| {
        let mut fields = line.splitn(3, ':');
        let (id, controllers, path) = (fields.next()?, fields.next()?, fields.next()?);
        let path = path.trim_start_matches('/');
        match (id, controllers) {
            ("0", "") => Some((
                root.join(path),
                "memory.max",
                "memory.current",
                "inactive_file",
            )),
            (_, controllers) if controllers.split(',').any(|name| name == "memory") => Some((
                root.join("memory").join(path),
                "memory.limit_in_bytes",
                "memory.usage_in_bytes",
                "total_inactive_file",
            )),
            _ => None,
        }
    })?;
    let read = |directory: &Path, name: &str| -> Option<u64> {
        std::fs::read_to_string(directory.join(name))
            .ok()?
            .trim()
            .parse()
            .ok()
    };
    let top = base.ancestors().take_while(|at| at.starts_with(root));
    top.filter_map(|at| {
        let limit = read(at, limit)?;
        let used = read(at, usage)?;
        let stat = std::fs::read_to_string(at.join("memory.stat")).unwrap_or_default();
        let cached = stat
            .lines()
            .find_map(|line| line.strip_prefix(cache)?.trim().parse::<u64>().ok())
            .unwrap_or(0);
        Some(limit.saturating_sub(used.saturating_sub(cached)))
    })
    .min()
    .filter(|room| *room > 0)
}
