mod lock;
mod reach;
mod tokens;
mod units;

use std::collections::{BTreeSet, HashSet};
use std::path::{Path, PathBuf};

use units::{Place, Unit};

pub struct Digests {
    pub host: String,
    pub device: String,
}

// Names resolve without types, so a reference reaches every item it could mean. Reaching too
// much retrains for nothing, and reaching too little would serve a stale space.
pub fn digests(
    root: &Path,
    entries: &[(&str, &str)],
    outside: &[(&str, &str)],
    device: &Path,
    lock: &str,
    manifest: &str,
) -> Digests {
    let files = rust_files(root);
    let (mut units, mut whole) = (Vec::new(), Vec::new());
    for (file, module) in &files {
        let source = std::fs::read_to_string(file).unwrap();
        let place = Place { file, module };
        match syn::parse_file(&source) {
            Ok(parsed) => units::walk(&parsed.items, &place, &mut units),
            Err(_) => {
                whole.push(units.len());
                units.push(Unit {
                    file: file.clone(),
                    module: module.clone(),
                    text: source,
                    ..Unit::default()
                });
            }
        }
    }
    let reached = reach::reach(&units, &files, entries, outside, whole);
    let coded = reached
        .iter()
        .map(|at| &units[*at])
        .filter(|unit| !unit.shape && unit.import.is_none())
        .map(|unit| unit.file.as_path())
        .collect::<HashSet<_>>();
    let [host, device] = [false, true].map(|on_device| {
        let group = reached
            .iter()
            .map(|at| &units[*at])
            .filter(|unit| !unit.shape || coded.contains(unit.file.as_path()))
            .filter(|unit| unit.file.starts_with(device) == on_device)
            .collect::<Vec<_>>();
        digest(&group, lock, manifest)
    });
    Digests { host, device }
}

fn rust_files(root: &Path) -> Vec<(PathBuf, Vec<String>)> {
    let mut found = Vec::new();
    let mut pending = vec![root.to_path_buf()];
    while let Some(path) = pending.pop() {
        if path.is_dir() {
            pending.extend(std::fs::read_dir(&path).unwrap().map(|e| e.unwrap().path()));
        } else if path.extension().is_some_and(|extension| extension == "rs") {
            let mut module = path
                .strip_prefix(root)
                .unwrap()
                .with_extension("")
                .iter()
                .map(|part| part.to_string_lossy().into_owned())
                .collect::<Vec<_>>();
            if matches!(module.as_slice(), [only] if only == "lib")
                || module.last().is_some_and(|last| last == "mod")
            {
                module.pop();
            }
            found.push((path, module));
        }
    }
    found.sort();
    found
}

fn digest(group: &[&Unit], lock: &str, manifest: &str) -> String {
    let mut parts = BTreeSet::new();
    let mut crates = HashSet::new();
    for unit in group {
        parts.insert(unit.text.clone().into_bytes());
        for file in &unit.includes {
            let mut part = file.to_string_lossy().into_owned().into_bytes();
            part.push(0);
            part.extend(std::fs::read(file).unwrap());
            parts.insert(part);
        }
        let firsts = unit
            .refs
            .paths
            .iter()
            .chain(&unit.import)
            .flat_map(|path| path.first());
        crates.extend(firsts.chain(&unit.refs.idents).cloned());
    }
    parts.extend(
        lock::locked(lock, manifest, &crates)
            .into_iter()
            .map(String::into_bytes),
    );
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for part in parts {
        for byte in part.into_iter().chain([0]) {
            hash = (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
        }
    }
    format!("{hash:016x}")
}
