use std::collections::{BTreeSet, HashSet};

struct Package<'a> {
    name: &'a str,
    version: &'a str,
    text: String,
    dependencies: Vec<&'a str>,
}

fn packages(lock: &str) -> Vec<Package<'_>> {
    lock.split("[[package]]")
        .skip(1)
        .map(|block| {
            let field = |key: &str| {
                block
                    .lines()
                    .find_map(|line| line.strip_prefix(key))
                    .map(|value| value.trim_matches('"'))
                    .unwrap_or_default()
            };
            let dependencies = block
                .split_once("dependencies = [")
                .map(|(_, rest)| rest.split(']').next().unwrap_or_default())
                .unwrap_or_default()
                .lines()
                .map(|line| line.trim().trim_end_matches(',').trim_matches('"'))
                .filter(|line| !line.is_empty())
                .collect();
            Package {
                name: field("name = "),
                version: field("version = "),
                text: ["name = ", "version = ", "source = ", "checksum = "]
                    .map(field)
                    .join(" "),
                dependencies,
            }
        })
        .collect()
}

fn named(name: &str, crates: &HashSet<String>) -> bool {
    crates.contains(&name.trim().replace('-', "_"))
}

// A crate the reached code never names cannot move the weights, so only the named crates, what
// they pull in, and how the manifest asks for them join the key.
pub fn locked(lock: &str, manifest: &str, crates: &HashSet<String>) -> Vec<String> {
    let packages = packages(lock);
    let mut pending = packages
        .iter()
        .enumerate()
        .filter(|(_, package)| named(package.name, crates))
        .map(|(at, _)| at)
        .collect::<Vec<_>>();
    let mut kept = BTreeSet::new();
    while let Some(at) = pending.pop() {
        if !kept.insert(at) {
            continue;
        }
        for wanted in &packages[at].dependencies {
            let mut words = wanted.split(' ');
            let (name, version) = (words.next().unwrap_or_default(), words.next());
            pending.extend(
                packages
                    .iter()
                    .enumerate()
                    .filter(|(_, p)| p.name == name && version.is_none_or(|v| p.version == v))
                    .map(|(at, _)| at),
            );
        }
    }
    let asked = manifest
        .lines()
        .map(str::trim)
        .filter(|line| {
            line.split_once('=')
                .is_some_and(|(name, _)| named(name, crates))
        })
        .map(str::to_string);
    kept.into_iter()
        .map(|at| packages[at].text.clone())
        .chain(asked)
        .collect()
}
