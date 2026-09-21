use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    let output = Command::new("git").args(args).output().ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_string())
}

fn main() {
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/index");
    println!("cargo:rerun-if-changed=src");

    let commit = git(&["rev-parse", "--short", "HEAD"]).unwrap_or_else(|| "unknown".to_string());
    let dirty = git(&["status", "--porcelain"]).is_some_and(|out| !out.is_empty());
    let stamp = if dirty {
        format!("{commit}-dirty")
    } else {
        commit
    };
    println!("cargo:rustc-env=ROSELLA_BUILD_COMMIT={stamp}");
    println!("cargo:rerun-if-changed=Cargo.lock");
    println!("cargo:rustc-env=ROSELLA_GENE_CALLER={}", gene_caller());
}

/// The marker annotation only changes when the gene caller does, so the cache key reads the
/// locked frugal. A registry source is one string for every version, so it cannot carry that.
fn gene_caller() -> String {
    let Ok(lock) = std::fs::read_to_string("Cargo.lock") else {
        return "unknown".to_string();
    };
    let mut frugal = false;
    let mut version = String::new();
    let mut source = String::new();
    for line in lock.lines() {
        if line.starts_with("name = ") {
            if frugal {
                break;
            }
            frugal = line.contains("\"frugal\"");
            continue;
        }
        if !frugal {
            continue;
        }
        if let Some(v) = line.strip_prefix("version = ") {
            version = v.trim_matches('"').to_string();
        }
        if let Some(s) = line.strip_prefix("source = ") {
            source = s.trim_matches('"').to_string();
        }
        if let Some(c) = line.strip_prefix("checksum = ") {
            return format!("{version}+{}", c.trim_matches('"'));
        }
    }
    if version.is_empty() && source.is_empty() {
        return "unknown".to_string();
    }
    format!("{version}+{source}")
}
