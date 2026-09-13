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
/// locked revision rather than rosella's own commit.
fn gene_caller() -> String {
    let Ok(lock) = std::fs::read_to_string("Cargo.lock") else {
        return "unknown".to_string();
    };
    let mut frugal = false;
    for line in lock.lines() {
        if line.starts_with("name = ") {
            frugal = line.contains("\"frugal\"");
        }
        if frugal && let Some(source) = line.strip_prefix("source = ") {
            return source.trim_matches('"').to_string();
        }
    }
    "unknown".to_string()
}
