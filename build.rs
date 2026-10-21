use std::path::{Path, PathBuf};
use std::process::Command;

const HMMER: &str = "vendor/hmmer";

fn git(args: &[&str]) -> Option<String> {
    let output = Command::new("git").args(args).output().ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_string())
}

fn main() {
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/logs/HEAD");
    println!("cargo:rerun-if-changed=.git/index");
    println!("cargo:rerun-if-changed=src");

    let commit = git(&["rev-parse", "--short", "HEAD"]).unwrap_or_else(|| "unknown".to_string());
    let dirty =
        git(&["status", "--porcelain", "--untracked-files=no"]).is_some_and(|out| !out.is_empty());
    let stamp = if dirty {
        format!("{commit}-dirty")
    } else {
        commit
    };
    println!("cargo:rustc-env=ROSELLA_BUILD_COMMIT={stamp}");
    println!("cargo:rerun-if-changed=Cargo.lock");
    println!("cargo:rustc-env=ROSELLA_GENE_CALLER={}", gene_caller());
    hmmer();
}

fn sources(directory: &Path) -> Vec<PathBuf> {
    let mut found = std::fs::read_dir(directory)
        .unwrap_or_else(|error| panic!("cannot read {}: {error}", directory.display()))
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|extension| extension == "c"))
        .collect::<Vec<_>>();
    found.sort();
    found
}

fn fingerprint(directory: &Path) -> String {
    let mut stack = vec![directory.to_path_buf()];
    let mut seen = Vec::new();
    while let Some(at) = stack.pop() {
        for entry in std::fs::read_dir(&at).unwrap() {
            let path = entry.unwrap().path();
            let meta = std::fs::metadata(&path).unwrap();
            match meta.is_dir() {
                true => stack.push(path),
                false => seen.push(format!(
                    "{} {} {:?}",
                    path.display(),
                    meta.len(),
                    meta.modified().ok()
                )),
            }
        }
    }
    seen.sort();
    seen.join("\n")
}

// This script reruns on every source edit, and HMMER takes most of a minute to compile, so its
// library is rebuilt only when the vendored tree changes.
fn hmmer() {
    println!("cargo:rerun-if-changed={HMMER}");
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let simd = match std::env::var("CARGO_CFG_TARGET_ARCH").unwrap().as_str() {
        "aarch64" => "impl_neon",
        "x86_64" => "impl_sse",
        other => panic!("rosella builds HMMER for x86_64 and aarch64, not {other}"),
    };
    let root = Path::new(HMMER);
    let includes = [
        root.join("easel"),
        root.join("src"),
        root.join("src").join(simd),
    ];
    let configure = |build: &mut cc::Build| {
        build
            .includes(&includes)
            .define("HAVE_CONFIG_H", None)
            .opt_level(3)
            .debug(false)
            .warnings(false)
            .flag_if_supported("-pthread")
            .cargo_metadata(false);
    };
    let stamp = out.join("hmmer.stamp");
    let print = format!("{simd}\n{}", fingerprint(root));
    if std::fs::read_to_string(&stamp).ok().as_deref() != Some(print.as_str()) {
        let mut library = cc::Build::new();
        configure(&mut library);
        for directory in [
            root.join("easel"),
            root.join("src"),
            root.join("src").join(simd),
        ] {
            library.files(sources(&directory));
        }
        library.compile("hmmer");
        std::fs::write(&stamp, print).unwrap();
    }
    let mut search = cc::Build::new();
    configure(&mut search);
    search
        .file("src/external/hmmer_search.c")
        .compile("rosella_hmmer");
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=rosella_hmmer");
    println!("cargo:rustc-link-lib=static=hmmer");
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
