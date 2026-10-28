use std::path::{Path, PathBuf};
use std::process::Command;

const MINIMAP2: &str = "vendor/minimap2";

fn git(args: &[&str]) -> Option<String> {
    let output = Command::new("git").args(args).output().ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_string())
}

fn stamped(commit: &str, dirty: bool) -> String {
    match dirty {
        true => format!("{commit}-dirty"),
        false => commit.to_string(),
    }
}

// A crates.io build has no .git, so it reads the commit cargo packed the crate from.
fn packaged() -> Option<String> {
    let info = std::fs::read_to_string(".cargo_vcs_info.json").ok()?;
    let commit = info.split("\"sha1\": \"").nth(1)?.get(..7)?;
    Some(stamped(commit, info.contains("\"dirty\": true")))
}

fn main() {
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/logs/HEAD");
    println!("cargo:rerun-if-changed=.git/index");
    println!("cargo:rerun-if-changed=src");

    let stamp = packaged().unwrap_or_else(|| {
        let commit =
            git(&["rev-parse", "--short", "HEAD"]).unwrap_or_else(|| "unknown".to_string());
        let dirty = git(&["status", "--porcelain", "--untracked-files=no"])
            .is_some_and(|out| !out.is_empty());
        stamped(&commit, dirty)
    });
    println!("cargo:rustc-env=ROSELLA_BUILD_COMMIT={stamp}");
    println!("cargo:rerun-if-changed=Cargo.lock");
    println!(
        "cargo:rustc-env=ROSELLA_TRAINING_DIGEST={}",
        training_digest()
    );
    println!("cargo:rustc-env=ROSELLA_GENE_CALLER={}", gene_caller());
    minimap2();
}

// The learned space cache is keyed on the code that shapes training, so a change anywhere else
// reuses a kept space and a change to training never reads one trained before it.
const TRAINING: [&str; 11] = [
    "src/learned",
    "src/kmers",
    "src/coverage",
    "src/tables.rs",
    "src/recover/inputs.rs",
    "src/recover/recover_engine/learned.rs",
    "src/digest.rs",
    "src/rows.rs",
    "src/defaults.rs",
    "src/lib.rs",
    "Cargo.lock",
];

fn training_digest() -> String {
    let mut files = Vec::new();
    for root in TRAINING {
        files_under(Path::new(root), &mut files);
    }
    files.sort();
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for file in files {
        let name = file.to_string_lossy().into_owned().into_bytes();
        let body = std::fs::read(&file).unwrap();
        for byte in name.into_iter().chain([0]).chain(body) {
            hash = (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
        }
    }
    format!("{hash:016x}")
}

fn files_under(path: &Path, found: &mut Vec<PathBuf>) {
    if path.is_dir() {
        for entry in std::fs::read_dir(path).unwrap() {
            files_under(&entry.unwrap().path(), found);
        }
    } else {
        found.push(path.to_path_buf());
    }
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

// The forced header sends minimap2's SAM to a pipe rosella reads, so the vendored tree stays as
// released. Its alignment kernels build once per instruction set, as its Makefile does.
fn minimap2() {
    println!("cargo:rerun-if-changed={MINIMAP2}");
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let arch = std::env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    let root = Path::new(MINIMAP2);
    let header = Path::new("src/external/minimap2_out.h");
    let zlib = std::env::var("DEP_Z_INCLUDE").ok();
    let configure = |build: &mut cc::Build| {
        build
            .include(root)
            .define("HAVE_KALLOC", None)
            .flag("-include")
            .flag(header)
            .opt_level(2)
            .debug(false)
            .warnings(false)
            .flag_if_supported("-pthread")
            .cargo_metadata(false);
        if let Some(zlib) = &zlib {
            build.include(zlib);
        }
        if arch == "aarch64" {
            build
                .include(root.join("sse2neon"))
                .define("_FILE_OFFSET_BITS", "64")
                .flag("-fsigned-char");
        }
    };
    let kernels: &[(&str, &[&str])] = match arch.as_str() {
        "x86_64" => &[
            ("sse41", &["-msse4.1", "-DKSW_CPU_DISPATCH"]),
            (
                "sse2",
                &[
                    "-msse2",
                    "-mno-sse4.1",
                    "-DKSW_CPU_DISPATCH",
                    "-DKSW_SSE2_ONLY",
                ],
            ),
        ],
        "aarch64" => &[("neon", &["-DKSW_SSE2_ONLY", "-D__SSE2__"])],
        other => panic!("rosella builds minimap2 for x86_64 and aarch64, not {other}"),
    };
    let stamp = out.join("minimap2.stamp");
    let print = format!(
        "{arch}\n{}\n{}",
        fingerprint(root),
        std::fs::read_to_string(header).unwrap()
    );
    if std::fs::read_to_string(&stamp).ok().as_deref() != Some(print.as_str()) {
        let simd =
            ["ksw2_extz2_sse.c", "ksw2_extd2_sse.c", "ksw2_exts2_sse.c"].map(|c| root.join(c));
        let dispatch = root.join("ksw2_dispatch.c");
        let mut library = cc::Build::new();
        configure(&mut library);
        for (name, flags) in kernels {
            let mut kernel = cc::Build::new();
            configure(&mut kernel);
            kernel.out_dir(out.join(name)).files(&simd);
            for flag in *flags {
                kernel.flag(flag);
            }
            if *name == "sse41" {
                kernel.file(&dispatch);
            }
            library.objects(kernel.compile_intermediates());
        }
        library
            .files(
                sources(root)
                    .into_iter()
                    .filter(|path| !simd.contains(path) && *path != dispatch),
            )
            .define("main", "rosella_minimap2_main")
            .compile("minimap2");
        std::fs::write(&stamp, print).unwrap();
    }
    cc::Build::new()
        .file("src/external/minimap2_run.c")
        .opt_level(2)
        .debug(false)
        .cargo_metadata(false)
        .compile("rosella_minimap2");
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=rosella_minimap2");
    println!("cargo:rustc-link-lib=static=minimap2");
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
