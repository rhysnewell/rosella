use crate::training::digests;

const LIB: &str = "pub mod train;\npub mod other;\npub mod device;\n";

const TRAIN: &str = r#"use crate::other::helper;

pub struct Net {
    pub width: usize,
}

impl Net {
    pub fn step(&self) -> f32 {
        helper(self.width)
    }
}

impl Drop for Net {
    fn drop(&mut self) {}
}

pub fn trained() -> f32 {
    log::info!("training");
    let net = Net { width: 3 };
    ndarray::Array2::<f32>::zeros((1, 1));
    net.step()
}

pub fn unrelated() -> u32 {
    7
}
"#;

const OTHER: &str = r#"pub fn helper(width: usize) -> f32 {
    width as f32 * 0.5
}

pub fn report() -> String {
    String::new()
}

pub struct Report;

impl Report {
    pub fn step(&self) -> u8 {
        1
    }
}
"#;

const DEVICE: &str = "pub fn open() -> u8 {\n    2\n}\n";

const LOCK: &str = r#"[[package]]
name = "ndarray"
version = "0.17.0"
dependencies = [
 "num-traits",
]

[[package]]
name = "num-traits"
version = "0.2.19"

[[package]]
name = "tongs"
version = "0.2.0"
"#;

const MANIFEST: &str = "[dependencies]\nndarray = \"0.17\"\ntongs = \"0.2\"\n";

fn keyed(edit: Option<(&str, &str, &str)>) -> (String, String) {
    let root = tempfile::tempdir().unwrap();
    let src = root.path().join("src");
    std::fs::create_dir_all(src.join("device")).unwrap();
    let mut lock = LOCK.to_string();
    let mut manifest = MANIFEST.to_string();
    for (file, body) in [
        ("lib.rs", LIB),
        ("train.rs", TRAIN),
        ("other.rs", OTHER),
        ("device/mod.rs", DEVICE),
    ] {
        let body = match edit {
            Some((edited, from, to)) if edited == file => {
                assert!(body.contains(from), "{from} is not in {file}");
                body.replacen(from, to, 1)
            }
            _ => body.to_string(),
        };
        std::fs::write(src.join(file), body).unwrap();
    }
    if let Some(("Cargo.lock", from, to)) = edit {
        lock = lock.replacen(from, to, 1);
    }
    if let Some(("Cargo.toml", from, to)) = edit {
        manifest = manifest.replacen(from, to, 1);
    }
    let file = |name: &str| src.join(name).to_string_lossy().into_owned();
    let (train, device) = (file("train.rs"), file("device/mod.rs"));
    let found = digests(
        &src,
        &[(&train, "trained"), (&device, "open")],
        &[],
        &src.join("device"),
        &lock,
        &manifest,
    );
    (found.host, found.device)
}

#[test]
fn the_key_moves_only_with_what_training_reaches() {
    let base = keyed(None);
    let cases = [
        ("train.rs", "    7\n", "    8\n", false),
        ("train.rs", "\"training\"", "\"learning\"", false),
        (
            "train.rs",
            "pub fn trained",
            "// why\npub fn trained",
            false,
        ),
        ("other.rs", "String::new()", "String::from(\"x\")", false),
        ("other.rs", "        1\n", "        2\n", false),
        ("Cargo.lock", "0.2.0", "0.3.0", false),
        ("Cargo.toml", "tongs = \"0.2\"", "tongs = \"0.3\"", false),
        ("other.rs", "* 0.5", "* 0.25", true),
        (
            "train.rs",
            "helper(self.width)",
            "helper(self.width + 1)",
            true,
        ),
        (
            "train.rs",
            "fn drop(&mut self) {}",
            "fn drop(&mut self) { self.width = 0; }",
            true,
        ),
        ("Cargo.lock", "0.2.19", "0.2.20", true),
        (
            "Cargo.toml",
            "ndarray = \"0.17\"",
            "ndarray = \"0.18\"",
            true,
        ),
    ];
    for (file, from, to, moves) in cases {
        let (host, device) = keyed(Some((file, from, to)));
        assert_eq!(host != base.0, moves, "{file}: {from} -> {to}");
        assert_eq!(
            device, base.1,
            "{file}: {from} -> {to} moved the device key"
        );
    }
    let (host, device) = keyed(Some(("device/mod.rs", "2", "3")));
    assert_eq!((host == base.0, device == base.1), (true, false));
}

#[test]
fn a_renamed_entry_stops_the_build() {
    let renamed =
        std::panic::catch_unwind(|| keyed(Some(("train.rs", "pub fn trained", "pub fn learned"))));
    assert!(renamed.is_err());
}
