use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicUsize, Ordering};

static NAMES: OnceLock<Vec<String>> = OnceLock::new();
static CALLS: AtomicUsize = AtomicUsize::new(0);

pub(crate) fn directory() -> Option<PathBuf> {
    std::env::var_os("ROSELLA_JOIN_REPORT").map(PathBuf::from)
}

pub fn remember(names: &[String]) {
    if directory().is_some() {
        let _ = NAMES.set(names.to_vec());
    }
}

pub(crate) struct Report {
    call: usize,
    tests: BufWriter<std::fs::File>,
    members: BufWriter<std::fs::File>,
    written: std::collections::HashSet<usize>,
    pieces: usize,
}

impl Report {
    pub(crate) fn open() -> Option<Self> {
        let directory = directory()?;
        std::fs::create_dir_all(&directory).ok()?;
        let process = std::process::id();
        let open = |name: &str| {
            OpenOptions::new()
                .create(true)
                .append(true)
                .open(directory.join(format!("{name}.{process}.tsv")))
                .map(BufWriter::new)
                .ok()
        };
        Some(Self {
            call: CALLS.fetch_add(1, Ordering::Relaxed),
            tests: open("joins")?,
            members: open("members")?,
            written: Default::default(),
            pieces: 0,
        })
    }

    fn names(&mut self, kind: char, id: usize, contigs: &[usize]) {
        let Some(names) = NAMES.get() else {
            return;
        };
        let joined = contigs
            .iter()
            .map(|contig| names[*contig].as_str())
            .collect::<Vec<_>>()
            .join(",");
        let _ = writeln!(self.members, "{}\t{kind}\t{id}\t{joined}", self.call);
    }

    pub(crate) fn piece(&mut self, contigs: &[usize]) -> usize {
        self.pieces += 1;
        self.names('p', self.pieces, contigs);
        self.pieces
    }

    pub(crate) fn test(&mut self, rival: usize, contigs: &[usize], fields: &[String]) {
        if self.written.insert(rival) {
            self.names('c', rival, contigs);
        }
        let _ = writeln!(
            self.tests,
            "{}\t{}\t{rival}\t{}",
            self.call,
            self.pieces,
            fields.join("\t")
        );
    }
}
