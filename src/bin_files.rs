use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::{File, OpenOptions};
use std::hash::Hash;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use needletail::parser::{LineEnding, write_fasta};

const TOO_MANY_OPEN_FILES: i32 = 24;

// A run can write more bins than the process may hold open, so hitting the limit closes the
// oldest half and a bin that comes back reopens in append mode.
pub struct BinFiles<K, F> {
    path_of: F,
    open: HashMap<K, BufWriter<File>>,
    order: VecDeque<K>,
    seen: HashSet<K>,
}

impl<K: Hash + Eq + Clone, F: Fn(&K) -> PathBuf> BinFiles<K, F> {
    pub fn new(path_of: F) -> Self {
        Self {
            path_of,
            open: HashMap::new(),
            order: VecDeque::new(),
            seen: HashSet::new(),
        }
    }

    pub fn write(&mut self, bin: &K, id: &[u8], sequence: &[u8]) -> Result<()> {
        if !self.open.contains_key(bin) {
            let file = self.open_file(bin)?;
            self.open.insert(bin.clone(), BufWriter::new(file));
            self.order.push_back(bin.clone());
            self.seen.insert(bin.clone());
        }
        let writer = self.open.get_mut(bin).expect("opened above");
        write_fasta(id, sequence, writer, LineEnding::Unix)?;
        Ok(())
    }

    fn open_file(&mut self, bin: &K) -> Result<File> {
        let path = (self.path_of)(bin);
        loop {
            match OpenOptions::new().append(true).create(true).open(&path) {
                Ok(file) => return Ok(file),
                Err(error) if error.raw_os_error() == Some(TOO_MANY_OPEN_FILES) => {
                    if self.open.is_empty() {
                        bail!("no file can be opened for {}: {error}", path.display());
                    }
                    self.close((self.open.len() / 2).max(1))?;
                }
                Err(error) => {
                    return Err(error).with_context(|| format!("opening {}", path.display()));
                }
            }
        }
    }

    fn close(&mut self, count: usize) -> Result<()> {
        for _ in 0..count {
            let Some(bin) = self.order.pop_front() else {
                break;
            };
            if let Some(mut writer) = self.open.remove(&bin) {
                writer
                    .flush()
                    .with_context(|| format!("flushing {}", (self.path_of)(&bin).display()))?;
            }
        }
        Ok(())
    }

    // Dropping a BufWriter flushes it and throws the error away, so a full disk would truncate
    // a bin silently.
    pub fn finish(mut self) -> Result<usize> {
        self.close(self.open.len())?;
        Ok(self.seen.len())
    }
}
