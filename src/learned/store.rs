use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use rand::{SeedableRng, rngs::StdRng};

use super::Learned;
use super::corpus::Scales;
use super::towers::Towers;

pub const FORMAT: &str = "rosella-space-1";

pub fn path(directory: &Path, key: &str) -> PathBuf {
    directory.join(format!(
        "space.{:016x}.bin",
        crate::digest::fold(key.as_bytes())
    ))
}

fn header(key: &str) -> String {
    format!("{FORMAT}\t{key}\n")
}

impl Learned {
    pub fn save(&self, path: &Path, key: &str) -> Result<()> {
        fs::create_dir_all(path.parent().unwrap_or(Path::new(".")))?;
        crate::report_sink::write_atomically(path, |file| {
            let mut sink = BufWriter::new(file);
            sink.write_all(header(key).as_bytes())?;
            let split = self.scales.split();
            for width in [split, self.scales.width() - split] {
                sink.write_all(&(width as u64).to_le_bytes())?;
            }
            for part in self.state() {
                for value in part {
                    sink.write_all(&value.to_le_bytes())?;
                }
            }
            sink.flush()?;
            Ok(())
        })
    }

    // The key is read back from the file as well as its name, so a name collision never loads
    // another run's space.
    pub fn load(path: &Path, key: &str) -> Result<Option<Self>> {
        let Ok(bytes) = fs::read(path) else {
            return Ok(None);
        };
        let Some(body) = bytes.strip_prefix(header(key).as_bytes()) else {
            return Ok(None);
        };
        let Some((widths, values)) = body.split_at_checked(16) else {
            bail!("{} is not a whole space", path.display());
        };
        let [composition, depth] = [0, 8].map(|at| {
            u64::from_le_bytes(widths[at..at + 8].try_into().expect("eight bytes")) as usize
        });
        let mut learned = Self::new(
            Towers::new(
                composition + depth,
                composition,
                &mut StdRng::seed_from_u64(0),
            ),
            Scales::blank(composition, depth),
        );
        let (values, tail) = values.as_chunks::<4>();
        let mut state = learned.state_mut();
        if values.len() != state.iter().map(|part| part.len()).sum::<usize>() || !tail.is_empty() {
            bail!("{} is not a whole space", path.display());
        }
        for (value, bytes) in state
            .iter_mut()
            .flat_map(|part| part.iter_mut())
            .zip(values)
        {
            *value = f32::from_le_bytes(*bytes);
        }
        Ok(Some(learned))
    }

    fn state(&self) -> Vec<&[f32]> {
        let mut out = self.scales.parts().to_vec();
        out.extend(self.net.state());
        out
    }

    fn state_mut(&mut self) -> Vec<&mut [f32]> {
        let mut out = Vec::from(self.scales.parts_mut());
        out.extend(self.net.state_mut());
        out
    }
}
