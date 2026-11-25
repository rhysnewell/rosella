use std::collections::HashMap;
use std::io::{BufRead, BufWriter, Write};
use std::sync::LazyLock;

use anyhow::{Context, Result, bail};
use flate2::{Compression, write::GzEncoder};
use ndarray::{Array2, ArrayView1, ArrayView2};

pub static VECTORS: LazyLock<Option<String>> =
    LazyLock::new(|| std::env::var("ROSELLA_VECTORS").ok());

// Given vectors already carry depth, so by default they are read alone, as their own binner reads
// them.
pub static WEIGHT: LazyLock<Option<f64>> =
    LazyLock::new(|| match std::env::var("ROSELLA_VECTORS_WEIGHT") {
        Ok(value) => Some(
            value
                .parse()
                .expect("ROSELLA_VECTORS_WEIGHT is not a number"),
        ),
        Err(_) => VECTORS.as_ref().map(|_| 0.0),
    });

pub static DUMP: LazyLock<Option<String>> =
    LazyLock::new(|| std::env::var("ROSELLA_SPACE_DUMP").ok());

pub fn dump(path: &str, names: &[String], rows: ArrayView2<f64>) -> Result<()> {
    let file = std::fs::File::create(path).with_context(|| format!("writing {path}"))?;
    let mut out = BufWriter::new(GzEncoder::new(file, Compression::default()));
    let columns = (0..rows.ncols()).map(|column| column.to_string());
    writeln!(out, "contig\t{}", columns.collect::<Vec<_>>().join("\t"))?;
    for (name, row) in names.iter().zip(rows.rows()) {
        write!(out, "{name}")?;
        for value in row {
            write!(out, "\t{value}")?;
        }
        writeln!(out)?;
    }
    out.into_inner()
        .map_err(|error| error.into_error())?
        .finish()?;
    Ok(())
}

pub struct Given {
    rows: HashMap<String, Vec<f64>>,
    width: usize,
}

impl Given {
    pub fn read(path: &str) -> Result<Self> {
        let reader =
            crate::get_file_reader(path).with_context(|| format!("reading vectors {path}"))?;
        let mut rows = HashMap::new();
        for line in reader.lines().skip(1) {
            let line = line?;
            let mut fields = line.split('\t');
            let name = fields.next().unwrap_or_default().to_string();
            let values = fields
                .map(str::parse)
                .collect::<Result<Vec<f64>, _>>()
                .with_context(|| format!("reading the vector of {name} in {path}"))?;
            rows.insert(name, values);
        }
        let width = rows.values().next().map_or(0, Vec::len);
        if let Some((name, _)) = rows.iter().find(|(_, values)| values.len() != width) {
            bail!("{name} in {path} is not {width} wide");
        }
        Ok(Self { rows, width })
    }

    pub fn rows(&self, names: &[String]) -> Result<Array2<f64>> {
        let mut table = Array2::zeros((names.len(), self.width));
        for (mut row, name) in table.rows_mut().into_iter().zip(names) {
            let Some(values) = self.rows.get(name) else {
                bail!("{name} has no given vector");
            };
            row.assign(&ArrayView1::from(values));
        }
        Ok(table)
    }
}
