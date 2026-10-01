use std::{
    collections::{HashMap, HashSet},
    io::BufRead,
    path::Path,
};

use anyhow::{Context, Result, anyhow};
use ndarray::{Array2, Axis};

use crate::external::coverm_engine::MappingMode;
use crate::get_file_reader;

/// One row per contig, ordered the same way in every field. A row is the per sample mean
/// and variance interleaved, which is the order the coverage distance reads it in.
pub struct CoverageTable {
    pub table: Array2<f64>,
    pub average_depths: Vec<f64>,
    pub contig_names: Vec<String>,
    pub contig_lengths: Vec<usize>,
    pub sample_names: Vec<String>,
}

impl CoverageTable {
    pub fn filter_by_length(&mut self, min_contig_size: usize) -> Result<HashSet<String>> {
        let indices_to_remove = self
            .contig_lengths
            .iter()
            .enumerate()
            .filter_map(|(index, length)| {
                if *length < min_contig_size {
                    Some(index)
                } else {
                    None
                }
            })
            .collect::<HashSet<_>>();

        self.filter_by_index(&indices_to_remove)
    }

    pub fn filter_by_index(
        &mut self,
        indices_to_remove: &HashSet<usize>,
    ) -> Result<HashSet<String>> {
        let removed = crate::rows::dropped_names(&self.contig_names, indices_to_remove);
        self.table = crate::rows::keep_rows(&self.table, indices_to_remove)?;
        self.average_depths = crate::rows::keep(&self.average_depths, indices_to_remove);
        self.contig_names = crate::rows::keep(&self.contig_names, indices_to_remove);
        self.contig_lengths = crate::rows::keep(&self.contig_lengths, indices_to_remove);
        Ok(removed)
    }

    /// Read a coverage table, taking the column layout from the run mode. CoverM emits a
    /// different table for short and long reads.
    pub fn from_file<P: AsRef<Path>>(file_path: P, mode: MappingMode) -> Result<Self> {
        Self::read(file_path, Some(Layout::of(mode)))
    }

    /// Read a table whose layout is taken from its own header. A `--coverage-file` rosella
    /// did not write itself can be either layout, and the run mode does not know which.
    pub fn from_any_file<P: AsRef<Path>>(file_path: P) -> Result<Self> {
        Self::read(file_path, None)
    }

    fn reader(file_path: &Path) -> Result<csv::Reader<Box<dyn BufRead>>> {
        let source = get_file_reader(file_path)
            .with_context(|| format!("reading the coverage table {}", file_path.display()))?;
        Ok(csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_reader(source))
    }

    /// The samples a table already holds, without reading its rows.
    pub fn sample_names_in<P: AsRef<Path>>(file_path: P) -> Result<Vec<String>> {
        let mut reader = Self::reader(file_path.as_ref())?;
        let headers = reader.headers()?.clone();
        Ok(Layout::detect(&headers, file_path.as_ref())?.sample_names(&headers))
    }

    fn read<P: AsRef<Path>>(file_path: P, layout: Option<Layout>) -> Result<Self> {
        let mut table = Vec::new();
        let mut contig_names = Vec::new();
        let mut contig_lengths = Vec::new();
        let mut average_depths = Vec::new();
        let sample_names = Self::visit(
            file_path.as_ref(),
            layout,
            |_| true,
            |row| {
                table.push(row.values);
                contig_names.push(row.name);
                contig_lengths.push(row.length);
                average_depths.push(row.average_depth);
            },
        )?;

        let table = Array2::from_shape_vec(
            (contig_names.len(), sample_names.len() * 2),
            table.into_iter().flatten().collect(),
        )?;

        Ok(Self {
            table,
            average_depths,
            contig_names,
            contig_lengths,
            sample_names,
        })
    }

    pub fn rows_named<P: AsRef<Path>>(file_path: P, names: &[&str]) -> Result<Array2<f64>> {
        let wanted = names
            .iter()
            .enumerate()
            .map(|(at, name)| (*name, at))
            .collect::<HashMap<_, _>>();
        let mut rows = vec![Vec::new(); names.len()];
        Self::visit(
            file_path.as_ref(),
            None,
            |name| wanted.contains_key(name),
            |row| rows[wanted[row.name.as_str()]] = row.values,
        )?;
        if let Some(at) = rows.iter().position(Vec::is_empty) {
            bail!("{} is not in {}", names[at], file_path.as_ref().display());
        }
        let width = rows.first().map_or(0, Vec::len);
        Ok(Array2::from_shape_vec((names.len(), width), rows.concat())?)
    }

    // A band reads a few rows of a table that can hold millions, so a row is named before its
    // numbers are parsed.
    fn visit(
        file_path: &Path,
        layout: Option<Layout>,
        keep: impl Fn(&str) -> bool,
        mut visit: impl FnMut(Row),
    ) -> Result<Vec<String>> {
        let mut reader = Self::reader(file_path)?;
        let headers = reader.headers()?.clone();
        let layout = match layout {
            Some(layout) => layout,
            None => Layout::detect(&headers, file_path)?,
        };
        let sample_names = layout.sample_names(&headers);
        for result in reader.records() {
            let record = result?;
            if record.get(0).is_some_and(|name| !keep(name)) {
                continue;
            }
            let row = layout.parse(&record)?;
            if row.values.len() != sample_names.len() * 2 {
                bail!(
                    "{} has {} samples in its header but {} mean and variance columns on \
                     contig {}",
                    file_path.display(),
                    sample_names.len(),
                    row.values.len(),
                    row.name
                );
            }
            visit(row);
        }
        Ok(sample_names)
    }

    pub fn merge(&mut self, other: Self) -> Result<()> {
        if self.contig_names != other.contig_names {
            return Err(anyhow!(
                "Cannot merge coverage tables with different contig names",
            ));
        }
        if let Some(name) = other
            .sample_names
            .iter()
            .find(|name| self.sample_names.contains(name))
        {
            bail!("two coverage tables both hold a sample named {name}");
        }

        self.sample_names.extend(other.sample_names);
        self.table.append(Axis(1), other.table.view())?;

        self.average_depths = self
            .table
            .axis_iter(Axis(0))
            .map(|row| row.iter().step_by(2).sum::<f64>() / self.sample_names.len() as f64)
            .collect();

        Ok(())
    }

    /// A resumed directory computes only the samples it is missing and appends them, so the
    /// columns come out in a different order from a fresh run and the folds over them differ.
    pub fn align_to(&mut self, wanted: &[&str]) {
        let at = self
            .sample_names
            .iter()
            .enumerate()
            .map(|(index, name)| (name.as_str(), index))
            .collect::<HashMap<_, _>>();
        let mut order = wanted
            .iter()
            .filter_map(|name| at.get(name).copied())
            .collect::<Vec<_>>();
        let asked = order.iter().copied().collect::<HashSet<_>>();
        order.extend((0..self.sample_names.len()).filter(|index| !asked.contains(index)));
        if order.iter().copied().eq(0..self.sample_names.len()) {
            return;
        }

        let mut table = Array2::zeros((self.table.nrows(), order.len() * 2));
        for (slot, source) in order.iter().enumerate() {
            table
                .column_mut(slot * 2)
                .assign(&self.table.column(source * 2));
            table
                .column_mut(slot * 2 + 1)
                .assign(&self.table.column(source * 2 + 1));
        }
        self.table = table;
        self.sample_names = order
            .iter()
            .map(|index| self.sample_names[*index].clone())
            .collect();
    }

    pub fn merge_many(coverage_tables: Vec<CoverageTable>) -> Result<Self> {
        let mut tables = coverage_tables.into_iter();
        let mut merged = tables
            .next()
            .ok_or_else(|| anyhow!("No coverage file or reads provided."))?;
        for table in tables {
            merged.merge(table)?;
        }
        Ok(merged)
    }

    /// Written in CoverM's contigName, contigLen, totalAvgDepth layout with each sample's mean
    /// and `-var` column after, so a later run reads it back through the same parser.
    pub fn write<P: AsRef<Path>>(&self, output_path: P) -> Result<()> {
        crate::report_sink::write_atomically(output_path.as_ref(), |file| self.rows_into(file))
    }

    fn rows_into(&self, file: &std::fs::File) -> Result<()> {
        let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_writer(file);

        writer.write_field("contigName")?;
        writer.write_field("contigLen")?;
        writer.write_field("totalAvgDepth")?;
        for sample_name in &self.sample_names {
            writer.write_field(sample_name)?;
            writer.write_field(format!("{}-var", sample_name))?;
        }
        writer.write_record(None::<&[u8]>)?;

        for (((contig_name, contig_length), average_depth), row) in self
            .contig_names
            .iter()
            .zip(self.contig_lengths.iter())
            .zip(self.average_depths.iter())
            .zip(self.table.axis_iter(Axis(0)))
        {
            let mut record = Vec::with_capacity(3 + self.sample_names.len() * 2);
            record.push(contig_name.to_string());
            record.push(format!("{}", contig_length));
            record.push(average_depth.to_string());
            record.extend(row.iter().map(f64::to_string));
            writer.write_record(record)?;
        }
        writer.flush()?;

        Ok(())
    }
}

/// CoverM writes a different table per `--methods` choice. The depth table carries one
/// length column for the contig; the long read triple repeats the length once per sample,
/// which is why the two cannot share a stride.
#[derive(Clone, Copy)]
enum Layout {
    SharedLength,
    PerSampleLength,
}

struct Row {
    name: String,
    length: usize,
    average_depth: f64,
    /// Per sample mean and variance, interleaved, which is the order the coverage distance
    /// reads a row in.
    values: Vec<f64>,
}

impl Layout {
    fn of(mode: MappingMode) -> Self {
        match mode {
            MappingMode::ShortBam | MappingMode::ShortRead => Self::SharedLength,
            MappingMode::LongBam | MappingMode::LongRead => Self::PerSampleLength,
        }
    }

    fn detect(headers: &csv::StringRecord, path: &Path) -> Result<Self> {
        match headers.get(0) {
            Some("contigName") => Ok(Self::SharedLength),
            Some("Contig") => Ok(Self::PerSampleLength),
            _ => bail!(
                "{} starts with neither the contigName column of CoverM's depth table nor the \
                 Contig column its length, trimmed_mean and variance methods write",
                path.display()
            ),
        }
    }

    fn sample_names(&self, headers: &csv::StringRecord) -> Vec<String> {
        match self {
            Self::SharedLength => headers
                .iter()
                .skip(3)
                .step_by(2)
                .map(|header| bam_stem(header).to_string())
                .collect(),
            Self::PerSampleLength => headers
                .iter()
                .skip(2)
                .step_by(3)
                .map(|header| bam_stem(header.split_whitespace().next().unwrap_or(header)))
                .map(str::to_string)
                .collect(),
        }
    }

    fn parse(&self, record: &csv::StringRecord) -> Result<Row> {
        let mut fields = record.iter();
        let name = fields
            .next()
            .ok_or_else(|| anyhow!("a coverage row has no contig name"))?
            .to_string();

        match self {
            Self::SharedLength => {
                let length = number(fields.next(), &name)? as usize;
                let average_depth = number(fields.next(), &name)?;
                let values = fields
                    .map(|field| number(Some(field), &name))
                    .collect::<Result<Vec<_>>>()?;
                Ok(Row {
                    name,
                    length,
                    average_depth,
                    values,
                })
            }
            Self::PerSampleLength => {
                let columns = fields.collect::<Vec<_>>();
                let mut length = 0;
                let mut values = Vec::with_capacity(columns.len() / 3 * 2);
                for (sample, triple) in columns.chunks(3).enumerate() {
                    if triple.len() < 3 {
                        bail!(
                            "contig {} has a trailing partial length, trimmed mean and \
                             variance triple",
                            name
                        );
                    }
                    if sample == 0 {
                        length = number(Some(triple[0]), &name)? as usize;
                    }
                    values.push(number(Some(triple[1]), &name)?);
                    values.push(number(Some(triple[2]), &name)?);
                }
                if values.is_empty() {
                    bail!("contig {} has no coverage columns", name);
                }
                let means = values.iter().step_by(2);
                let average_depth = means.clone().sum::<f64>() / means.count() as f64;
                Ok(Row {
                    name,
                    length,
                    average_depth,
                    values,
                })
            }
        }
    }
}

/// CoverM names a mapped column `{reference}/{sample}.bam` and a BAM column `{sample}`, so
/// a read or BAM path names its sample the way a column does only once stemmed the same way.
pub(crate) fn bam_stem(path: &str) -> &str {
    let name = Path::new(path)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or(path);
    name.strip_suffix(".bam").unwrap_or(name)
}

fn number(field: Option<&str>, contig: &str) -> Result<f64> {
    let field = field.ok_or_else(|| anyhow!("contig {} has too few columns", contig))?;
    field
        .trim()
        .parse::<f64>()
        .map_err(|_| anyhow!("contig {} has `{}` where a number belongs", contig, field))
}
