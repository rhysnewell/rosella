use std::{collections::HashSet, path::Path};

use anyhow::{Result, anyhow};
use ndarray::{Array, Array2, Axis};

use crate::external::coverm_engine::MappingMode;

/// One row per contig, ordered the same way in every field. A row is the per sample mean
/// and variance interleaved, which is the order `metabat` reads it in.
pub struct CoverageTable {
    pub table: Array2<f64>,
    pub average_depths: Vec<f64>,
    pub contig_names: Vec<String>,
    pub contig_lengths: Vec<usize>,
    pub sample_names: Vec<String>,
    pub output_path: String,
}

impl CoverageTable {
    pub fn new(
        table: Array2<f64>,
        average_depths: Vec<f64>,
        contig_names: Vec<String>,
        contig_lengths: Vec<usize>,
        sample_names: Vec<String>,
        output_path: String,
    ) -> Self {
        Self {
            table,
            average_depths,
            contig_names,
            contig_lengths,
            sample_names,
            output_path,
        }
    }

    pub fn filter_by_length(&mut self, min_contig_size: usize) -> Result<HashSet<String>> {
        // find the indices of the contigs that are too small
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

    pub fn clear_variances(&mut self) {
        for mut column in self.table.axis_iter_mut(Axis(1)).skip(1).step_by(2) {
            column.fill(0.0);
        }
    }

    pub fn filter_by_index(
        &mut self,
        indices_to_remove: &HashSet<usize>,
    ) -> Result<HashSet<String>> {
        // remove the contigs from the table
        let new_table = self
            .table
            .axis_iter(Axis(0))
            .enumerate()
            .filter_map(|(index, row)| {
                if indices_to_remove.contains(&index) {
                    None
                } else {
                    Some(row)
                }
            })
            .flat_map(|row| row.to_vec());
        let new_n_rows = self.table.nrows() - indices_to_remove.len();
        self.table =
            Array::from_iter(new_table).into_shape_with_order((new_n_rows, self.table.ncols()))?;

        // remove the contigs from the average depths
        self.average_depths = self
            .average_depths
            .iter()
            .enumerate()
            .filter_map(|(index, depth)| {
                if indices_to_remove.contains(&index) {
                    None
                } else {
                    Some(*depth)
                }
            })
            .collect::<Vec<_>>();

        let filtered_contig_names = self
            .contig_names
            .iter()
            .enumerate()
            .filter_map(|(index, name)| {
                if indices_to_remove.contains(&index) {
                    Some(name.clone())
                } else {
                    None
                }
            })
            .collect::<HashSet<_>>();
        // remove the contigs from the contig names
        self.contig_names = self
            .contig_names
            .iter()
            .enumerate()
            .filter_map(|(index, name)| {
                if indices_to_remove.contains(&index) {
                    None
                } else {
                    Some(name.clone())
                }
            })
            .collect::<Vec<_>>();

        // remove the contigs from the contig lengths
        self.contig_lengths = self
            .contig_lengths
            .iter()
            .enumerate()
            .filter_map(|(index, length)| {
                if indices_to_remove.contains(&index) {
                    None
                } else {
                    Some(*length)
                }
            })
            .collect::<Vec<_>>();

        Ok(filtered_contig_names)
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

    /// The samples a table already holds, without reading its rows.
    pub fn sample_names_in<P: AsRef<Path>>(file_path: P) -> Result<Vec<String>> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(&file_path)?;
        let headers = reader.headers()?.clone();
        Ok(Layout::detect(&headers, file_path.as_ref())?.sample_names(&headers))
    }

    fn read<P: AsRef<Path>>(file_path: P, layout: Option<Layout>) -> Result<Self> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(&file_path)?;

        let headers = reader.headers()?.clone();
        let layout = match layout {
            Some(layout) => layout,
            None => Layout::detect(&headers, file_path.as_ref())?,
        };
        let sample_names = layout.sample_names(&headers);

        let mut table = Vec::new();
        let mut contig_names = Vec::new();
        let mut contig_lengths = Vec::new();
        let mut average_depths = Vec::new();
        for result in reader.records() {
            let row = layout.parse(&result?)?;
            if row.values.len() != sample_names.len() * 2 {
                bail!(
                    "{} has {} samples in its header but {} mean and variance columns on \
                     contig {}",
                    file_path.as_ref().display(),
                    sample_names.len(),
                    row.values.len(),
                    row.name
                );
            }
            table.push(row.values);
            contig_names.push(row.name);
            contig_lengths.push(row.length);
            average_depths.push(row.average_depth);
        }

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
            output_path: file_path.as_ref().to_string_lossy().to_string(),
        })
    }

    pub fn merge(&mut self, other: Self) -> Result<()> {
        if &self.contig_names != &other.contig_names {
            return Err(anyhow!(
                "Cannot merge coverage tables with different contig names",
            ));
        }

        // extend the sample names
        self.sample_names.extend(other.sample_names);
        // merge the tables along the columns
        self.table.append(Axis(1), other.table.view())?;

        // recalculate average depths
        self.average_depths = self
            .table
            .axis_iter(Axis(0))
            .map(|row| row.iter().step_by(2).sum::<f64>() / self.sample_names.len() as f64)
            .collect();

        Ok(())
    }

    /// merge multiple coverage tables into one
    /// Make sure that the tables have the same contig names
    /// in the same order.
    /// We also want to be aware of sample name order.
    /// We will also need to recalculate average depths
    pub fn merge_many(coverage_tables: Vec<CoverageTable>) -> Result<Self> {
        let mut merged_table: Option<CoverageTable> = None;
        for coverage_table in coverage_tables.into_iter() {
            match &mut merged_table {
                Some(table) => {
                    table.merge(coverage_table)?;
                }
                None => {
                    merged_table = Some(coverage_table);
                }
            }
        }

        Ok(merged_table.unwrap())
    }

    /// Write the coverage table to a file
    /// The file will be a tab delimited file with the following columns:
    /// contig_name, contig_length, sample1_coverage, sample1_variance, sample2_coverage, sample2_variance, ...
    /// The first row will be a header row with the sample names
    pub fn write<P: AsRef<Path>>(&mut self, output_path: P) -> Result<()> {
        self.set_output_path(output_path.as_ref().to_string_lossy().to_string());
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(output_path)?;

        // write header row
        writer.write_field("contigName")?;
        writer.write_field("contigLen")?;
        writer.write_field("totalAvgDepth")?;
        for sample_name in &self.sample_names {
            writer.write_field(format!("{}", sample_name))?;
            writer.write_field(format!("{}-var", sample_name))?;
        }
        writer.write_record(None::<&[u8]>)?;

        // write table
        for (((contig_name, contig_length), average_depth), row) in self
            .contig_names
            .iter()
            .zip(self.contig_lengths.iter())
            .zip(self.average_depths.iter())
            .zip(self.table.axis_iter(Axis(0)))
        {
            let mut record = Vec::with_capacity(3 + self.sample_names.len() * 2);
            record.push(format!("{}", contig_name));
            record.push(format!("{}", contig_length));
            record.push(format!("{:.3}", average_depth));
            for value in row {
                record.push(format!("{:.3}", value));
            }
            writer.write_record(record)?;
        }
        writer.flush()?;

        Ok(())
    }

    pub fn set_output_path(&mut self, output_path: String) {
        self.output_path = output_path;
    }
}

/// CoverM writes a different table per `--methods` choice. `metabat` carries one length
/// column for the contig; the long read triple repeats the length once per sample, which
/// is why the two cannot share a stride.
#[derive(Clone, Copy)]
enum Layout {
    Metabat,
    PerSampleLength,
}

struct Row {
    name: String,
    length: usize,
    average_depth: f64,
    /// Per sample mean and variance, interleaved, which is the order `metabat` reads a
    /// coverage row in.
    values: Vec<f64>,
}

impl Layout {
    fn of(mode: MappingMode) -> Self {
        match mode {
            MappingMode::ShortBam | MappingMode::ShortRead => Self::Metabat,
            MappingMode::LongBam | MappingMode::LongRead => Self::PerSampleLength,
        }
    }

    fn detect(headers: &csv::StringRecord, path: &Path) -> Result<Self> {
        match headers.get(0) {
            Some("contigName") => Ok(Self::Metabat),
            Some("Contig") => Ok(Self::PerSampleLength),
            _ => bail!(
                "{} starts with neither the contigName column CoverM's metabat method writes \
                 nor the Contig column its length, trimmed_mean and variance methods write",
                path.display()
            ),
        }
    }

    fn sample_names(&self, headers: &csv::StringRecord) -> Vec<String> {
        match self {
            Self::Metabat => headers.iter().skip(3).step_by(2).map(bam_stem).collect(),
            Self::PerSampleLength => headers
                .iter()
                .skip(2)
                .step_by(3)
                .map(|header| bam_stem(header.split_whitespace().next().unwrap_or(header)))
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
            Self::Metabat => {
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

/// CoverM names a mapped column `{reference}/{sample}.bam` and a BAM column `{sample}`.
fn bam_stem(header: &str) -> String {
    let name = header.rsplit('/').next().unwrap_or(header);
    name.strip_suffix(".bam").unwrap_or(name).to_string()
}

fn number(field: Option<&str>, contig: &str) -> Result<f64> {
    let field = field.ok_or_else(|| anyhow!("contig {} has too few columns", contig))?;
    field
        .trim()
        .parse::<f64>()
        .map_err(|_| anyhow!("contig {} has `{}` where a number belongs", contig, field))
}
