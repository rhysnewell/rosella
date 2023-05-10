use std::{path::Path, collections::HashSet};

use anyhow::{Result, anyhow};
use ndarray::{Array, Array2, prelude::*};

use crate::external::coverm_engine::MappingMode;


pub struct CoverageTable {
    pub table: Array2<f64>, // rows are contigs, columns are coverage and variance. number of columns is twice the number of samples.
                            // the first column is the coverage, the second is the variance, the third is the coverage, the fourth is the variance, etc.
    pub average_depths: Vec<f64>, // the average of the coverage values in each row for a contig. Order is identical to the order of the rows of table.
    pub contig_names: Vec<String>, // same length as the rows of table. Order is identical to the order of the rows of table.
    pub contig_lengths: Vec<usize>, // same length as the rows of table. Order is identical to the order of the rows of table.
    pub sample_names: Vec<String>, // half the length of the columns of table. Order is identical to the order of the columns of table.
}

impl CoverageTable {
    pub fn new(
        table: Array2<f64>,
        average_depths: Vec<f64>,
        contig_names: Vec<String>,
        contig_lengths: Vec<usize>,
        sample_names: Vec<String>,
    ) -> Self {
        Self {
            table,
            average_depths,
            contig_names,
            contig_lengths,
            sample_names,
        }
    }

    pub fn filter(&mut self, min_contig_size: usize) -> Result<()> {
        // find the indices of the contigs that are too small
        let indices_to_remove = self.contig_lengths
            .iter()
            .enumerate()
            .filter_map(|(index, length)| {
                if *length < min_contig_size {
                    Some(index)
                } else {
                    None
                }
            }).collect::<HashSet<_>>();
        
        // remove the contigs from the table
        let new_table = self.table
            .axis_iter(Axis(0))
            .enumerate()
            .filter_map(|(index, row)| {
                if indices_to_remove.contains(&index) {
                    None
                } else {
                    Some(row)
                }
            }).flat_map(|row| row.to_vec());
        let new_n_rows = self.table.nrows() - indices_to_remove.len();
        self.table = Array::from_iter(new_table).into_shape((new_n_rows, self.table.ncols()))?;
        
        // remove the contigs from the average depths
        self.average_depths = self.average_depths
            .iter()
            .enumerate()
            .filter_map(|(index, depth)| {
                if indices_to_remove.contains(&index) {
                    None
                } else {
                    Some(*depth)
                }
            }).collect::<Vec<_>>();
        
        // remove the contigs from the contig names
        self.contig_names = self.contig_names
            .iter()
            .enumerate()
            .filter_map(|(index, name)| {
                if indices_to_remove.contains(&index) {
                    None
                } else {
                    Some(name.clone())
                }
            }).collect::<Vec<_>>();
        
        // remove the contigs from the contig lengths
        self.contig_lengths = self.contig_lengths
            .iter()
            .enumerate()
            .filter_map(|(index, length)| {
                if indices_to_remove.contains(&index) {
                    None
                } else {
                    Some(*length)
                }
            }).collect::<Vec<_>>();

        Ok(())
    }

    /// read a coverage table from a file
    /// we specify the mode as a parameter as coverm has
    /// different output formats for different modes depending on 
    /// short/long read inputs
    pub fn from_file<P: AsRef<Path>>(file_path: P, mode: MappingMode) -> Result<Self> {
        match mode {
            MappingMode::ShortBam | MappingMode::ShortRead => {
                let mut reader = csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .has_headers(true)
                    .from_path(file_path)?;

                let mut table = Vec::new();
                let mut contig_names = Vec::new();
                let mut contig_lengths = Vec::new();
                let mut average_depths = Vec::new();
                let mut sample_names = Vec::new();

                // get sample name from header
                let headers = reader.headers()?;
                for header in headers.iter().skip(3).step_by(2) {
                    if header.contains("/") {
                        // when performing read mapping, coverm sets the column name to
                        // {reference}/{sample}.bam
                        let mut sample_name = header.split("/").last().unwrap().to_string();
                        sample_name = sample_name.replace(".bam", "");
                        sample_names.push(sample_name);
                    } else {
                        // when using BAM files, coverm sets the column name to
                        // {sample}
                        sample_names.push(header.to_string());
                    }
                }

                for result in reader.records() {
                    let record = result?;
                    let mut record_iter = record.iter();

                    let contig_name = record_iter.next().unwrap().to_string();
                    let contig_length = record_iter.next().unwrap().parse::<usize>()?;
                    let average_depth = record_iter.next().unwrap().parse::<f64>()?;

                    let mut coverage = Vec::new();
                    let mut variance = Vec::new();

                    for (i, value) in record_iter.enumerate() {
                        if i % 2 == 0 {
                            coverage.push(value.parse::<f64>()?);
                        } else {
                            variance.push(value.parse::<f64>()?);
                        }
                    }

                    table.push(coverage);
                    table.push(variance);

                    contig_names.push(contig_name);
                    contig_lengths.push(contig_length);
                    average_depths.push(average_depth);
                }

                let table = Array2::from_shape_vec(
                    (contig_names.len(), sample_names.len() * 2),
                    table.into_iter().flatten().collect(),
                )?;

                Ok(
                    Self {
                        table,
                        average_depths,
                        contig_names,
                        contig_lengths,
                        sample_names,
                    }
                )
            },
            MappingMode::LongBam | MappingMode::LongRead => {
                // long read/bam output is different.
                // the first column is still the contig name
                // the second column is the contig length
                // the third column is the sample coverage
                // the fourth column is the sample variance
                // but then the fifth column is the contig length again, just
                // reclalculated for the second sample. So we need to ignore every extra
                // length column
                let mut reader = csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .has_headers(true)
                    .from_path(file_path)?;

                let mut table = Vec::new();
                let mut contig_names = Vec::new();
                let mut contig_lengths = Vec::new();
                // we need to calculate average depths ourselves after collecting the table
                let mut average_depths = Vec::new();
                let mut sample_names = Vec::new();

                // get sample name from header
                let headers = reader.headers()?;
                // skip 2 and then step by 3 to bypase length columns
                for header in headers.iter().skip(2).step_by(3) {
                    // split on white space to get rid of "Mean" or "Variance" in header name
                    let mut sample_name = header.split_whitespace().next().unwrap().to_string();
                    if header.contains("/") {
                        // when performing read mapping, coverm sets the column name to
                        // {reference}/{sample}.bam
                        sample_name = sample_name.split("/").last().unwrap().to_string();
                        sample_name = sample_name.replace(".bam", "");
                        sample_names.push(sample_name);
                    } else {
                        // when using BAM files, coverm sets the column name to
                        // {sample}
                        sample_names.push(sample_name.to_string());
                    }
                }

                for result in reader.records() {
                    let record = result?;
                    let mut record_iter = record.iter();

                    let contig_name = record_iter.next().unwrap().to_string();
                    let contig_length = record_iter.next().unwrap().parse::<usize>()?;

                    let mut coverage = Vec::new();
                    let mut variance = Vec::new();

                    for (i, value) in record_iter.enumerate() {
                        if i % 2 == 0 {
                            coverage.push(value.parse::<f64>()?);
                        } else {
                            variance.push(value.parse::<f64>()?);
                        }
                    }

                    let average_depth = coverage.iter().sum::<f64>() / coverage.len() as f64;

                    table.push(coverage);
                    table.push(variance);

                    contig_names.push(contig_name);
                    contig_lengths.push(contig_length);
                    average_depths.push(average_depth);
                }

                let table = Array2::from_shape_vec(
                    (contig_names.len(), sample_names.len() * 2),
                    table.into_iter().flatten().collect(),
                )?;

                Ok(
                    Self {
                        table,
                        average_depths,
                        contig_names,
                        contig_lengths,
                        sample_names,
                    }
                )
            }
        }
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
                },
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
    pub fn write<P: AsRef<Path>>(&self, output_path: P) -> Result<()> {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(output_path)?;

        // write header row
        writer.write_field("contigName")?;
        writer.write_field("contigLength")?;
        writer.write_field("totalAvgDepth")?;
        for sample_name in &self.sample_names {
            writer.write_field(format!("{}", sample_name))?;
            writer.write_field(format!("{}-var", sample_name))?;
        }
        writer.write_record(None::<&[u8]>)?;

        // write table
        for (((contig_name, contig_length), average_depth), row) in 
            self.contig_names.iter().zip(
                self.contig_lengths.iter()
            ).zip(
                self.average_depths.iter()
            ).zip(
                self.table.axis_iter(Axis(0)))
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
}