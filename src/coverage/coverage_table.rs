use std::path::Path;

use anyhow::Result;
use ndarray::Array2;

use crate::external::coverm_engine::MappingMode;


pub struct CoverageTable {
    pub table: Array2<f32>, // rows are contigs, columns are coverage and variance. number of columns is twice the number of samples.
                            // the first column is the coverage, the second is the variance, the third is the coverage, the fourth is the variance, etc.
    pub average_depths: Vec<f32>, // the average of the coverage values in each row for a contig. Order is identical to the order of the rows of table.
    pub contig_names: Vec<String>, // same length as the rows of table. Order is identical to the order of the rows of table.
    pub contig_lengths: Vec<usize>, // same length as the rows of table. Order is identical to the order of the rows of table.
    pub sample_names: Vec<String>, // half the length of the columns of table. Order is identical to the order of the columns of table.
}

impl CoverageTable {
    pub fn new(
        table: Array2<f32>,
        average_depths: Vec<f32>,
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
                for (i, header) in headers.iter().skip(3).step_by(2).enumerate() {
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
                    let average_depth = record_iter.next().unwrap().parse::<f32>()?;

                    let mut coverage = Vec::new();
                    let mut variance = Vec::new();

                    for (i, value) in record_iter.enumerate() {
                        if i % 2 == 0 {
                            coverage.push(value.parse::<f32>()?);
                        } else {
                            variance.push(value.parse::<f32>()?);
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
                for (i, header) in headers.iter().skip(2).step_by(3).enumerate() {
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
                            coverage.push(value.parse::<f32>()?);
                        } else {
                            variance.push(value.parse::<f32>()?);
                        }
                    }

                    let average_depth = coverage.iter().sum::<f32>() / coverage.len() as f32;

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

    pub fn merge(coverage_tables: Vec<CoverageTable>) -> Self {
        // TODO: Rewrite this
        
    }
}