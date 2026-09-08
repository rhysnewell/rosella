use std::collections::HashMap;
use std::io::Read;

use anyhow::{Result, bail};

const FEATURES_GZ: &[u8] = include_bytes!("../../data/checkm2_features.tsv.gz");

pub const METADATA: usize = 22;

/// A derived column is the share of its group's genes the bin holds. Pathways and categories
/// count a gene once however often it is seen; modules sum the raw counts and divide by the
/// definition length before any gene is dropped, which is the trained model's arithmetic.
pub struct Group {
    pub length: usize,
    pub columns: Vec<u32>,
}

pub struct Tables {
    pub kos: HashMap<String, u32>,
    pub gene_count: usize,
    pub pathways: Vec<Group>,
    pub modules: Vec<Group>,
    pub categories: Vec<Group>,
}

fn groups(lines: &mut std::str::Lines, count: usize) -> Result<Vec<Group>> {
    (0..count)
        .map(|_| {
            let Some((length, columns)) = lines.next().and_then(|line| line.split_once('\t'))
            else {
                bail!("the feature table ended inside a group");
            };
            Ok(Group {
                length: length.parse()?,
                columns: columns
                    .split(',')
                    .filter(|field| !field.is_empty())
                    .map(|field| field.parse::<u32>())
                    .collect::<Result<_, _>>()?,
            })
        })
        .collect()
}

impl Tables {
    pub fn load() -> Result<Self> {
        let mut text = String::new();
        flate2::read::GzDecoder::new(FEATURES_GZ).read_to_string(&mut text)?;
        let mut lines = text.lines();

        let mut kos = HashMap::new();
        let mut gene_count = 0;
        let (mut pathways, mut modules, mut categories) = (Vec::new(), Vec::new(), Vec::new());

        while let Some(header) = lines.next() {
            let Some((tag, count)) = header.strip_prefix('#').and_then(|h| h.split_once('\t'))
            else {
                continue;
            };
            let count = count.parse::<usize>()?;
            match tag {
                "metadata" => {
                    lines.next();
                    if count != METADATA {
                        bail!("the model takes {METADATA} metadata columns, the table has {count}");
                    }
                }
                "kos" => {
                    gene_count = count;
                    for index in 0..count {
                        let Some(ko) = lines.next() else {
                            bail!("the feature table ended inside the gene list");
                        };
                        kos.insert(ko.to_string(), index as u32);
                    }
                }
                "pathways" => pathways = groups(&mut lines, count)?,
                "modules" => modules = groups(&mut lines, count)?,
                "categories" => categories = groups(&mut lines, count)?,
                _ => bail!("{tag} is not a feature group"),
            }
        }
        if gene_count == 0 {
            bail!("the feature table holds no genes");
        }
        Ok(Self {
            kos,
            gene_count,
            pathways,
            modules,
            categories,
        })
    }

    pub fn width(&self) -> usize {
        METADATA
            + self.gene_count
            + self.pathways.len()
            + self.modules.len()
            + self.categories.len()
    }

    /// Columns follow the order the models were trained on: metadata, genes, pathways,
    /// modules, categories.
    pub fn fill(&self, metadata: &[f64; METADATA], genes: &[f64], into: &mut Vec<f64>) {
        into.clear();
        into.extend_from_slice(metadata);
        into.extend_from_slice(genes);
        for group in &self.pathways {
            into.push(share(genes, group, true));
        }
        for group in &self.modules {
            into.push(share(genes, group, false));
        }
        for group in &self.categories {
            into.push(share(genes, group, true));
        }
    }
}

fn share(genes: &[f64], group: &Group, presence: bool) -> f64 {
    if group.length == 0 {
        return 0.0;
    }
    let total = group
        .columns
        .iter()
        .map(|column| match genes[*column as usize] {
            count if presence && count > 1.0 => 1.0,
            count => count,
        })
        .sum::<f64>();
    total / group.length as f64
}
