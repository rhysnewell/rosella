use anyhow::{Result, bail};

pub struct Tree {
    split_feature: Vec<u32>,
    threshold: Vec<f64>,
    left: Vec<i32>,
    right: Vec<i32>,
    leaf_value: Vec<f64>,
}

/// The trained objective is `regression sqrt`, so the boosted sum is the square root of the
/// answer and the sign has to survive the squaring.
pub struct Booster {
    trees: Vec<Tree>,
}

fn numbers<T: std::str::FromStr>(line: &str) -> Result<Vec<T>> {
    line.split_whitespace()
        .map(|field| {
            field
                .parse::<T>()
                .map_err(|_| anyhow::anyhow!("{field} is not a number"))
        })
        .collect()
}

impl Booster {
    pub fn parse(model: &str) -> Result<Self> {
        let mut trees = Vec::new();
        let mut fields: std::collections::HashMap<&str, &str> = std::collections::HashMap::new();
        let mut open = false;

        let close = |fields: &mut std::collections::HashMap<&str, &str>,
                     trees: &mut Vec<Tree>|
         -> Result<()> {
            let leaf_value = numbers(fields.get("leaf_value").copied().unwrap_or_default())?;
            let tree = match fields.get("split_feature") {
                Some(split) => Tree {
                    split_feature: numbers(split)?,
                    threshold: numbers(fields["threshold"])?,
                    left: numbers(fields["left_child"])?,
                    right: numbers(fields["right_child"])?,
                    leaf_value,
                },
                None => Tree {
                    split_feature: Vec::new(),
                    threshold: Vec::new(),
                    left: Vec::new(),
                    right: Vec::new(),
                    leaf_value,
                },
            };
            fields.clear();
            trees.push(tree);
            Ok(())
        };

        for line in model.lines() {
            let line = line.trim();
            if let Some(rest) = line.strip_prefix("Tree=") {
                let _ = rest;
                if open {
                    close(&mut fields, &mut trees)?;
                }
                open = true;
            } else if line == "end of trees" {
                if open {
                    close(&mut fields, &mut trees)?;
                    open = false;
                }
            } else if open {
                if let Some((key, value)) = line.split_once('=') {
                    fields.insert(key, value);
                }
            }
        }
        if open {
            close(&mut fields, &mut trees)?;
        }
        if trees.is_empty() {
            bail!("the model holds no trees");
        }
        Ok(Self { trees })
    }

    pub fn predict(&self, features: &[f64]) -> f64 {
        let raw = self
            .trees
            .iter()
            .map(|tree| tree.walk(features))
            .sum::<f64>();
        raw.signum() * raw * raw
    }
}

impl Tree {
    fn walk(&self, features: &[f64]) -> f64 {
        if self.split_feature.is_empty() {
            return self.leaf_value.first().copied().unwrap_or(0.0);
        }
        let mut node = 0usize;
        loop {
            let feature = self.split_feature[node] as usize;
            let value = features.get(feature).copied().unwrap_or(0.0);
            let next = match value <= self.threshold[node] {
                true => self.left[node],
                false => self.right[node],
            };
            if next < 0 {
                return self.leaf_value[!next as usize];
            }
            node = next as usize;
        }
    }
}
