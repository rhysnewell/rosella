use anyhow::Result;

use crate::{
    embedding::metrics::DistanceSettings,
    seeds::Seeds,
};

pub fn seeds(seed: u64, overrides: &crate::cli::SeedOverrides) -> Seeds {
    Seeds {
        knn: overrides.knn.unwrap_or(seed),
        sample: overrides.sample.unwrap_or(seed),
        partition: overrides.partition.unwrap_or(seed),
    }
}

pub fn distance_settings(distance: &crate::cli::DistanceParams) -> Result<DistanceSettings> {
    Ok(DistanceSettings {
        presence_fraction: distance.presence_fraction,
        aggregate_weight: None,
    })
}

pub const DISSOLVE_NAMES: [&str; 2] = ["on", "off"];

pub fn dissolve(choice: &str) -> bool {
    choice != "off"
}
