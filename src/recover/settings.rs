use anyhow::Result;

use crate::{
    embedding::metrics::DistanceSettings,
    seeds::Seeds,
};

pub fn seeds(params: &crate::cli::SeedParams) -> Seeds {
    Seeds {
        knn: params.knn.unwrap_or(params.seed),
        sample: params.sample.unwrap_or(params.seed),
        partition: params.partition.unwrap_or(params.seed),
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
