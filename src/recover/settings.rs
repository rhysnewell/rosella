use crate::{embedding::metrics::DistanceSettings, seeds::Seeds};

pub fn seeds(params: &crate::cli::SeedParams) -> Seeds {
    Seeds {
        seed: params.seed,
        knn: params.knn.unwrap_or(crate::defaults::KNN_SEED),
        partition: params.partition.unwrap_or(crate::defaults::PARTITION_SEED),
    }
}

pub fn distance_settings() -> DistanceSettings {
    DistanceSettings {
        presence_fraction: crate::tuning::PRESENCE_FRACTION,
        aggregate_weight: None,
    }
}

pub const DISSOLVE_NAMES: [&str; 2] = ["on", "off"];

pub fn dissolve(choice: &str) -> bool {
    choice != "off"
}
