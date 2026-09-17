use crate::{embedding::metrics::DistanceSettings, seeds::Seeds};

pub fn seeds(params: &crate::cli::SeedParams) -> Seeds {
    Seeds {
        seed: params.seed,
        knn: params.knn.unwrap_or(params.seed),
        partition: params.partition.unwrap_or(params.seed),
    }
}

pub fn distance_settings(params: &crate::cli::binning::DistanceParams) -> DistanceSettings {
    DistanceSettings {
        presence_fraction: crate::tuning::PRESENCE_FRACTION,
        aggregate_weight: None,
        calibrate: params.calibrate_composition,
    }
}

pub const DISSOLVE_NAMES: [&str; 2] = ["on", "off"];

pub fn dissolve(choice: &str) -> bool {
    choice != "off"
}

pub const HOLD_NAMES: [&str; 4] = ["bars", "size", "tier", "complete"];

pub fn hold(choice: &str) -> crate::refine::dissolve::Hold {
    use crate::refine::dissolve::Hold;
    match choice {
        "size" => Hold::Size,
        "tier" => Hold::Tier,
        "complete" => Hold::Complete,
        _ => Hold::Bars,
    }
}
