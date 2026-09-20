use crate::{embedding::metrics::DistanceSettings, seeds::Seeds};

pub fn seeds(params: &crate::cli::SeedParams) -> Seeds {
    Seeds {
        seed: params.seed,
        knn: params.knn.unwrap_or(crate::defaults::KNN_SEED),
        partition: params.partition.unwrap_or(crate::defaults::PARTITION_SEED),
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

pub const HOLD_NAMES: [&str; 3] = ["bars", "size", "tier"];

pub fn hold(choice: &str) -> crate::refine::dissolve::Hold {
    use crate::refine::dissolve::Hold;
    match choice {
        "size" => Hold::Size,
        "tier" => Hold::Tier,
        _ => Hold::Bars,
    }
}

pub const RUNG_WALK_NAMES: [&str; 2] = ["walk", "break"];

pub fn rung_walk(choice: &str) -> crate::refine::dissolve::RungWalk {
    use crate::refine::dissolve::RungWalk;
    match choice {
        "break" => RungWalk::Break,
        _ => RungWalk::Walk,
    }
}

pub const DUPLICATE_NAMES: [&str; 2] = ["hits", "carriers"];

pub fn duplicates(choice: &str) -> crate::markers::Duplicates {
    use crate::markers::Duplicates;
    match choice {
        "carriers" => Duplicates::Carriers,
        _ => Duplicates::Hits,
    }
}

pub const CONSERVE_NAMES: [&str; 2] = ["off", "on"];

pub fn conserve(choice: &str) -> crate::refine::dissolve::Conserve {
    use crate::refine::dissolve::Conserve;
    match choice {
        "on" => Conserve::On,
        _ => Conserve::Off,
    }
}

pub const NOVELTY_NAMES: [&str; 2] = ["strict", "gain"];

pub fn novelty(choice: &str) -> crate::refine::join::Novelty {
    use crate::refine::join::Novelty;
    match choice {
        "gain" => Novelty::Gain,
        _ => Novelty::Strict,
    }
}

