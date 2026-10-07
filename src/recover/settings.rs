use anyhow::Result;

use crate::cli::RecoverArgs;
use crate::cli::binning::Unit;
use crate::cli::rescue::Switch;
use crate::clustering::graph_partition::Partition;
use crate::embedding::metrics::likelihood::CoverageDistance;
use crate::embedding::{features::FeatureSettings, knn::KnnSettings, metrics::DistanceSettings};
use crate::recover::recover_engine::{Stage, parse_order};
use crate::refine::rung::Bars;
use crate::seeds::Seeds;

pub fn seeds(params: &crate::cli::SeedParams) -> Seeds {
    Seeds {
        seed: params.seed,
        knn: params.knn.unwrap_or(crate::defaults::KNN_SEED),
        partition: params.partition.unwrap_or(crate::defaults::PARTITION_SEED),
    }
}

pub fn distance_settings(variance_unit: Option<f64>) -> DistanceSettings {
    DistanceSettings {
        presence_fraction: crate::tuning::PRESENCE_FRACTION,
        aggregate_weight: None,
        variance_unit,
        evidence: None,
    }
}

pub fn knn_settings(params: &crate::cli::GraphParams) -> KnnSettings {
    KnnSettings {
        neighbours: params.n_neighbours,
        candidates: params.candidates(),
    }
}

pub fn feature_settings(params: &crate::cli::GraphParams) -> FeatureSettings {
    FeatureSettings {
        distance: distance_settings(params.given_unit()),
        link_weight: params.assembly_graph_weight as f32,
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Lengths {
    pub cutoff: usize,
    pub attach_floor: Option<usize>,
    pub max_bin: usize,
}

#[derive(Debug, Clone)]
pub struct CycleSettings {
    pub refine_rounds: usize,
    pub partition_seeds: usize,
    pub dissolve: bool,
    pub dissolve_rounds: usize,
    pub dissolve_passes: usize,
    pub join: bool,
    pub order: Vec<Stage>,
    pub split: bool,
    pub peel: bool,
    pub passengers: bool,
    pub replicons: bool,
}

#[derive(Debug, Clone)]
pub struct Settings {
    pub knn: KnnSettings,
    pub seeds: Seeds,
    pub bars: Bars,
    pub lengths: Lengths,
    pub partition: Partition,
    pub cycle: CycleSettings,
    pub features: FeatureSettings,
    pub derive_unit: bool,
    pub coverage_distance: CoverageDistance,
}

impl Settings {
    pub fn of(args: &RecoverArgs) -> Result<Self> {
        Ok(Self {
            knn: knn_settings(&args.graph),
            seeds: seeds(&args.seeds),
            bars: Bars {
                min_bin_size: args.binning.min_bin_size,
                completeness: args.rescue.min_completeness,
                contamination: args.rescue.max_contamination,
                worth: args.rescue.worth_contamination,
            },
            lengths: Lengths {
                cutoff: args.binning.cutoff(),
                attach_floor: args.attach_floor(),
                max_bin: args.binning.max_bin_size,
            },
            partition: args.binning.partition,
            cycle: CycleSettings {
                refine_rounds: usize::from(!args.no_refine),
                partition_seeds: args.rescue.partition_seeds as usize,
                dissolve: args.rescue.dissolve == Switch::On,
                dissolve_rounds: args.rescue.dissolve_rounds as usize,
                dissolve_passes: args.rescue.dissolve_passes as usize,
                join: args.join,
                order: parse_order(&args.rescue.stage_order)?,
                split: args.rescue.split == Switch::On,
                peel: args.rescue.peel == Switch::On,
                passengers: args.rescue.passengers == Switch::On,
                replicons: args.rescue.replicons == Switch::On,
            },
            features: feature_settings(&args.graph),
            derive_unit: args.graph.variance_unit == Some(Unit::Auto),
            coverage_distance: args.graph.coverage_distance,
        })
    }
}
