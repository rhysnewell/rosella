use crate::cli::RunArgs;
use crate::clustering::graph_partition::Partition;
use crate::embedding::{features::FeatureSettings, knn::KnnSettings, metrics::DistanceSettings};
use crate::refine::rung::Bars;
use crate::seeds::Seeds;

pub fn seeds(params: &crate::cli::SeedParams) -> Seeds {
    Seeds {
        seed: params.seed,
        knn: crate::defaults::KNN_SEED,
        partition: crate::defaults::PARTITION_SEED,
    }
}

pub fn distance_settings() -> DistanceSettings {
    DistanceSettings {
        presence_fraction: crate::tuning::PRESENCE_FRACTION,
        aggregate_weight: None,
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
        distance: distance_settings(),
        link_weight: params.assembly_graph_weight as f32,
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Lengths {
    pub cutoff: usize,
    pub attach_floor: Option<usize>,
    pub min_bin: usize,
    pub max_bin: usize,
}

#[derive(Debug, Clone)]
pub struct CycleSettings {
    pub refine_rounds: usize,
    pub partition_seeds: usize,
    pub dissolve_rounds: usize,
    pub dissolve_passes: usize,
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
    pub epochs: usize,
    pub patience: usize,
    pub stretches: usize,
    pub fast_train: bool,
    pub marker_steps: Option<f64>,
    pub table_marker_steps: bool,
    pub gpu: bool,
    pub marker_cache: Option<std::path::PathBuf>,
    pub space_cache: Option<std::path::PathBuf>,
    pub split_contamination: Option<f64>,
}

impl Settings {
    pub fn of(args: &RunArgs, partition_seeds: usize, split_contamination: Option<f64>) -> Self {
        Self {
            knn: knn_settings(&args.graph),
            seeds: seeds(&args.seeds),
            bars: Bars {
                genome: crate::refine::rung::FALLBACK_GENOME,
                completeness: args.rescue.min_completeness,
                contamination: args.rescue.max_contamination,
                worth: args.rescue.worth_contamination,
            },
            lengths: Lengths {
                cutoff: args.binning.cutoff(),
                attach_floor: args.binning.attach_floor(),
                min_bin: args.binning.min_bin_size(),
                max_bin: args.binning.max_bin_size,
            },
            partition: args.binning.partition,
            cycle: CycleSettings {
                refine_rounds: usize::from(!args.rescue.no_refine),
                partition_seeds,
                dissolve_rounds: args.rescue.dissolve_rounds as usize,
                dissolve_passes: args.rescue.dissolve_passes as usize,
            },
            features: feature_settings(&args.graph),
            epochs: args.binning.epochs as usize,
            patience: args.binning.patience as usize,
            stretches: args.binning.stretches as usize,
            fast_train: args.binning.fast_train,
            marker_steps: (args.binning.marker_steps > 0.0).then_some(args.binning.marker_steps),
            table_marker_steps: args.binning.table_marker_steps,
            gpu: args.binning.gpu,
            marker_cache: args.markers.marker_cache.as_ref().map(Into::into),
            space_cache: args.binning.space_cache.as_ref().map(Into::into),
            split_contamination,
        }
    }

    // A later run on the same assembly reads the markers and the learned space from here.
    pub fn keep_in(mut self, folder: std::path::PathBuf) -> Self {
        self.marker_cache.get_or_insert_with(|| folder.clone());
        self.space_cache.get_or_insert(folder);
        self
    }
}
