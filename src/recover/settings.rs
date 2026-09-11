use anyhow::Result;

use crate::{
    embedding::{
        manifold::GraphWeights,
        metrics::{Combination, CompositionMetric, CoverageAggregation, DistanceSettings},
        umap::EmbedOverrides,
    },
    seeds::Seeds,
};

pub fn embed_overrides(overrides: &crate::cli::EmbeddingOverrides) -> EmbedOverrides {
    EmbedOverrides {
        a: overrides.umap_a,
        b: overrides.umap_b,
        min_dist: overrides.min_dist,
        spread: overrides.spread,
        knn_candidates: overrides.knn_candidates,
        graph_weights: GraphWeights::parse(&overrides.graph_weights)
            .expect("clap restricts the value"),
    }
}

pub fn seeds(seed: u64, overrides: &crate::cli::SeedOverrides) -> Seeds {
    Seeds {
        knn: overrides.knn.unwrap_or(seed),
        sample: overrides.sample.unwrap_or(seed),
        partition: overrides.partition.unwrap_or(seed),
    }
}

pub fn distance_settings(distance: &crate::cli::DistanceParams) -> Result<DistanceSettings> {
    Ok(DistanceSettings {
        aggregation: CoverageAggregation::parse(&distance.coverage_aggregation)
            .expect("clap restricts the value"),
        length_scaled_variance: distance.length_scaled_variance,
        combination: Combination::parse(&distance.distance_combination)
            .expect("clap restricts the value"),
        presence_fraction: distance.presence_fraction,
        composition: CompositionMetric::parse(&distance.composition_metric)
            .expect("clap restricts the value"),
        composition_scale: 1.0,
        aggregate_weight: None,
    })
}

/// No contig length statistic separates the two arms: the pool is worth +16 t1 on the highest
/// N50 single set and -55 on CAMI I high. Total assembly bp splits the sign on all 14 measured.
const LARGE_ASSEMBLY_BP: usize = 350_000_000;

pub const DISSOLVE_NAMES: [&str; 3] = ["auto", "on", "off"];

pub fn dissolve(choice: &str, lengths: &[usize]) -> bool {
    match choice {
        "on" => true,
        "off" => false,
        _ => {
            let total = lengths.iter().sum::<usize>();
            let on = total < LARGE_ASSEMBLY_BP;
            log::info!(
                "Assembly {} Mbp past the filter, rescue pool {}",
                total / 1_000_000,
                if on { "on" } else { "off" }
            );
            on
        }
    }
}
