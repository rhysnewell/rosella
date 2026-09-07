use anyhow::Result;

use crate::{
    embedding::{
        manifold::GraphWeights,
        metrics::{
            Combination, CompositionMetric, CoverageAggregation, CoverageBand, DistanceSettings,
            Views,
        },
        spectral::SpectralInit,
        umap::EmbedOverrides,
    },
    kmers::kmer_counting::KmerFrequencyTable,
    seeds::Seeds,
};

pub fn embed_overrides(overrides: &crate::cli::EmbeddingOverrides) -> EmbedOverrides {
    EmbedOverrides {
        a: overrides.umap_a,
        b: overrides.umap_b,
        min_dist: overrides.min_dist,
        spread: overrides.spread,
        n_components: overrides.n_components,
        n_epochs: overrides.n_epochs,
        length_weight: overrides.length_weight,
        spectral_init: SpectralInit::parse(&overrides.spectral_init)
            .expect("clap restricts the value"),
        report_preservation: overrides.report_preservation,
        knn_candidates: overrides.knn_candidates,
        graph_weights: GraphWeights::parse(&overrides.graph_weights)
            .expect("clap restricts the value"),
    }
}

pub fn seeds(seed: u64, overrides: &crate::cli::SeedOverrides) -> Seeds {
    Seeds {
        knn: overrides.knn.unwrap_or(seed),
        init: overrides.init.unwrap_or(seed),
        layout: overrides.layout.unwrap_or(seed),
        sample: overrides.sample.unwrap_or(seed),
        partition: overrides.partition.unwrap_or(seed),
    }
}

pub fn distance_settings(distance: &crate::cli::DistanceParams) -> Result<DistanceSettings> {
    let views = Views::parse(&distance.embedding_views).ok_or_else(|| {
        anyhow::anyhow!("--embedding-views combined cannot be listed beside another view")
    })?;
    Ok(DistanceSettings {
        aggregation: CoverageAggregation::parse(&distance.coverage_aggregation)
            .expect("clap restricts the value"),
        length_scaled_variance: distance.length_scaled_variance,
        views,
        aggregate_weight: distance.aggregate_weight,
        combination: Combination::parse(&distance.distance_combination)
            .expect("clap restricts the value"),
        presence_fraction: distance.presence_fraction,
        coverage_band: CoverageBand::parse(&distance.coverage_band)
            .expect("clap restricts the value"),
        composition: CompositionMetric::parse(&distance.composition_metric)
            .expect("clap restricts the value"),
        composition_scale: 1.0,
    })
}

/// The composition metric picks the table it compares, so both engines transform through here
/// rather than each choosing its own.
pub fn transform_table(
    table: &mut KmerFrequencyTable,
    metric: CompositionMetric,
    contig_lengths: &[usize],
) -> Result<()> {
    match metric {
        CompositionMetric::Hellinger => {
            crate::kmers::transform::sqrt_frequencies(&mut table.kmer_table);
            Ok(())
        }
        CompositionMetric::TetraZ => {
            crate::kmers::transform::tetra_z(&mut table.kmer_table, table.kmer_size)
        }
        _ => table.clr(contig_lengths),
    }
}
