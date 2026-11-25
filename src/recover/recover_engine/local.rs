use std::sync::{Arc, LazyLock};

use anyhow::Result;
use log::debug;
use ndarray::ArrayView1;
use rayon::prelude::*;

use crate::embedding::features::{FeatureSettings, row_slice};
use crate::embedding::metrics::{Centred, LocalWeights, rho_between};
use crate::recover::recover_engine::RecoverEngine;

pub(super) static LOCAL: LazyLock<bool> =
    LazyLock::new(|| std::env::var_os("ROSELLA_LOCAL_WEIGHT").is_some());

const SAMPLE: usize = 2_000;
// Two halves of one contig read as foreign once in twenty (Fisher 1925).
const OWN: f64 = 0.95;

impl RecoverEngine {
    // Where a contig's depth neighbours fall into several composition groups, depth cannot part
    // them, so its pairs lean on composition by the share of those neighbours that are foreign.
    pub(super) fn local_weights(
        &self,
        weight: f64,
        contigs: &[usize],
    ) -> Result<Option<LocalWeights>> {
        let Some(foreign) = self.foreign_past(contigs)? else {
            return Ok(None);
        };
        let mut depth = self.settings.features.clone();
        depth.distance.aggregate_weight = Some(1.0);
        depth.distance.local = None;
        let knn = self.tables.features(&depth).knn_of(
            contigs,
            self.settings.knn.neighbours,
            self.settings.knn.candidates,
            self.settings.seeds.knn,
            "local_weight",
        );
        let composition = self
            .tables
            .features(&FeatureSettings {
                distance: depth.distance.composition_only(),
                ..depth.clone()
            })
            .prepared(contigs);
        let shares = (0..contigs.len())
            .into_par_iter()
            .map(|row| {
                let near = knn
                    .indices
                    .row(row)
                    .iter()
                    .filter(|at| **at != u32::MAX)
                    .map(|at| *at as usize)
                    .collect::<Vec<_>>();
                let far = near
                    .iter()
                    .filter(|at| composition.distance(row, **at) > foreign)
                    .count();
                far as f64 / near.len().max(1) as f64
            })
            .collect::<Vec<_>>();
        let mut local = vec![weight; self.tables.tnf.kmer_table.nrows()];
        for (contig, share) in contigs.iter().zip(&shares) {
            local[*contig] = weight * (1.0 - share);
        }
        debug!(
            "Local weight: mean foreign share {:.3} of depth neighbours past {foreign:.3}, weight {weight:.3}",
            shares.iter().sum::<f64>() / shares.len().max(1) as f64
        );
        Ok(Some(LocalWeights(Arc::new(local))))
    }

    fn foreign_past(&self, contigs: &[usize]) -> Result<Option<f64>> {
        let lengths = self.lengths();
        let pool = contigs
            .iter()
            .copied()
            .filter(|contig| lengths[*contig] >= 2 * self.settings.lengths.cutoff)
            .collect::<Vec<_>>();
        if pool.is_empty() {
            return Ok(None);
        }
        let chosen = crate::seeds::sample_positions(pool.len(), SAMPLE, self.settings.seeds.seed)
            .into_iter()
            .map(|at| pool[at])
            .collect::<Vec<_>>();
        let names = chosen
            .iter()
            .map(|contig| self.tables.coverage.contig_names[*contig].as_str())
            .collect::<Vec<_>>();
        let mut halves = self.tables.halves.composition(&names)?;
        if let Some(space) = &self.tables.learned {
            let coverage = chosen
                .iter()
                .map(|contig| row_slice(&self.tables.coverage.table, *contig))
                .collect::<Vec<_>>();
            let sizes = chosen
                .iter()
                .map(|contig| lengths[*contig])
                .collect::<Vec<_>>();
            let Some(embedded) = space.halves(halves, &coverage, &sizes, self.settings.seeds.seed)
            else {
                return Ok(None);
            };
            halves = embedded;
        }
        let [first, second] = halves;
        let mut own = (0..chosen.len())
            .map(|row| apart(first.row(row), second.row(row)))
            .collect::<Vec<_>>();
        own.sort_by(f64::total_cmp);
        Ok(Some(own[((own.len() - 1) as f64 * OWN) as usize]))
    }
}

fn apart(first: ArrayView1<f64>, second: ArrayView1<f64>) -> f64 {
    let [first, second] = [first, second].map(|row| Centred::new(&row.to_vec()));
    rho_between(&first, &second)
}
