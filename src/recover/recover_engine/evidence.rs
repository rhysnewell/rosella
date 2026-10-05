use std::sync::Arc;

use anyhow::Result;
use log::{info, warn};
use rand::{Rng, SeedableRng, rngs::StdRng};

use crate::embedding::features::row_slice;
use crate::embedding::metrics::likelihood::{CoverageDistance, Evidence};
use crate::embedding::metrics::{Centred, DistanceSettings};
use crate::embedding::weight::composition::{self, Halves};
use crate::embedding::weight::depth::{self, POOL, QUERIES};
use crate::recover::recover_engine::RecoverEngine;

const MATES: usize = 10;

impl RecoverEngine {
    pub(super) fn evidence(&self, contigs: &[usize]) -> Result<Option<Evidence>> {
        if self.settings.coverage_distance == CoverageDistance::Overlap {
            return Ok(None);
        }
        let _timer = crate::timing::scope("weight");
        let lengths = self.lengths();
        let long = contigs
            .iter()
            .copied()
            .filter(|contig| lengths[*contig] >= self.settings.lengths.cutoff)
            .collect::<Vec<_>>();
        let pool = crate::seeds::sample_positions(
            long.len(),
            POOL.min(long.len()),
            self.settings.seeds.seed,
        )
        .into_iter()
        .map(|at| long[at])
        .collect::<Vec<_>>();
        let rows = |table| {
            pool.iter()
                .map(|contig| row_slice(table, *contig))
                .collect::<Vec<_>>()
        };
        let split = pool
            .iter()
            .copied()
            .filter(|contig| lengths[*contig] >= 2 * self.settings.lengths.cutoff)
            .take(QUERIES)
            .collect::<Vec<_>>();
        let names = split
            .iter()
            .map(|contig| self.names()[*contig].as_str())
            .collect::<Vec<_>>();
        let [first, second] = self
            .tables
            .halves
            .composition(&names, self.tables.tnf.kmer_size())?;
        let halves = Halves {
            first: (0..split.len()).map(|at| row_slice(&first, at)).collect(),
            second: (0..split.len()).map(|at| row_slice(&second, at)).collect(),
            lengths: split.iter().map(|contig| lengths[*contig] as f64).collect(),
        };
        let model = depth::fit(
            &rows(&self.tables.coverage.table),
            &pool
                .iter()
                .map(|contig| lengths[*contig] as f64)
                .collect::<Vec<_>>(),
            &rows(&self.tables.tnf.kmer_table),
            &halves,
            self.settings.seeds.seed,
        );
        match &model {
            None => {
                warn!("Too few long contigs to fit the depth model, so depth is scored by overlap.")
            }
            Some(model) => {
                let fits = model
                    .samples
                    .iter()
                    .map(|sample| match sample {
                        Some(sample) => format!(
                            "{:.0}:{:.3}:{:.2}",
                            sample.read_length,
                            sample.spread.sqrt(),
                            sample.absent.exp()
                        ),
                        None => "-".to_string(),
                    })
                    .collect::<Vec<_>>();
                info!(
                    "Depth model from {} contigs, prior odds of a different genome {:.2}, read \
                     length, log depth spread and absent share per sample [{}].",
                    pool.len(),
                    model.prior.exp(),
                    fits.join(", ")
                );
            }
        }
        Ok(model)
    }

    pub(super) fn refit(
        &self,
        near_complete: &[Vec<usize>],
        distance: &DistanceSettings,
    ) -> Option<DistanceSettings> {
        let evidence = distance.evidence.as_ref()?;
        let model = evidence.composition.as_ref()?;
        let lengths = self.lengths();
        let bins = near_complete
            .iter()
            .map(|bin| {
                bin.iter()
                    .copied()
                    .filter(|contig| lengths[*contig] >= self.settings.lengths.cutoff)
                    .collect::<Vec<_>>()
            })
            .filter(|bin| bin.len() > 1)
            .collect::<Vec<_>>();
        let members = bins
            .iter()
            .enumerate()
            .flat_map(|(at, bin)| bin.iter().map(move |contig| (at, *contig)))
            .collect::<Vec<_>>();
        let seed = self.settings.seeds.seed;
        let mut rng = StdRng::seed_from_u64(seed);
        let centred = |contig: usize| Centred::new(row_slice(&self.tables.tnf.kmer_table, contig));
        let mut points = Vec::new();
        for at in crate::seeds::sample_positions(members.len(), QUERIES.min(members.len()), seed) {
            let (bin, contig) = members[at];
            let own = centred(contig);
            for _ in 0..MATES {
                let mate = bins[bin][rng.random_range(0..bins[bin].len())];
                if mate != contig {
                    points.push(composition::point(
                        &own,
                        &centred(mate),
                        lengths[contig] as f64,
                        lengths[mate] as f64,
                    ));
                }
            }
        }
        let refitted = Evidence {
            composition: Some(composition::refit(model, &points)?),
            ..(**evidence).clone()
        };
        Some(DistanceSettings {
            evidence: Some(Arc::new(refitted)),
            ..distance.clone()
        })
    }
}
