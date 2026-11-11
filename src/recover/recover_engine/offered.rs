use anyhow::Result;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::knn::KnnGraph;
use crate::recover::ladder::{Judge, combine_in};
use crate::recover::recover_engine::RecoverEngine;
use crate::recover::recover_engine::weight::{Settled, Worth};
use crate::refine::splitter::Refined;

pub const PARTITION: &str = "partition";

fn spaced() -> bool {
    std::env::var_os("ROSELLA_HOME_SPACE").is_some()
}

impl RecoverEngine {
    // Finished bin sets stand in for a partition, so their contest is the first pass. The run's
    // own set comes first and teaches the markers, as given bins do.
    pub(super) fn offered_partition(&self, contigs: &[usize]) -> Result<Settled> {
        let long = self.long_end(contigs);
        let arms = self.candidate_arms(&contigs[..long]);
        let own = || arms[0].cluster_map.values().map(Vec::as_slice);
        let spaced = match spaced() {
            true => Some(self.weighed(&arms[0], contigs)?),
            false => None,
        };
        self.quality
            .learn_checkm(own(), self.settings.bars.quality());
        self.quality.learn_multicopy(own());
        let voters = (0..arms.len()).collect::<Vec<_>>();
        let space = spaced.as_ref().map(|(_, _, knn)| knn);
        let partitioning = self.contest("contest", &arms, &voters, contigs, space);
        let (weight, graph, knn) = match spaced {
            Some(found) => found,
            None => self.weighed(&partitioning, contigs)?,
        };
        self.learn_multicopy(&partitioning, contigs);
        Ok(Settled {
            graph,
            knn,
            worth: Worth {
                kept: self.pass_worth(&partitioning, contigs),
                spread: 0.0,
            },
            partitioning,
            weight,
            parked: contigs[long..].to_vec(),
        })
    }

    // The refine cycle reads markers alone, which cannot see a genome's own marker-free contigs,
    // so its bins go back to the vote beside every offered set.
    pub(super) fn revote(&self, refined: &mut Refined, contigs: &[usize], knn: &KnnGraph) {
        let long = &contigs[..self.long_end(contigs)];
        let mut arms = vec![Partitioning {
            cluster_map: std::mem::take(&mut refined.bins),
            outliers: std::mem::take(&mut refined.unbinned),
            score: None,
            arm: self.settings.partition,
            seed: self.settings.seeds.partition,
        }];
        arms.extend(self.candidate_arms(long));
        // The run's own finished set speaks with the bins refined from it, not as one more binner.
        let voters = std::iter::once(0)
            .chain(0..arms.len() - 1)
            .collect::<Vec<_>>();
        let space = spaced().then_some(knn);
        let combined = self.contest("revote", &arms, &voters, contigs, space);
        refined.bins = combined.cluster_map;
        refined.unbinned = combined.outliers;
    }

    fn contest(
        &self,
        name: &str,
        arms: &[Partitioning],
        voters: &[usize],
        contigs: &[usize],
        space: Option<&KnnGraph>,
    ) -> Partitioning {
        let features = self.features();
        let counted = crate::markers::Fair(&self.quality);
        let judge = Judge {
            quality: &counted,
            contigs,
            lengths: self.lengths(),
            bars: self.settings.bars,
            depth: &features,
            markers: &self.quality,
        };
        let combined = combine_in(arms, voters, &judge, space);
        if let Some(report) = &self.stages {
            let names = &self.tables.coverage.contig_names;
            report.combined(name, arms, &combined, contigs, names);
        }
        combined.partitioning
    }
}
