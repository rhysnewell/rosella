use std::path::Path;

use anyhow::Result;
use log::debug;

use crate::embedding::{KNN_ASSEMBLY, knn::KnnGraph};
use crate::recover::recover_engine::RecoverEngine;

impl RecoverEngine {
    pub(super) fn probed(&mut self, contigs: &[usize], knn: &KnnGraph) -> Result<bool> {
        if let Some(path) = self.reports.knn_report.clone() {
            self.load(contigs.len())?;
            self.write_knn_report(contigs, &path)?;
            debug!("Wrote the kNN report to {}.", path.display());
            return Ok(true);
        }
        if let Some(paths) = &self.reports.reach_report {
            let groups =
                crate::refine::oracle::read_groups(&paths[0].to_string_lossy(), self.names())?;
            crate::embedding::reach::write(
                &paths[1],
                &groups,
                &self.features(),
                knn,
                self.lengths(),
                self.names(),
            )?;
            debug!("Wrote the reach report to {}.", paths[1].display());
            return Ok(true);
        }
        Ok(self.reports.partition_report.is_some())
    }

    fn write_knn_report(&self, contigs: &[usize], path: &Path) -> Result<()> {
        let knn = self.settings.knn;
        let seed = self.settings.seeds.knn;
        let combined =
            self.features()
                .knn_of(contigs, knn.neighbours, knn.candidates, seed, KNN_ASSEMBLY);
        let rho = self
            .features()
            .with_distance(self.settings.features.distance.composition_only())
            .knn_of(contigs, knn.neighbours, knn.candidates, seed, KNN_ASSEMBLY);
        let names = contigs
            .iter()
            .map(|index| self.names()[*index].as_str())
            .collect::<Vec<_>>();
        crate::embedding::knn::write_report(&[("combined", combined), ("rho", rho)], &names, path)
    }
}
