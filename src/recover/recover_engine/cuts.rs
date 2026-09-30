use std::collections::{HashMap, HashSet};

use log::warn;

use crate::recover::bin_writer::Published;
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::cut_report::{CutLog, owners};

impl RecoverEngine {
    pub(super) fn publish_traced(
        &self,
        bins: HashMap<usize, HashSet<usize>>,
        outliers: HashSet<usize>,
        cuts: Option<CutLog>,
    ) -> Published {
        let (Some(path), Some(mut log)) = (&self.cut_report, cuts) else {
            return self.publish(bins, outliers);
        };
        let before = bins.values().cloned().collect::<Vec<_>>();
        let published = self.publish(bins, outliers);
        let lengths = &self.coverage_table.contig_lengths;
        let after = owners(
            published
                .bins
                .iter()
                .map(|(label, members)| (*label, members)),
        );
        log.diff("publish", &before, &after, |contig| lengths[contig]);
        let finals = published
            .bins
            .iter()
            .map(|(label, members)| {
                let mut members = members.iter().copied().collect::<Vec<_>>();
                members.sort_unstable();
                (*label, members)
            })
            .collect();
        if let Err(error) = log.write(
            path,
            &finals,
            &self.quality,
            self.worth,
            lengths,
            &self.coverage_table.contig_names,
        ) {
            warn!("No cut report at {}: {error}", path.display());
        }
        published
    }
}
