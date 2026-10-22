use std::collections::{BTreeMap, HashSet};

use anyhow::Result;
use log::warn;

use super::bin_writer::Published;
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::cut_report::CutLog;
use crate::refine::owners::owners;

impl RecoverEngine {
    pub(super) fn publish_traced(
        &self,
        bins: BTreeMap<usize, Vec<usize>>,
        outliers: HashSet<usize>,
        cuts: Option<CutLog>,
    ) -> Result<Published> {
        let (Some(path), Some(mut log)) = (&self.reports.cut_report, cuts) else {
            return self.publish(bins, outliers);
        };
        let before = bins.values().cloned().collect::<Vec<_>>();
        let published = self.publish(bins, outliers)?;
        let lengths = &self.tables.coverage.contig_lengths;
        let after = owners(
            published
                .bins
                .iter()
                .map(|(label, members)| (*label, members)),
        );
        log.diff("publish", &before, &after, |contig| lengths[contig]);
        if let Err(error) = log.write(
            path,
            &published.bins,
            &self.quality,
            self.settings.bars.worth,
            lengths,
            &self.tables.coverage.contig_names,
        ) {
            warn!("Cut report: {error} ({})", path.display());
        }
        Ok(published)
    }
}
