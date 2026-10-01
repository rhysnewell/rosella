use std::collections::{BTreeMap, HashSet};

use anyhow::Result;
use rayon::prelude::*;

use crate::quality::Scorer;
use crate::recover::recover_engine::RecoverEngine;
use crate::recover::recover_engine::attach::Searched;

impl RecoverEngine {
    pub(super) fn write_attach_report(
        &self,
        path: &std::path::Path,
        bins: &BTreeMap<usize, Vec<usize>>,
        searched: &[Searched],
        refused: &HashSet<usize>,
    ) -> Result<()> {
        use std::io::Write;
        let rows = searched
            .iter()
            .flat_map(|band| {
                band.proposals
                    .iter()
                    .map(move |proposal| (proposal, band.bar))
            })
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|((contig, best), bar)| {
                let Some((bin, share)) = best else {
                    return format!(
                        "{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t{bar:.3}",
                        self.coverage_table.contig_names[*contig],
                        self.coverage_table.contig_lengths[*contig]
                    );
                };
                let rest = &bins[bin];
                let lengths = &self.coverage_table.contig_lengths;
                let anchor = rest.iter().max_by_key(|at| (lengths[**at], **at)).copied();
                let mut with = rest.clone();
                with.insert(with.partition_point(|at| at < contig), *contig);
                let flag =
                    |seen: Option<bool>| seen.map_or("NA", |seen| if seen { "1" } else { "0" });
                format!(
                    "{}\t{}\t{bin}\t{}\t{share:.3}\t{}\t{}\t{:.3}\t{}\t{}\t{bar:.3}",
                    self.coverage_table.contig_names[*contig],
                    self.coverage_table.contig_lengths[*contig],
                    anchor.map_or("NA", |at| self.coverage_table.contig_names[at].as_str()),
                    flag(self.quality.repeats(&with, *contig)),
                    flag(self.quality.repeats_any(&with, *contig)),
                    self.quality.score(rest).completeness / 100.0,
                    u8::from(refused.contains(bin)),
                    flag(self.quality.repeats_in_place(&with, *contig)),
                )
            })
            .collect::<Vec<_>>();
        let mut sink = std::io::BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            sink,
            "contig\tlength\tbin\tanchor\tshare\trepeats_whole\trepeats_any\tcomplete\trefused\t\
             repeats_place\tbar"
        )?;
        for row in rows {
            writeln!(sink, "{row}")?;
        }
        sink.flush()?;
        Ok(())
    }
}
