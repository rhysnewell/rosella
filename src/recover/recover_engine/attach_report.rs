use std::collections::{HashMap, HashSet};

use anyhow::Result;
use rayon::prelude::*;

use crate::quality::Scorer;
use crate::recover::recover_engine::RecoverEngine;
use crate::recover::recover_engine::attach::Searched;

impl RecoverEngine {
    pub(super) fn write_attach_report(
        &self,
        path: &std::path::Path,
        bins: &HashMap<usize, HashSet<usize>>,
        searched: &[Searched],
        refused: &HashSet<usize>,
    ) -> Result<()> {
        use std::io::Write;
        let members = bins
            .iter()
            .map(|(bin, contigs)| {
                let mut contigs = contigs.iter().copied().collect::<Vec<_>>();
                contigs.sort_unstable();
                (*bin, contigs)
            })
            .collect::<HashMap<_, _>>();
        let rows = searched
            .iter()
            .flat_map(|band| {
                band.proposals
                    .iter()
                    .enumerate()
                    .map(move |(at, proposal)| {
                        let chance = band.chances.as_ref().map(|chances| chances[at]);
                        (proposal, chance, band.taken)
                    })
            })
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|((contig, best), chance, taken)| {
                let chance = chance.map_or("NA".to_string(), |chance| format!("{chance:.3}"));
                let taken = u8::from(taken);
                let Some((bin, share)) = best else {
                    return format!(
                        "{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t{chance}\tNA\t{taken}",
                        self.coverage_table.contig_names[*contig],
                        self.coverage_table.contig_lengths[*contig]
                    );
                };
                let rest = &members[bin];
                let lengths = &self.coverage_table.contig_lengths;
                let anchor = rest.iter().max_by_key(|at| (lengths[**at], **at)).copied();
                let mut with = rest.clone();
                with.insert(with.partition_point(|at| at < contig), *contig);
                let flag =
                    |seen: Option<bool>| seen.map_or("NA", |seen| if seen { "1" } else { "0" });
                format!(
                    "{}\t{}\t{bin}\t{}\t{share:.3}\t{}\t{}\t{:.3}\t{}\t{chance}\t{}\t{taken}",
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
            "contig\tlength\tbin\tanchor\tshare\trepeats_whole\trepeats_any\tcomplete\trefused\tchance\t\
             repeats_place\ttaken"
        )?;
        for row in rows {
            writeln!(sink, "{row}")?;
        }
        Ok(())
    }
}
