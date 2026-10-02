use std::{
    collections::{BTreeMap, HashMap, HashSet},
    path,
};

use anyhow::{Result, bail};
use log::{info, warn};
use ndarray::s;

use crate::bin_files::BinFiles;
use crate::quality::bases::Bases;
use crate::recover::recover_engine::{RecoverEngine, UNBINNED};

pub(super) const REPLICON_PREFIX: &str = "replicon_";
pub(super) const BIN_PREFIX: &str = "rosella_bin_";

#[derive(Default)]
pub(super) struct Written {
    pub(super) bases: HashMap<usize, Bases>,
    pub(super) labels: HashMap<usize, String>,
}

pub(super) struct Published {
    pub(super) bins: BTreeMap<usize, Vec<usize>>,
    pub(super) replicons: Vec<usize>,
    leftover: Vec<usize>,
}

impl RecoverEngine {
    pub(super) fn publish(
        &self,
        mut bins: BTreeMap<usize, Vec<usize>>,
        outliers: HashSet<usize>,
    ) -> Published {
        let departures = self.quality.departures(
            bins.values().map(Vec::as_slice),
            self.tnf_table.kmer_table.view(),
            self.coverage_table.table.slice(s![.., ..;2]),
        );
        let mut leftover = outliers.into_iter().collect::<Vec<_>>();
        for contigs in bins.values_mut() {
            contigs.retain(|contig| {
                departures.replicons.binary_search(contig).is_err()
                    && departures.passengers.binary_search(contig).is_err()
            });
        }
        leftover.extend(departures.passengers);
        let replicons = departures.replicons;
        bins.retain(|_, contigs| {
            let bp = contigs
                .iter()
                .map(|contig| self.coverage_table.contig_lengths[*contig])
                .sum::<usize>();
            if bp < self.min_bin_size {
                leftover.extend(contigs.iter().copied());
            }
            bp >= self.min_bin_size
        });
        Published {
            bins,
            replicons,
            leftover,
        }
    }

    // Keyed on contig name rather than position, because the clustering indexes the coverage
    // table and the length filter has already shortened it.
    fn placements<'a>(
        &'a self,
        published: &Published,
    ) -> HashMap<&'a str, (usize, Option<Target>)> {
        let name = |contig: &usize| self.coverage_table.contig_names[*contig].as_str();
        published
            .bins
            .iter()
            .flat_map(|(bin, contigs)| {
                contigs
                    .iter()
                    .map(move |contig| (contig, Some(Target::Bin(*bin))))
            })
            .chain(
                published
                    .replicons
                    .iter()
                    .enumerate()
                    .map(|(at, contig)| (contig, Some(Target::Replicon(at + 1)))),
            )
            .chain(published.leftover.iter().map(|contig| (contig, None)))
            .map(|(contig, target)| (name(contig), (*contig, target)))
            .collect()
    }

    pub(super) fn write_clusters(&self, published: &Published) -> Result<Written> {
        let placed = self.placements(published);
        let mut held = Written::default();
        let directory = path::Path::new(&self.output_directory);
        let mut files = BinFiles::new(|target: &Target| {
            directory.join(format!(
                "{BIN_PREFIX}{}.{}",
                target.name(),
                crate::defaults::FASTA_EXTENSION
            ))
        });

        let mut singles = 0;
        let mut unrecognised = 0;
        let mut read = 0;
        let mut written = 0;
        let progress = crate::progress::spinning(crate::progress::Stage::WritingBins);

        crate::kmers::pipelined(
            &self.assembly,
            |record| Ok((record.id().to_vec(), record.seq().into_owned())),
            |chunk| {
                for (id, sequence) in chunk {
                    read += 1;
                    let found = placed.get(crate::contig_id(&id)?);
                    let long = sequence.len() >= self.min_contig_size;
                    let target = match found {
                        Some((_, Some(target))) if long => *target,
                        _ => {
                            unrecognised += usize::from(long && found.is_none());
                            self.leftover(sequence.len(), &mut singles)
                        }
                    };

                    if let Some((contig, placement)) = found {
                        if placement.is_some() {
                            held.bases.insert(*contig, Bases::count(&sequence));
                        }
                        if self.reports_markers(*contig) {
                            held.labels
                                .insert(*contig, format!("{BIN_PREFIX}{}", target.name()));
                        }
                    }
                    files.write(&target, &id, &sequence)?;
                    written += 1;
                    if written % PROGRESS_EVERY == 0 {
                        progress.set_message(format!("{written} contigs"));
                    }
                }
                Ok(())
            },
        )?;
        progress.finish_and_clear();
        let n_bins = files.finish()?;

        if unrecognised > 0 {
            warn!(
                "{} contigs in the assembly were not in the coverage table and went unbinned",
                unrecognised
            );
        }
        if written != read {
            bail!(
                "{} of {} assembly contigs were written. Every contig belongs in a bin, in \
                 rosella_bin_unbinned or in rosella_bin_small_unbinned, so a shortfall means \
                 contigs were dropped",
                written,
                read
            );
        }
        info!(
            "Wrote {written} contigs into {n_bins} bins in {}.",
            self.output_directory
        );

        Ok(held)
    }

    fn leftover(&self, contig_length: usize, singles: &mut usize) -> Target {
        if contig_length >= self.min_bin_size {
            *singles += 1;
            return Target::Single(*singles);
        }
        if contig_length < self.min_contig_size {
            return Target::Small;
        }
        Target::Unbinned
    }
}

const PROGRESS_EVERY: usize = 4096;

#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
enum Target {
    Bin(usize),
    Replicon(usize),
    Single(usize),
    Unbinned,
    Small,
}

impl Target {
    fn name(self) -> String {
        match self {
            Self::Bin(bin) => bin.to_string(),
            Self::Replicon(at) => format!("{REPLICON_PREFIX}{at}"),
            Self::Single(at) => format!("single_contig_{at}"),
            Self::Unbinned => UNBINNED.to_string(),
            Self::Small => format!("small_{UNBINNED}"),
        }
    }
}
