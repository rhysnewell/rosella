use std::{
    collections::{BTreeMap, HashMap, HashSet},
    path,
};

use anyhow::Result;
use log::{debug, info, warn};
use ndarray::s;

use crate::bin_files::BinFiles;
use crate::kmers::kmer_counting::holds_clean_kmer;
use crate::quality::Scorer;
use crate::quality::bases::Bases;
use crate::recover::circular;
use crate::recover::fold::{Role, folds};
use crate::recover::recover_engine::{RecoverEngine, UNBINNED};

pub(super) const REPLICON_PREFIX: &str = "replicon_";
pub(super) const EUKARYOTE_PREFIX: &str = "eukaryote_";

#[derive(Default)]
pub(super) struct Written {
    pub(super) bases: HashMap<usize, Bases>,
    pub(super) labels: HashMap<usize, String>,
}

pub(super) struct Published {
    pub(super) bins: BTreeMap<usize, Vec<usize>>,
    pub(super) replicons: Vec<usize>,
    pub(super) origins: Vec<Option<usize>>,
    eukaryotes: HashSet<usize>,
    leftover: Vec<usize>,
}

impl Published {
    pub(super) fn genome(&self, bin: usize) -> Target {
        if self.eukaryotes.contains(&bin) {
            return Target::Eukaryote(bin);
        }
        Target::Bin(bin)
    }
}

impl RecoverEngine {
    pub(super) fn publish(
        &self,
        mut bins: BTreeMap<usize, Vec<usize>>,
        outliers: HashSet<usize>,
    ) -> Result<Published> {
        let names = &self.tables.coverage.contig_names;
        let mut departures = self.quality.departures(
            bins.values().map(Vec::as_slice),
            |members| {
                let members = members
                    .iter()
                    .map(|contig| names[*contig].as_str())
                    .collect::<Vec<_>>();
                self.tables
                    .halves
                    .whole(&members, self.tables.tnf.kmer_size())
            },
            self.tables.coverage.table.slice(s![.., ..;2]),
        )?;
        if !self.settings.cycle.passengers {
            departures.passengers.clear();
        }
        if !self.settings.cycle.replicons {
            departures.replicons.clear();
        }
        let mut leftover = outliers.into_iter().collect::<Vec<_>>();
        let mut origins = vec![None; departures.replicons.len()];
        for (bin, contigs) in bins.iter_mut() {
            contigs.retain(|contig| match departures.replicons.binary_search(contig) {
                Ok(at) => {
                    origins[at] = Some(*bin);
                    false
                }
                Err(_) => departures.passengers.binary_search(contig).is_err(),
            });
        }
        leftover.extend(departures.passengers);
        let replicons = departures.replicons;
        bins.retain(|_, contigs| {
            let bp = contigs
                .iter()
                .map(|contig| self.tables.coverage.contig_lengths[*contig])
                .sum::<usize>();
            if bp < self.settings.bars.min_bin_size {
                leftover.extend(contigs.iter().copied());
            }
            bp >= self.settings.bars.min_bin_size
        });
        let eukaryotes = self
            .quality
            .eukaryotic(
                bins.values().map(Vec::as_slice),
                self.settings.bars.quality(),
            )
            .into_iter()
            .zip(bins.keys())
            .filter_map(|(eukaryotic, bin)| eukaryotic.then_some(*bin))
            .collect::<HashSet<_>>();
        if !eukaryotes.is_empty() {
            debug!("Eukaryotic bins: {}", eukaryotes.len());
        }
        let moved = self.fold(&mut bins, &replicons, &eukaryotes);
        let origins = origins
            .into_iter()
            .map(|bin| {
                let bin = bin?;
                let bin = moved.get(&bin).copied().unwrap_or(bin);
                bins.contains_key(&bin).then_some(bin)
            })
            .collect();
        Ok(Published {
            bins,
            replicons,
            origins,
            eukaryotes,
            leftover,
        })
    }

    fn fold(
        &self,
        bins: &mut BTreeMap<usize, Vec<usize>>,
        replicons: &[usize],
        eukaryotes: &HashSet<usize>,
    ) -> HashMap<usize, usize> {
        let keys = bins.keys().copied().collect::<Vec<_>>();
        let groups = bins
            .values()
            .map(Vec::as_slice)
            .chain(replicons.iter().map(std::slice::from_ref))
            .collect::<Vec<_>>();
        let roles = bins
            .iter()
            .map(|(bin, members)| Role::of(self.quality.score(members), eukaryotes.contains(bin)))
            .chain(replicons.iter().map(|_| Role::Other))
            .collect::<Vec<_>>();
        let folds = folds(
            &groups,
            &roles,
            self.tables.coverage.table.slice(s![.., ..;2]),
            &self.tables.coverage.contig_lengths,
        );
        for (donor, host) in &folds {
            let moved = bins.remove(&keys[*donor]).unwrap_or_default();
            let into = bins.entry(keys[*host]).or_default();
            into.extend(moved);
            into.sort_unstable();
        }
        if !folds.is_empty() {
            debug!("Folded marker-free bins: {}", folds.len());
        }
        folds
            .into_iter()
            .map(|(donor, host)| (keys[donor], keys[host]))
            .collect()
    }

    // Keyed on contig name rather than position, because the clustering indexes the coverage
    // table and the length filter has already shortened it.
    fn placements<'a>(
        &'a self,
        published: &Published,
    ) -> HashMap<&'a str, (usize, Option<Target>)> {
        let name = |contig: &usize| self.tables.coverage.contig_names[*contig].as_str();
        published
            .bins
            .iter()
            .flat_map(|(bin, contigs)| {
                let target = published.genome(*bin);
                contigs.iter().map(move |contig| (contig, Some(target)))
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

    pub(super) fn write_clusters(&self, published: &Published, floor: usize) -> Result<Written> {
        let placed = self.placements(published);
        let mut held = Written::default();
        let directory = path::Path::new(&self.output_directory);
        let mut files = BinFiles::new(|target: &Target| {
            directory.join(format!(
                "{}{}.{}",
                self.prefix,
                target.name(),
                crate::defaults::FASTA_EXTENSION
            ))
        });
        let universe = self.given.as_ref().map(|given| &given.universe);

        let kmer_size = self.tables.tnf.kmer_size();
        let mut singles = 0;
        let mut circles = 0;
        let mut unrecognised = 0;
        let mut written = 0;
        let progress = crate::progress::spinning(crate::progress::Stage::WritingBins);

        crate::kmers::pipelined(
            &self.tables.assembly,
            |record| Ok((record.id().to_vec(), record.seq().into_owned())),
            |chunk| {
                for (id, sequence) in chunk {
                    let name = crate::contig_id(&id)?;
                    if universe.is_some_and(|universe| !universe.contains(name)) {
                        continue;
                    }
                    let found = placed.get(name);
                    let long = sequence.len() >= floor;
                    let target = match found {
                        Some((_, Some(target))) if long => *target,
                        _ => {
                            unrecognised += usize::from(
                                long && found.is_none() && holds_clean_kmer(&sequence, kmer_size),
                            );
                            match self.circular(name, &id, &sequence) {
                                true => {
                                    circles += 1;
                                    Target::Circular(circles)
                                }
                                false => self.leftover(sequence.len(), floor, &mut singles),
                            }
                        }
                    };

                    if let Some((contig, placement)) = found {
                        if placement.is_some() {
                            held.bases.insert(*contig, Bases::count(&sequence));
                        }
                        if self.reports_markers(*contig) {
                            held.labels
                                .insert(*contig, format!("{}{}", self.prefix, target.name()));
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
        crate::progress::done(&progress);
        let n_bins = files.finish()?;

        if unrecognised > 0 {
            warn!("Contigs not in the coverage table: {unrecognised} (unbinned)");
        }
        if let Some(missing) = universe.map(|universe| universe.len().saturating_sub(written))
            && missing > 0
        {
            warn!("Bin contigs not in the assembly: {missing}");
        }
        info!(
            "Bins: {n_bins} in {} ({written} contigs)",
            self.output_directory
        );

        Ok(held)
    }

    fn circular(&self, name: &str, header: &[u8], sequence: &[u8]) -> bool {
        self.tables.circles.contains(name) || circular::marked(header) || circular::closes(sequence)
    }

    fn leftover(&self, contig_length: usize, floor: usize, singles: &mut usize) -> Target {
        if contig_length >= self.settings.bars.min_bin_size {
            *singles += 1;
            return Target::Single(*singles);
        }
        if contig_length < floor {
            return Target::Small;
        }
        Target::Unbinned
    }
}

const PROGRESS_EVERY: usize = 4096;

#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub(super) enum Target {
    Bin(usize),
    Eukaryote(usize),
    Replicon(usize),
    Circular(usize),
    Single(usize),
    Unbinned,
    Small,
}

impl Target {
    pub(super) fn name(self) -> String {
        match self {
            Self::Bin(bin) => bin.to_string(),
            Self::Eukaryote(bin) => format!("{EUKARYOTE_PREFIX}{bin}"),
            Self::Replicon(at) => format!("{REPLICON_PREFIX}{at}"),
            Self::Circular(at) => format!("circular_{at}"),
            Self::Single(at) => format!("single_contig_{at}"),
            Self::Unbinned => UNBINNED.to_string(),
            Self::Small => format!("small_{UNBINNED}"),
        }
    }
}
