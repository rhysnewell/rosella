use std::collections::{HashMap, HashSet};

use anyhow::Result;
use log::info;
use rayon::prelude::*;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::{
    Graph,
    knn::{KnnGraph, nearest_in},
};
use crate::quality::Scorer;
use crate::recover::recover_engine::RecoverEngine;

// One complete, clean bin in squared marker worth. A smaller gain cannot add a bin.
const ONE_GENOME: f64 = 100.0 * 100.0;

pub(super) struct Nearest {
    first: usize,
    knn: KnnGraph,
}

struct Evidence {
    bin: usize,
    repeats: bool,
    complete: f64,
}

impl RecoverEngine {
    // Gold says short contigs in the partition move foreign long contigs home but add no bins,
    // and the short contigs a partition places are mostly foreign. So they wait for attach.
    pub(super) fn partitioned(
        &mut self,
        contigs: &[usize],
    ) -> Result<(Graph, KnnGraph, Partitioning)> {
        let lengths = &self.coverage_table.contig_lengths;
        let long = contigs.partition_point(|contig| lengths[*contig] >= self.cutoff);
        if long == contigs.len() {
            return self.weighted_partition(contigs);
        }
        let (graph, knn, settled) = self.weighted_partition(&contigs[..long])?;
        let kept = self.pass_worth(&settled, contigs);
        let bar = self.worth_spread.max(ONE_GENOME);
        let share = self.quality.hit_count(&contigs[long..]) as f64
            / self.quality.hit_count(&contigs[..long]).max(1) as f64;
        let reach = kept * ((1.0 + share).powi(2) - 1.0);
        info!(
            "{} shorter contigs carry {share:.3} of the long contigs' markers, worth at most \
             {reach:.0} against a bar of {bar:.0}.",
            contigs.len() - long
        );
        self.parked = contigs[long..].to_vec();
        if self.attach_given || reach > bar {
            self.nearest = Some(self.nearest_long(&knn, contigs, long));
        }
        Ok((graph, knn, settled))
    }

    fn nearest_long(&self, knn: &KnnGraph, contigs: &[usize], long: usize) -> Nearest {
        let _timer = crate::timing::scope("nearest");
        let prepared = self.features().prepared(contigs);
        Nearest {
            first: contigs[long],
            knn: nearest_in(
                knn,
                contigs.len() - long,
                knn.indices.ncols(),
                self.knn_candidates,
                self.seeds.knn,
                |short, base| prepared.distance(long + short, base),
            ),
        }
    }

    pub(super) fn attach(
        &self,
        bins: &mut HashMap<usize, HashSet<usize>>,
        unbinned: &mut HashSet<usize>,
        parked: &[usize],
        nearest: &Nearest,
        long_graph: &KnnGraph,
    ) {
        let _timer = crate::timing::scope("attach");
        let among = self.among(parked);
        let knn = merge(
            &nearest.knn,
            |row| parked[row] - nearest.first,
            &among,
            parked.len(),
            nearest.first,
        );
        let bin_of = bins
            .iter()
            .flat_map(|(bin, members)| members.iter().map(move |contig| (*contig, *bin)))
            .collect::<HashMap<_, _>>();
        let proposals = parked
            .iter()
            .copied()
            .zip(best_bins(&knn, nearest.first, &bin_of))
            .collect::<Vec<_>>();
        let chances = self
            .attach_calibrate
            .then(|| self.home_chances(&bin_of, &proposals, long_graph, &among, nearest.first))
            .flatten();
        let joins = proposals
            .iter()
            .enumerate()
            .filter_map(|(at, (contig, best))| {
                let (bin, share) = (*best)?;
                let keep = match &chances {
                    Some(chances) => chances[at] > 0.5,
                    None => share > 0.5,
                };
                keep.then_some((*contig, bin))
            })
            .collect::<Vec<_>>();
        let refused = self.refused(bins, &joins);
        if let Some(path) = &self.attach_report
            && let Err(error) =
                self.write_attach_report(path, bins, &proposals, chances.as_deref(), &refused)
        {
            log::warn!("Could not write {}: {error}", path.display());
        }
        let mut taken = 0;
        for (contig, bin) in joins.iter().filter(|(_, bin)| !refused.contains(bin)) {
            unbinned.remove(contig);
            bins.entry(*bin).or_default().insert(*contig);
            taken += 1;
        }
        info!(
            "{} of {} parked short contigs sit in one bin's neighbourhood. {} bins refuse theirs \
             on marker evidence and {taken} join.",
            joins.len(),
            parked.len(),
            refused.len(),
        );
    }

    // A foreign contig repeats a marker as often as its bin is complete and an own one almost never,
    // so each bin weighs its own fills against its repeats at the run's contamination weight.
    fn refused(
        &self,
        bins: &HashMap<usize, HashSet<usize>>,
        joins: &[(usize, usize)],
    ) -> HashSet<usize> {
        let evidence = self.marker_evidence(bins, joins);
        let repeats = evidence.iter().filter(|seen| seen.repeats).count() as f64;
        let complete = evidence.iter().map(|seen| seen.complete).sum::<f64>();
        let foreign = (repeats / complete).min(1.0);
        let mut trade = HashMap::<usize, f64>::new();
        for seen in &evidence {
            *trade.entry(seen.bin).or_default() += match seen.repeats {
                true => -self.worth,
                false => 1.0 - foreign * (1.0 - seen.complete),
            };
        }
        info!(
            "Short contigs' markers repeat {repeats} times against {complete:.1} if all were \
             foreign, a foreign share of {foreign:.2}."
        );
        trade
            .into_iter()
            .filter(|(_, gain)| *gain <= 0.0)
            .map(|(bin, _)| bin)
            .collect()
    }

    fn write_attach_report(
        &self,
        path: &std::path::Path,
        bins: &HashMap<usize, HashSet<usize>>,
        proposals: &[(usize, Option<(usize, f32)>)],
        chances: Option<&[f64]>,
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
        let rows = proposals
            .par_iter()
            .enumerate()
            .map(|(at, (contig, best))| {
                let chance =
                    chances.map_or("NA".to_string(), |chances| format!("{:.3}", chances[at]));
                let Some((bin, share)) = best else {
                    return format!(
                        "{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t{chance}\tNA",
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
                    "{}\t{}\t{bin}\t{}\t{share:.3}\t{}\t{}\t{:.3}\t{}\t{chance}\t{}",
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
             repeats_place"
        )?;
        for row in rows {
            writeln!(sink, "{row}")?;
        }
        Ok(())
    }

    fn marker_evidence(
        &self,
        bins: &HashMap<usize, HashSet<usize>>,
        joins: &[(usize, usize)],
    ) -> Vec<Evidence> {
        let members = bins
            .iter()
            .map(|(bin, contigs)| {
                let mut contigs = contigs.iter().copied().collect::<Vec<_>>();
                contigs.sort_unstable();
                (*bin, contigs)
            })
            .collect::<HashMap<_, _>>();
        joins
            .par_iter()
            .filter(|(contig, _)| self.quality.hit_count(&[*contig]) > 0)
            .filter_map(|(contig, bin)| {
                let rest = &members[bin];
                let mut with = rest.clone();
                with.insert(with.partition_point(|at| at < contig), *contig);
                Some(Evidence {
                    bin: *bin,
                    repeats: self.quality.repeats(&with, *contig)?,
                    complete: self.quality.score(rest).completeness / 100.0,
                })
            })
            .collect()
    }

    fn among(&self, parked: &[usize]) -> KnnGraph {
        self.features().knn_of(
            parked,
            self.n_neighbours,
            self.seeds,
            self.knn_candidates,
            crate::embedding::KNN_ATTACH,
        )
    }
}

// The long half of each neighbourhood is already known from the nearest search, so only the
// parked contigs are searched among themselves and the cost follows their count.
pub(super) fn merge(
    long: &KnnGraph,
    long_row: impl Fn(usize) -> usize,
    among: &KnnGraph,
    rows: usize,
    first: usize,
) -> KnnGraph {
    let width = long.indices.ncols();
    let mut merged = KnnGraph {
        indices: ndarray::Array2::from_elem((rows, width), u32::MAX),
        dists: ndarray::Array2::from_elem((rows, width), f32::INFINITY),
    };
    for row in 0..rows {
        let own = long_row(row);
        let mut both = long
            .indices
            .row(own)
            .iter()
            .zip(long.dists.row(own))
            .map(|(at, distance)| (*at, *distance))
            .chain(
                among
                    .indices
                    .row(row)
                    .iter()
                    .zip(among.dists.row(row))
                    .filter(|(at, _)| **at != u32::MAX)
                    .map(|(at, distance)| ((first + *at as usize) as u32, *distance)),
            )
            .filter(|(at, distance)| *at != u32::MAX && distance.is_finite())
            .collect::<Vec<_>>();
        both.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        for (slot, (at, distance)) in both.into_iter().take(width).enumerate() {
            merged.indices[[row, slot]] = at;
            merged.dists[[row, slot]] = distance;
        }
    }
    merged
}

// Only long contigs vouch for a bin, over a neighbourhood that counts every contig. Letting
// short ones vouch fills sink bins through chains of short contigs.
pub(super) fn best_bins(
    knn: &KnnGraph,
    first: usize,
    bin_of: &HashMap<usize, usize>,
) -> Vec<Option<(usize, f32)>> {
    let (sigmas, rhos) = crate::embedding::fuzzy::scales(knn.dists.view(), knn.indices.ncols());
    (0..knn.indices.nrows())
        .into_par_iter()
        .map(|row| {
            let mut mass = HashMap::<usize, f32>::new();
            let mut total = 0.0;
            for (neighbour, distance) in knn.indices.row(row).iter().zip(knn.dists.row(row)) {
                if *neighbour == u32::MAX {
                    break;
                }
                let weight = crate::embedding::fuzzy::membership(*distance, rhos[row], sigmas[row]);
                total += weight;
                let neighbour = *neighbour as usize;
                if neighbour < first
                    && let Some(bin) = bin_of.get(&neighbour)
                {
                    *mass.entry(*bin).or_default() += weight;
                }
            }
            mass.into_iter()
                .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)))
                .map(|(bin, held)| (bin, held / total.max(f32::MIN_POSITIVE)))
        })
        .collect()
}
