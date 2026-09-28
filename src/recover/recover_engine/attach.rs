use std::collections::{HashMap, HashSet};
use std::ops::Range;

use anyhow::Result;
use log::info;
use rayon::prelude::*;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::{
    Graph,
    knn::{KnnGraph, nearest_in},
};
use crate::quality::Scorer;
use crate::recover::floor_walk::bands;
use crate::recover::recover_engine::RecoverEngine;
use plan::Plan;

mod plan;

// One complete, clean bin in squared marker worth. A smaller gain cannot add a bin.
const ONE_GENOME: f64 = 100.0 * 100.0;

struct Evidence {
    bin: usize,
    repeats: Option<bool>,
    in_place: bool,
    complete: f64,
}

pub(super) type Proposal = (usize, Option<(usize, f32)>);

pub(super) struct Searched {
    pub(super) proposals: Vec<Proposal>,
    pub(super) chances: Option<Vec<f64>>,
    pub(super) taken: bool,
}

#[derive(Default)]
struct Walked {
    searched: Vec<Searched>,
    joins: Vec<(usize, usize)>,
    evidence: Vec<Evidence>,
    floor: usize,
}

struct Down<'a> {
    order: &'a [usize],
    long_graph: &'a KnnGraph,
    bin_of: &'a HashMap<usize, usize>,
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
        self.kept = self.pass_worth(&settled, contigs);
        self.parked = contigs[long..].to_vec();
        Ok((graph, knn, settled))
    }

    pub(super) fn attach(
        &mut self,
        bins: &mut HashMap<usize, HashSet<usize>>,
        unbinned: &mut HashSet<usize>,
        long_graph: &KnnGraph,
    ) -> Result<()> {
        let mut order = std::mem::take(&mut self.parked);
        if order.is_empty() {
            return Ok(());
        }
        let lengths = &self.coverage_table.contig_lengths;
        let bin_of = bins
            .iter()
            .flat_map(|(bin, members)| members.iter().map(move |contig| (*contig, *bin)))
            .collect::<HashMap<_, _>>();
        let (spans, plan) = match self.attach_given {
            true => (std::iter::once(0..order.len()).collect(), Plan::given()),
            false => {
                order.sort_by_key(|contig| (std::cmp::Reverse(lengths[*contig]), *contig));
                let spans = bands(
                    &order.iter().map(|at| lengths[*at]).collect::<Vec<_>>(),
                    long_graph.indices.nrows(),
                );
                let plan = self.plan(&order, &spans, long_graph, &bin_of, bins)?;
                (spans, plan)
            }
        };
        let down = Down {
            order: &order,
            long_graph,
            bin_of: &bin_of,
        };
        let walked = match plan.searched {
            0 => Walked {
                floor: self.cutoff,
                ..Walked::default()
            },
            searched => {
                let among = {
                    let _timer = crate::timing::scope("attach");
                    self.among(&order[..spans[searched - 1].end])
                };
                self.walk_down(&down, bins, &spans, &plan.taken, &among)
            }
        };
        if !self.attach_given {
            self.min_contig_size = walked.floor;
        }
        let refused = self.refused(&walked.evidence);
        if let Some(path) = &self.attach_report
            && let Err(error) = self.write_attach_report(path, bins, &walked.searched, &refused)
        {
            log::warn!("Could not write {}: {error}", path.display());
        }
        let mut taken = 0;
        for (contig, bin) in walked
            .joins
            .iter()
            .filter(|(_, bin)| !refused.contains(bin))
        {
            unbinned.remove(contig);
            bins.entry(*bin).or_default().insert(*contig);
            taken += 1;
        }
        info!(
            "{} of {} parked short contigs from {} bp sit in one bin's neighbourhood. {} bins \
             refuse theirs on marker evidence and {taken} join.",
            walked.joins.len(),
            order.len(),
            walked.floor,
            refused.len(),
        );
        Ok(())
    }

    fn walk_down(
        &self,
        down: &Down,
        bins: &HashMap<usize, HashSet<usize>>,
        spans: &[Range<usize>],
        taken: &[bool],
        among: &KnnGraph,
    ) -> Walked {
        let lengths = &self.coverage_table.contig_lengths;
        let mut walked = Walked {
            floor: self.cutoff,
            ..Walked::default()
        };
        for (span, taken) in spans.iter().zip(taken) {
            let band = &down.order[span.clone()];
            let (proposals, chances) = self.propose(band, down, among, down.order, span.clone());
            let joined = joining(&proposals, chances.as_deref());
            walked.evidence.extend(self.marker_evidence(bins, &joined));
            walked.joins.extend(joined);
            walked.searched.push(Searched {
                proposals,
                chances,
                taken: *taken,
            });
            walked.floor = lengths[band[band.len() - 1]];
        }
        walked
    }

    fn propose(
        &self,
        band: &[usize],
        down: &Down,
        among: &KnnGraph,
        among_contigs: &[usize],
        rows: Range<usize>,
    ) -> (Vec<Proposal>, Option<Vec<f64>>) {
        let first = down.long_graph.indices.nrows();
        let nearest = self.nearest_long(down.long_graph, band);
        let _timer = crate::timing::scope("attach");
        let own = KnnGraph {
            indices: among
                .indices
                .slice(ndarray::s![rows.clone(), ..])
                .to_owned(),
            dists: among.dists.slice(ndarray::s![rows, ..]).to_owned(),
        };
        let knn = merge(
            &nearest,
            |row| row,
            &own,
            band.len(),
            |at| among_contigs[at],
        );
        let proposals = band
            .iter()
            .copied()
            .zip(best_bins(&knn, first, down.bin_of))
            .collect::<Vec<_>>();
        let chances = self
            .attach_calibrate
            .then(|| {
                self.home_chances(
                    down.bin_of,
                    &proposals,
                    down.long_graph,
                    among,
                    among_contigs,
                )
            })
            .flatten();
        (proposals, chances)
    }

    fn nearest_long(&self, knn: &KnnGraph, band: &[usize]) -> KnnGraph {
        let _timer = crate::timing::scope("nearest");
        let first = knn.indices.nrows();
        let indices = (0..first).chain(band.iter().copied()).collect::<Vec<_>>();
        let prepared = self.features().prepared(&indices);
        nearest_in(
            knn,
            band.len(),
            knn.indices.ncols(),
            self.knn_candidates,
            self.seeds.knn,
            prepared.shifted(first),
        )
    }

    // A foreign contig repeats a marker as often as its bin is complete and an own one almost never,
    // so each bin weighs its own fills against its repeats at the run's contamination weight.
    fn refused(&self, evidence: &[Evidence]) -> HashSet<usize> {
        let whole = evidence
            .iter()
            .filter_map(|seen| Some((seen, seen.repeats?)))
            .collect::<Vec<_>>();
        let repeats = whole.iter().filter(|(_, repeats)| *repeats).count() as f64;
        let complete = whole.iter().map(|(seen, _)| seen.complete).sum::<f64>();
        let foreign = (repeats / complete).min(1.0);
        let mut trade = HashMap::<usize, f64>::new();
        for (seen, repeats) in whole {
            *trade.entry(seen.bin).or_default() += match repeats {
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
                    repeats: self.quality.repeats(&with, *contig),
                    in_place: self.quality.repeats_in_place(&with, *contig)?,
                    complete: self.quality.score(rest).completeness / 100.0,
                })
            })
            .collect()
    }

    fn among(&self, band: &[usize]) -> KnnGraph {
        self.features().knn_of(
            band,
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
    among_contig: impl Fn(usize) -> usize,
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
                    .map(|(at, distance)| (among_contig(*at as usize) as u32, *distance)),
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

fn joining(proposals: &[Proposal], chances: Option<&[f64]>) -> Vec<(usize, usize)> {
    proposals
        .iter()
        .enumerate()
        .filter_map(|(at, (contig, best))| {
            let (bin, share) = (*best)?;
            let keep = match chances {
                Some(chances) => chances[at] > 0.5,
                None => share > 0.5,
            };
            keep.then_some((*contig, bin))
        })
        .collect()
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
