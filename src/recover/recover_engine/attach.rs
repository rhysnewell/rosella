use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::ops::Range;

use anyhow::Result;
use log::debug;
use ndarray::Axis;
use rayon::prelude::*;

use crate::embedding::knn::{KnnGraph, nearest_in};
use crate::progress::{self, Stage};
use crate::quality::Scorer;
use crate::recover::floor_walk::bands;
use crate::recover::offered_votes::{Ballot, ballots};
use crate::recover::recover_engine::RecoverEngine;
use crate::recover::recover_engine::weight::Worth;
use crate::refine::owners::{credit, owners};

mod plan;

// One complete, clean bin in squared marker worth. A smaller gain cannot add a bin.
const ONE_GENOME: f64 = 100.0 * 100.0;

struct Evidence {
    contig: usize,
    bin: usize,
    repeats: Option<bool>,
    in_place: bool,
    complete: f64,
}

pub(super) type Proposal = (usize, Option<(usize, f32)>);

pub(super) struct Searched {
    pub(super) proposals: Vec<Proposal>,
    pub(super) bar: f32,
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
    ballots: &'a [Ballot],
}

impl Down<'_> {
    fn open(&self, span: Range<usize>) -> Vec<usize> {
        span.filter(|at| self.ballots.get(*at).is_none_or(Ballot::can_join))
            .collect()
    }
}

impl RecoverEngine {
    // Rows under the cutoff are read only as the walk reaches their band, since a walk that stops
    // early would otherwise hold millions of short contigs it never looks at.
    pub(super) fn load(&mut self, through: usize) -> Result<()> {
        let loaded = self.tables.tnf.kmer_table.nrows();
        if through <= loaded {
            return Ok(());
        }
        let (composition, coverage) = self.tables.read_rows(loaded..through)?;
        let composition = match &self.tables.learned {
            Some(learned) => {
                let _timer = crate::timing::scope("learned_embed");
                learned.embed(composition.view(), coverage.view())
            }
            None => composition,
        };
        self.tables
            .coverage
            .table
            .append(Axis(0), coverage.view())?;
        self.tables
            .tnf
            .kmer_table
            .append(Axis(0), composition.view())?;
        let names = &self.tables.coverage.contig_names[loaded..through];
        self.tables.tnf.contig_names.extend_from_slice(names);
        Ok(())
    }

    pub(super) fn attach(
        &mut self,
        bins: &mut BTreeMap<usize, Vec<usize>>,
        unbinned: &mut HashSet<usize>,
        mut order: Vec<usize>,
        worth: Worth,
        long_graph: &KnnGraph,
    ) -> Result<usize> {
        let given = self.settings.lengths.attach_floor;
        if order.is_empty() {
            return Ok(given.unwrap_or(0));
        }
        let progress = progress::hold(progress::spinning(Stage::AttachingContigs));
        progress.set_message(format!("{} short contigs", order.len()));
        let lengths = &self.tables.coverage.contig_lengths;
        let bin_of = owners(bins.iter().map(|(bin, members)| (*bin, members)));
        if given.is_none() {
            order.sort_by_key(|contig| (std::cmp::Reverse(lengths[*contig]), *contig));
        }
        let spans: Vec<Range<usize>> = match given {
            Some(_) => std::iter::once(0..order.len()).collect(),
            None => bands(
                &order.iter().map(|at| lengths[*at]).collect::<Vec<_>>(),
                long_graph.indices.nrows(),
            ),
        };
        let ballots = ballots(
            &order,
            &spans,
            lengths,
            self.settings.lengths.cutoff,
            &self.candidates,
            &bin_of,
        );
        let down = Down {
            order: &order,
            long_graph,
            bin_of: &bin_of,
            ballots: &ballots,
        };
        let bars = match given {
            Some(_) => vec![0.5],
            None => self.plan(&down, &spans, bins, worth)?,
        };
        let walked = match bars.len() {
            0 => Walked {
                floor: self.settings.lengths.cutoff,
                ..Walked::default()
            },
            searched => {
                let among = {
                    let _timer = crate::timing::scope("attach");
                    self.among(&order[..spans[searched - 1].end])
                };
                self.walk_down(&down, bins, &spans, &bars, &among)
            }
        };
        let refused = self.refused(&walked.evidence);
        if let Some(path) = &self.reports.attach_report
            && let Err(error) = self.write_attach_report(path, bins, &walked.searched, &refused)
        {
            log::warn!("Attach report: {error} ({})", path.display());
        }
        let mut taken = 0;
        for (contig, bin) in walked
            .joins
            .iter()
            .filter(|(_, bin)| !refused.contains(bin))
        {
            unbinned.remove(contig);
            bins.entry(*bin).or_default().push(*contig);
            taken += 1;
        }
        for contigs in bins.values_mut() {
            contigs.sort_unstable();
        }
        debug!(
            "Attach: {} of {} parked contigs from {} bp placed, {} bins refused, {taken} joined",
            walked.joins.len(),
            order.len(),
            walked.floor,
            refused.len(),
        );
        progress::done(&progress);
        Ok(given.unwrap_or(walked.floor))
    }

    fn walk_down(
        &self,
        down: &Down,
        bins: &BTreeMap<usize, Vec<usize>>,
        spans: &[Range<usize>],
        bars: &[f32],
        among: &KnnGraph,
    ) -> Walked {
        let lengths = &self.tables.coverage.contig_lengths;
        let mut walked = Walked {
            floor: self.settings.lengths.cutoff,
            ..Walked::default()
        };
        for (span, bar) in spans.iter().zip(bars) {
            let band = &down.order[span.clone()];
            let proposals = self.propose(&down.open(span.clone()), down, among);
            let joined = joining(&proposals, *bar);
            walked.evidence.extend(self.marker_evidence(bins, &joined));
            walked.joins.extend(joined);
            walked.searched.push(Searched {
                proposals,
                bar: *bar,
            });
            walked.floor = band
                .iter()
                .map(|contig| lengths[*contig])
                .fold(walked.floor, usize::min);
        }
        walked
    }

    fn propose(&self, open: &[usize], down: &Down, among: &KnnGraph) -> Vec<Proposal> {
        let first = down.long_graph.indices.nrows();
        let band = open.iter().map(|at| down.order[*at]).collect::<Vec<_>>();
        let nearest = self.nearest_long(down.long_graph, &band);
        let _timer = crate::timing::scope("attach");
        let own = KnnGraph {
            indices: among.indices.select(Axis(0), open),
            dists: among.dists.select(Axis(0), open),
        };
        let knn = merge(&nearest, &own, band.len(), |at| down.order[at]);
        band.into_iter()
            .zip(best_bins(&knn, first, down.bin_of, |row| {
                down.ballots.get(open[row])
            }))
            .collect()
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
            self.settings.knn.candidates,
            self.settings.seeds.knn,
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
        let foreign = match repeats > 0.0 {
            true => (repeats / complete).min(1.0),
            false => 0.0,
        };
        let mut trade = HashMap::<usize, f64>::new();
        for (seen, repeats) in whole {
            *trade.entry(seen.bin).or_default() += match repeats {
                true => -self.settings.bars.worth,
                false => 1.0 - foreign * (1.0 - seen.complete),
            };
        }
        debug!(
            "Attach markers: {repeats} repeats against {complete:.1} if all foreign (foreign share {foreign:.2})"
        );
        trade
            .into_iter()
            .filter(|(_, gain)| *gain <= 0.0)
            .map(|(bin, _)| bin)
            .collect()
    }

    fn marker_evidence(
        &self,
        bins: &BTreeMap<usize, Vec<usize>>,
        joins: &[(usize, usize)],
    ) -> Vec<Evidence> {
        let marked = joins
            .par_iter()
            .filter(|(contig, _)| self.quality.hit_count(&[*contig]) > 0)
            .collect::<Vec<_>>();
        let complete = marked
            .iter()
            .map(|(_, bin)| *bin)
            .collect::<BTreeSet<_>>()
            .into_par_iter()
            .map(|bin| (bin, self.quality.score(&bins[&bin]).completeness / 100.0))
            .collect::<HashMap<_, _>>();
        marked
            .into_par_iter()
            .filter_map(|(contig, bin)| {
                let mut with = bins[bin].clone();
                with.insert(with.partition_point(|at| at < contig), *contig);
                Some(Evidence {
                    contig: *contig,
                    bin: *bin,
                    repeats: self.quality.repeats(&with, *contig),
                    in_place: self.quality.repeats_in_place(&with, *contig)?,
                    complete: complete[bin],
                })
            })
            .collect()
    }

    fn among(&self, band: &[usize]) -> KnnGraph {
        self.features().knn_of(
            band,
            self.settings.knn.neighbours,
            self.settings.knn.candidates,
            self.settings.seeds.knn,
            crate::embedding::KNN_ATTACH,
        )
    }
}

// The long half of each neighbourhood is already known from the nearest search, so only the
// parked contigs are searched among themselves and the cost follows their count.
pub(super) fn merge(
    long: &KnnGraph,
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
        let mut both = long
            .indices
            .row(row)
            .iter()
            .zip(long.dists.row(row))
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

fn joining(proposals: &[Proposal], bar: f32) -> Vec<(usize, usize)> {
    proposals
        .iter()
        .filter_map(|(contig, best)| {
            let (bin, share) = (*best)?;
            (share > bar).then_some((*contig, bin))
        })
        .collect()
}

// Only long contigs vouch for a bin, over a neighbourhood that counts every contig. Letting
// short ones vouch fills sink bins through chains of short contigs. An offered set that bins all
// of a band weighs as much as the whole neighbourhood.
pub(super) fn best_bins<'a>(
    knn: &KnnGraph,
    first: usize,
    bin_of: &HashMap<usize, usize>,
    ballot: impl Fn(usize) -> Option<&'a Ballot> + Sync,
) -> Vec<Option<(usize, f32)>> {
    let (sigmas, rhos) = crate::embedding::fuzzy::scales(knn.dists.view(), knn.indices.ncols());
    (0..knn.indices.nrows())
        .into_par_iter()
        .map(|row| {
            let mut mass = Vec::<(usize, f32)>::new();
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
                    credit(&mut mass, *bin, weight);
                }
            }
            let scale = total.max(f32::MIN_POSITIVE);
            let say = ballot(row).map_or(0.0, |ballot| {
                for (bin, weight) in &ballot.bins {
                    credit(&mut mass, *bin, scale * weight);
                }
                ballot.say
            });
            mass.into_iter()
                .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)))
                .map(|(bin, held)| (bin, held / (scale * (1.0 + say))))
        })
        .collect()
}
