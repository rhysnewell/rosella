use std::collections::{BTreeMap, BTreeSet};

use anyhow::Result;
use log::debug;
use rayon::prelude::*;

use crate::clustering::clusterer::Partitioning;
use crate::embedding::weight::{Contigs, NEIGHBOURS, STEPS, centre, recall};
use crate::embedding::{Graph, features::row_slice, knn::KnnGraph};
use crate::learned::{Space, labelled::Labelled};
use crate::progress::{self, Bar, Stage};
use crate::quality::{Scorer, edge_spread};
use crate::recover::recover_engine::RecoverEngine;
use crate::refine::ranking::sorted;
use crate::refine::rung::FALLBACK_GENOME;

const SAMPLE: usize = 2_000;
// Each pass is a neighbour search and a partition, so a weight that wanders is cut off here.
const SETTLE_PASSES: usize = 8;
// Two plateaus with the same mean can differ in the last bits, which ran a pass twice.
const SAME_WEIGHT: f64 = 1e-9;

pub(super) struct Worth {
    pub kept: f64,
    pub spread: f64,
}

pub(super) struct Settled {
    pub graph: Graph,
    pub knn: KnnGraph,
    pub partitioning: Partitioning,
    pub weight: Option<f64>,
    pub parked: Vec<usize>,
    pub worth: Worth,
}

type Pass = (Graph, KnnGraph, Partitioning);

struct Opening {
    start: Option<f64>,
    partitioned: Pass,
}

struct Walk {
    weight: f64,
    spread: f64,
}

fn passing(label: String) -> Bar {
    let bar = progress::hold(progress::spinning(Stage::Partitioning));
    bar.set_message(label);
    bar
}

impl RecoverEngine {
    // Gold says short contigs in the partition move foreign long contigs home but add no bins,
    // and the short contigs a partition places are mostly foreign. So they wait for attach.
    pub(super) fn partitioned(&self, contigs: &[usize]) -> Result<Settled> {
        let long = self.long_end(contigs);
        let mut settled = self.weighted_partition(&contigs[..long])?;
        if long < contigs.len() {
            settled.worth.kept = self.pass_worth(&settled.partitioning, contigs);
            settled.parked = contigs[long..].to_vec();
        }
        Ok(settled)
    }

    fn opening(&self, contigs: &[usize]) -> Result<Opening> {
        let held = passing("sample pass".to_string());
        let (graph, knn) = self.embed(&self.features(), contigs);
        let first = self.partition_all(&graph, contigs)?;
        progress::done(&held);
        let start = self.derived_weight(&self.near_complete(&first, contigs))?;
        let partitioned = match start {
            None => (graph, knn, first),
            Some(weight) => self.pass(contigs, weight, 0)?,
        };
        Ok(Opening { start, partitioned })
    }

    // Given bins stand in for the first partition. Their near complete bins set the weight as a
    // first pass's would, and their short contigs wait for attach as a partition's do.
    pub(super) fn given_partition(
        &self,
        contigs: &[usize],
        bins: &BTreeMap<usize, Vec<usize>>,
    ) -> Result<Settled> {
        let long = self.long_end(contigs);
        let partitioning = self.given_long(bins, long);
        let (weight, graph, knn) = self.weighed(&partitioning, contigs)?;
        let given = || partitioning.cluster_map.values().map(Vec::as_slice);
        self.quality
            .learn_checkm(given(), self.settings.bars.quality());
        self.quality.learn_multicopy(given());
        Ok(Settled {
            graph,
            knn,
            worth: Worth {
                kept: self.pass_worth(&partitioning, contigs),
                spread: 0.0,
            },
            partitioning,
            weight,
            parked: contigs[long..].to_vec(),
        })
    }

    pub(super) fn weighed(
        &self,
        partitioning: &Partitioning,
        contigs: &[usize],
    ) -> Result<(Option<f64>, Graph, KnnGraph)> {
        let weight = self.derived_weight(&self.near_complete(partitioning, contigs))?;
        let mut settings = self.settings.features.clone();
        settings.distance.aggregate_weight = weight;
        let long = &contigs[..self.long_end(contigs)];
        let (graph, knn) = self.embed(&self.tables.features(&settings), long);
        Ok((weight, graph, knn))
    }

    // Rows run long first, so every row under `long` is a contig the partition would take.
    fn given_long(&self, bins: &BTreeMap<usize, Vec<usize>>, long: usize) -> Partitioning {
        Partitioning {
            cluster_map: bins
                .iter()
                .filter_map(|(bin, members)| {
                    let held = members
                        .iter()
                        .copied()
                        .filter(|contig| *contig < long)
                        .collect::<Vec<_>>();
                    (!held.is_empty()).then_some((*bin, held))
                })
                .collect(),
            outliers: Vec::new(),
            score: None,
            arm: self.settings.partition,
            seed: self.settings.seeds.partition,
        }
    }

    pub(super) fn given_labels(&self, contigs: &[usize]) -> Vec<Vec<usize>> {
        self.given.as_ref().map_or_else(Vec::new, |given| {
            let partitioning = self.given_long(&given.bins, contigs.len());
            self.near_complete(&partitioning, contigs)
        })
    }

    fn labelled(&self) -> Option<&Labelled> {
        match &self.tables.learned {
            Some(Space::Labelled(labelled)) => Some(labelled),
            _ => None,
        }
    }

    fn weighted_partition(&self, contigs: &[usize]) -> Result<Settled> {
        let Opening {
            start,
            mut partitioned,
        } = self.opening(contigs)?;
        let mut worth = Worth {
            kept: 0.0,
            spread: 0.0,
        };
        let weight = match start {
            None => None,
            Some(start) => {
                let settled = self.settle(contigs, start, &mut partitioned)?;
                worth.spread = settled.spread;
                Some(settled.weight)
            }
        };
        let (graph, knn, partitioning) = partitioned;
        // Learned after the weight is fixed, so the fit never reads it.
        self.learn_multicopy(&partitioning, contigs);
        Ok(Settled {
            graph,
            knn,
            partitioning,
            weight,
            parked: Vec::new(),
            worth,
        })
    }

    fn settle(&self, contigs: &[usize], start: f64, partitioned: &mut Pass) -> Result<Walk> {
        let mut near_complete = self.near_complete(&partitioned.2, contigs);
        let mut trace = vec![(start, self.pass_worth(&partitioned.2, contigs))];
        let mut best = trace[0];
        let mut seen = vec![partitioned.2.clone()];
        while trace.len() < SETTLE_PASSES
            && let Some(weight) = self.derived_weight(&near_complete)?
            && trace
                .iter()
                .all(|(used, _)| (used - weight).abs() > SAME_WEIGHT)
        {
            let next = self.pass(contigs, weight, trace.len())?;
            seen.push(next.2.clone());
            let worth = self.pass_worth(&next.2, contigs);
            let spread = self.spread(&partitioned.2, &next.2, contigs);
            // The one standard error rule (Breiman et al. 1984): an edge inside its own noise
            // keeps the pass already held.
            let clears = worth - best.1 > spread;
            debug!(
                "Coverage weight {weight:.3}: worth {:+.0}, spread {spread:.0}, {}",
                worth - best.1,
                if clears { "kept" } else { "refused" }
            );
            trace.push((weight, worth));
            // A refused pass says nothing about where the weight belongs, so the walk ends there.
            if !clears {
                break;
            }
            best = (weight, worth);
            near_complete = self.near_complete(&next.2, contigs);
            *partitioned = next;
        }
        debug!(
            "Coverage weight: kept {:.3} (weight and worth per pass {trace:.3?})",
            best.0
        );
        if seen.len() > 1 && std::env::var_os("ROSELLA_PASS_CONTEST").is_some() {
            partitioned.2 = self.contest_passes(&seen, contigs);
        }
        Ok(Walk {
            weight: best.0,
            spread: best.1 - trace.iter().map(|(_, worth)| *worth).fold(best.1, f64::min),
        })
    }

    fn pass(&self, contigs: &[usize], weight: f64, at: usize) -> Result<Pass> {
        let held = passing(format!("pass {}", at + 1));
        let mut settings = self.settings.features.clone();
        settings.distance.aggregate_weight = Some(weight);
        let (graph, knn) = self.embed(&self.tables.features(&settings), contigs);
        let partitioning = self.partition_all(&graph, contigs)?;
        progress::done(&held);
        Ok((graph, knn, partitioning))
    }

    // Squared so one whole genome outweighs its markers split across two bins. Long contigs only,
    // so a pass holding short contigs is judged on the markers the long-only pass had.
    pub(super) fn pass_worth(&self, partitioning: &Partitioning, contigs: &[usize]) -> f64 {
        partitioning
            .cluster_map
            .par_iter()
            .map(|(_, members)| {
                let long = self.long(members.iter().map(|at| contigs[*at]));
                let worth = self.quality.score(&long).score(self.settings.bars.worth);
                worth.max(0.0).powi(2)
            })
            .collect::<Vec<_>>()
            .into_iter()
            .sum()
    }

    // A bin both passes hold reads the same under any resample of the markers, so it cancels.
    fn spread(&self, kept: &Partitioning, next: &Partitioning, contigs: &[usize]) -> f64 {
        let bins = |partitioning: &Partitioning| {
            partitioning
                .cluster_map
                .values()
                .map(|members| self.long(members.iter().map(|at| contigs[*at])))
                .filter(|members| !members.is_empty())
                .collect::<BTreeSet<_>>()
        };
        let (before, after) = (bins(kept), bins(next));
        let points = |bin: &Vec<usize>| self.quality.points(bin, self.settings.bars.worth);
        edge_spread(
            after
                .difference(&before)
                .map(|bin| (1.0, points(bin)))
                .chain(before.difference(&after).map(|bin| (-1.0, points(bin)))),
        )
    }

    fn long(&self, members: impl Iterator<Item = usize>) -> Vec<usize> {
        let lengths = &self.tables.coverage.contig_lengths;
        sorted(members.filter(|contig| lengths[*contig] >= self.settings.lengths.cutoff))
    }

    pub(super) fn learn_multicopy(&self, partitioning: &Partitioning, contigs: &[usize]) {
        let bins = partitioning
            .cluster_map
            .values()
            .map(|members| members.iter().map(|at| contigs[*at]).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        self.quality.learn_multicopy(bins.iter().map(Vec::as_slice));
    }

    fn near_complete(&self, partitioning: &Partitioning, contigs: &[usize]) -> Vec<Vec<usize>> {
        let bars = self.settings.bars.quality();
        partitioning
            .cluster_map
            .par_iter()
            .map(|(_, members)| sorted(members.iter().map(|at| contigs[*at])))
            .filter(|members| self.quality.score(members).clears(bars))
            .collect()
    }

    // Bases over completeness is a genome's size, so the smallest near complete bin measures the
    // smallest genome the run can vouch for.
    pub(super) fn smallest_genome(&self, partitioning: &Partitioning, contigs: &[usize]) -> usize {
        let lengths = self.lengths();
        let smallest = self
            .near_complete(partitioning, contigs)
            .iter()
            .map(|members| {
                let bp = members.iter().map(|contig| lengths[*contig]).sum::<usize>();
                100.0 * bp as f64 / self.quality.score(members).completeness
            })
            .min_by(f64::total_cmp);
        smallest.map_or(FALLBACK_GENOME, |genome| genome as usize)
    }

    // Every contig split in two must leave halves long enough to be binned in their own right.
    fn judged(&self, near_complete: &[Vec<usize>]) -> Option<Vec<usize>> {
        let lengths = &self.tables.coverage.contig_lengths;
        let pool = sorted(
            near_complete
                .iter()
                .flatten()
                .copied()
                .filter(|contig| lengths[*contig] >= 2 * self.settings.lengths.cutoff)
                .filter(|contig| self.labelled().is_none_or(|learned| learned.holds(*contig))),
        );
        let order =
            crate::seeds::sample_positions(pool.len(), pool.len(), self.settings.seeds.seed);
        if order.len() <= NEIGHBOURS {
            debug!(
                "Coverage weight: unchanged ({} near complete bins, {} contigs)",
                near_complete.len(),
                order.len()
            );
            return None;
        }
        Some(order.iter().map(|at| pool[*at]).collect())
    }

    // One draw of SAMPLE contigs moved marine's weight from 0.575 to 0.68 with the seed. The curve
    // is averaged over the whole pool in blocks of the size the halves test is built for.
    fn derived_weight(&self, near_complete: &[Vec<usize>]) -> Result<Option<f64>> {
        let _timer = crate::timing::scope("weight");
        let Some(pool) = self.judged(near_complete) else {
            return Ok(None);
        };
        let blocks = pool.len().div_ceil(SAMPLE);
        let mut total = [0.0; STEPS];
        let mut counted = 0;
        for block in pool.chunks(pool.len().div_ceil(blocks)) {
            if let Some(curve) = self.recall_curve(block)? {
                for (sum, value) in total.iter_mut().zip(curve) {
                    *sum += value * block.len() as f64;
                }
                counted += block.len();
            }
        }
        let weight = (counted > 0).then(|| centre(&total.map(|sum| sum / counted as f64)));
        debug!(
            "Coverage weight: {} ({} contigs in {blocks} blocks, {} near complete bins)",
            weight.map_or("unchanged".to_string(), |weight| format!("{weight:.3}")),
            pool.len(),
            near_complete.len()
        );
        Ok(weight)
    }

    fn recall_curve(&self, chosen: &[usize]) -> Result<Option<[f64; STEPS]>> {
        let names = chosen
            .iter()
            .map(|contig| self.tables.coverage.contig_names[*contig].as_str())
            .collect::<Vec<_>>();
        let mut halves = self.tables.halves.composition(&names)?;
        let coverage = rows(&self.tables.coverage.table, chosen.iter().copied());
        let lengths = chosen
            .iter()
            .map(|contig| self.lengths()[*contig])
            .collect::<Vec<_>>();
        if let Some(space) = &self.tables.learned {
            halves = space.halves(halves, &coverage, &lengths, self.settings.seeds.seed);
        }
        let [first, second] = halves;
        let contigs = Contigs {
            coverage,
            whole: rows(&self.tables.tnf.kmer_table, chosen.iter().copied()),
            first: rows(&first, 0..chosen.len()),
            second: rows(&second, 0..chosen.len()),
        };
        Ok(recall(
            &contigs,
            &self.settings.features.distance,
            self.settings.seeds.seed,
        ))
    }
}

fn rows(table: &ndarray::Array2<f64>, at: impl Iterator<Item = usize>) -> Vec<&[f64]> {
    at.map(|row| row_slice(table, row)).collect()
}
