use std::collections::{BTreeMap, HashMap, HashSet};

use anyhow::Result;
use log::warn;

use crate::clustering::clusterer::{HDBSCANResult, placed_once};
use crate::embedding::features::ContigFeatures;
use crate::quality::ContigQuality;
use crate::refine::select::{Selection, ranked, remaining};

pub const DISSOLVE_SCOPE_NAMES: [&str; 2] = ["fused", "all"];

const MIN_NEIGHBOURS: usize = 2;

pub const DEFAULT_COMPLETENESS: f64 = 90.0;
pub const DEFAULT_CONTAMINATION: f64 = 5.0;

/// Completeness walks down and contamination up together, so the loop takes the genomes it is
/// sure of first and only then the ones it is not.
const QUALITY_LADDER: [(f64, f64); 5] = [
    (1.0, 1.0),
    (0.89, 1.0),
    (0.78, 2.0),
    (0.67, 2.0),
    (0.56, 3.0),
];

/// Floor as a share of the gap between the bin floor and genome scale, and the duplication bar
/// as a multiple of its setting. Rung zero is the fixed bar the single pass always used.
const LADDER: [(f64, f64); 5] = [(1.0, 1.0), (0.75, 1.0), (0.5, 1.0), (0.5, 2.0), (0.5, 3.0)];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DissolveScope {
    Fused,
    All,
}

impl DissolveScope {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "fused" => Some(Self::Fused),
            "all" => Some(Self::All),
            _ => None,
        }
    }
}

/// A round searches the pool that the rounds before it left, so the graph is rebuilt over
/// material that shrinks, and k walks down with it.
#[derive(Debug, Clone, Copy)]
pub struct RoundParams {
    pub n_neighbours: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct DissolveSettings {
    pub min_bin_size: usize,
    pub genome_floor: Option<usize>,
    pub duplication_bar: f64,
    pub min_contigs: usize,
    pub scope: DissolveScope,
    pub rounds: usize,
    pub ladder: bool,
    pub n_neighbours: usize,
    pub completeness: f64,
    pub contamination: f64,
    pub improve: bool,
    pub select: Selection,
}

/// What the pool took, what it refused and where the refusals went, in contigs and bases.
/// The bin count alone cannot say whether an arm has headroom left.
#[derive(Debug, Default, Clone, Copy)]
pub struct DissolveLedger {
    pub dissolved_small: usize,
    pub dissolved_duplicated: usize,
    pub dissolved_clean: usize,
    pub dissolved_bp: usize,
    pub pool_contigs: usize,
    pub pool_bp: usize,
    pub rounds: usize,
    pub rung: usize,
    pub proposed: usize,
    pub refused_small: usize,
    pub refused_duplicated: usize,
    pub refused_worse: usize,
    pub noise: usize,
    pub promoted: usize,
    pub adopted_contigs: usize,
    pub adopted_bp: usize,
    pub returned_contigs: usize,
    pub returned_bp: usize,
    pub emptied: usize,
    pub left_contigs: usize,
    pub left_bp: usize,
}

impl std::fmt::Display for DissolveLedger {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "dissolved {} small, {} duplicated and {} clean bins holding {} bp; pool {} contigs \
             {} bp; {} rounds ending at rung {}; proposed {} clusters, refused {} small, {} \
             duplicated and {} no better, {} noise; promoted {} bins adopting {} contigs {} bp; \
             returned {} contigs {} bp, emptied {} bins; left {} contigs {} bp unbinned",
            self.dissolved_small,
            self.dissolved_duplicated,
            self.dissolved_clean,
            self.dissolved_bp,
            self.pool_contigs,
            self.pool_bp,
            self.rounds,
            self.rung,
            self.proposed,
            self.refused_small,
            self.refused_duplicated,
            self.refused_worse,
            self.noise,
            self.promoted,
            self.adopted_contigs,
            self.adopted_bp,
            self.returned_contigs,
            self.returned_bp,
            self.emptied,
            self.left_contigs,
            self.left_bp
        )
    }
}

/// A piece smaller than the run's own genome scale is a shard of one, and relaxing that to let
/// more of the pool through cost more bins than it recovered.
fn floor_for(settings: DissolveSettings) -> usize {
    settings
        .genome_floor
        .unwrap_or(settings.min_bin_size)
        .max(settings.min_bin_size)
}

fn rung_of(settings: DissolveSettings, top: usize, rung: usize) -> (usize, f64) {
    let (share, multiple) = LADDER[rung];
    let floor = settings.min_bin_size
        + (share * top.saturating_sub(settings.min_bin_size) as f64).round() as usize;
    (floor, settings.duplication_bar * multiple)
}

fn quality_rung(settings: DissolveSettings, rung: usize) -> (f64, f64) {
    let (share, multiple) = QUALITY_LADDER[rung];
    (
        settings.completeness * share,
        settings.contamination * multiple,
    )
}

fn over_bar(features: &ContigFeatures, contigs: &[usize], bar: f64) -> bool {
    features
        .sketches()
        .and_then(|sketches| sketches.duplication(contigs))
        .is_some_and(|duplication| duplication > bar)
}

enum Verdict {
    Adopt,
    TooSmall,
    Duplicated,
}

fn judge(
    features: &ContigFeatures,
    quality: Option<&ContigQuality>,
    contigs: &[usize],
    rung: Rung,
) -> Verdict {
    if features.bin_size(contigs) < rung.floor {
        return Verdict::TooSmall;
    }
    match quality {
        Some(quality) => {
            let held = quality.score(contigs);
            match held.completeness >= rung.completeness && held.contamination <= rung.contamination
            {
                true => Verdict::Adopt,
                false => Verdict::Duplicated,
            }
        }
        None => match over_bar(features, contigs, rung.bar) {
            true => Verdict::Duplicated,
            false => Verdict::Adopt,
        },
    }
}

#[derive(Debug, Clone, Copy)]
struct Rung {
    floor: usize,
    bar: f64,
    completeness: f64,
    contamination: f64,
}

/// Completeness is its own size test, so a small genome is not held to the run's genome scale
/// once the gene families can say it is whole.
fn rung(settings: DissolveSettings, top: usize, at: usize, scored: bool) -> Rung {
    let (floor, bar) = rung_of(settings, top, at);
    let (completeness, contamination) = quality_rung(settings, at);
    Rung {
        floor: match scored {
            true => settings.min_bin_size,
            false => floor,
        },
        bar,
        completeness,
        contamination,
    }
}

fn sorted(contigs: HashSet<usize>) -> Vec<usize> {
    let mut contigs = contigs.into_iter().collect::<Vec<_>>();
    contigs.sort_unstable();
    contigs
}

fn bases(features: &ContigFeatures, contigs: &HashSet<usize>) -> usize {
    contigs.iter().map(|contig| features.length(*contig)).sum()
}

/// The eject arm has already stripped what it could before this runs, so a bin still over the
/// bar is one it could not fix, which is the bin least worth trusting as it stands.
fn dissolving(
    features: &ContigFeatures,
    bins: &BTreeMap<usize, Vec<usize>>,
    top: usize,
    settings: DissolveSettings,
    ledger: &mut DissolveLedger,
) -> Vec<(usize, Vec<usize>)> {
    let mut dissolving = Vec::new();
    for (bin_id, contigs) in bins.iter() {
        let small = features.bin_size(contigs) < top;
        let duplicated = !small && over_bar(features, contigs, settings.duplication_bar);
        if settings.scope == DissolveScope::Fused && !small && !duplicated {
            continue;
        }
        if small {
            ledger.dissolved_small += 1;
        } else if duplicated {
            ledger.dissolved_duplicated += 1;
        } else {
            ledger.dissolved_clean += 1;
        }
        ledger.dissolved_bp += features.bin_size(contigs);
        dissolving.push((*bin_id, contigs.clone()));
    }
    dissolving
}

struct Pot<'a> {
    features: &'a ContigFeatures<'a>,
    quality: Option<&'a ContigQuality>,
    origin: HashMap<usize, usize>,
    held: HashMap<usize, (f64, usize)>,
}

impl Pot<'_> {
    /// A cluster that takes the greater part of a bin has to be the better bin, or the loop
    /// trades a whole genome for a piece of one.
    fn improves(&self, contigs: &[usize]) -> bool {
        let mut taken: HashMap<usize, usize> = HashMap::new();
        for contig in contigs {
            if let Some(bin) = self.origin.get(contig) {
                *taken.entry(*bin).or_default() += self.features.length(*contig);
            }
        }
        let Some(quality) = self.quality else {
            return true;
        };
        let candidate = quality.score(contigs).score();
        taken
            .into_iter()
            .all(|(bin, bases)| match self.held.get(&bin) {
                Some((score, whole)) if bases * 2 >= *whole => candidate > *score,
                _ => true,
            })
    }

    fn adopt(&self, clusters: &[Vec<usize>], rung: Rung, improve: bool) -> Taken {
        let mut taken = Taken::default();
        for contigs in clusters {
            match judge(self.features, self.quality, contigs, rung) {
                Verdict::Adopt if improve && !self.improves(contigs) => taken.refused_worse += 1,
                Verdict::Adopt => taken.clusters.push(contigs.clone()),
                Verdict::TooSmall => taken.refused_small += 1,
                Verdict::Duplicated => taken.refused_duplicated += 1,
            }
        }
        taken
    }
}

#[derive(Default)]
struct Taken {
    clusters: Vec<Vec<usize>>,
    refused_small: usize,
    refused_duplicated: usize,
    refused_worse: usize,
}

fn neighbours_for(settings: DissolveSettings, round: usize) -> RoundParams {
    RoundParams {
        n_neighbours: settings
            .n_neighbours
            .checked_shr(round as u32)
            .unwrap_or(0)
            .max(MIN_NEIGHBOURS),
    }
}

fn search(
    pot: &Pot,
    pool: &mut HashSet<usize>,
    settings: DissolveSettings,
    top: usize,
    ledger: &mut DissolveLedger,
    partition: impl Fn(&HashSet<usize>, RoundParams) -> Result<HDBSCANResult>,
) -> Vec<Vec<usize>> {
    if settings.select == Selection::Ranked {
        return best_of_every_round(pot, pool, settings, top, ledger, partition);
    }
    let mut promoted = Vec::new();
    for round in 0..settings.rounds.max(1) {
        if pool.len() < settings.min_contigs {
            break;
        }
        let params = neighbours_for(settings, round);
        let result = match partition(pool, params) {
            Ok(result) => result,
            Err(error) => {
                warn!("Could not re-embed the pool: {error}");
                break;
            }
        };
        ledger.rounds += 1;
        ledger.noise = result.outliers.len();
        let clusters = result
            .cluster_map
            .into_values()
            .map(sorted)
            .collect::<Vec<_>>();

        let scored = pot.quality.is_some();
        let mut taken = pot.adopt(
            &clusters,
            rung(settings, top, ledger.rung, scored),
            settings.improve,
        );
        while taken.clusters.is_empty() && settings.ladder && ledger.rung + 1 < LADDER.len() {
            ledger.rung += 1;
            taken = pot.adopt(
                &clusters,
                rung(settings, top, ledger.rung, scored),
                settings.improve,
            );
        }
        ledger.proposed += clusters.len();
        ledger.refused_small += taken.refused_small;
        ledger.refused_duplicated += taken.refused_duplicated;
        ledger.refused_worse += taken.refused_worse;
        if taken.clusters.is_empty() {
            break;
        }
        for contigs in &taken.clusters {
            for contig in contigs {
                pool.remove(contig);
            }
        }
        promoted.extend(taken.clusters);
    }
    promoted
}

fn proposals(
    pool: &HashSet<usize>,
    settings: DissolveSettings,
    ledger: &mut DissolveLedger,
    partition: impl Fn(&HashSet<usize>, RoundParams) -> Result<HDBSCANResult>,
) -> Vec<Vec<usize>> {
    let mut candidates = Vec::new();
    for round in 0..settings.rounds.max(1) {
        let result = match partition(pool, neighbours_for(settings, round)) {
            Ok(result) => result,
            Err(error) => {
                warn!("Could not re-embed the pool: {error}");
                break;
            }
        };
        ledger.rounds += 1;
        ledger.noise = result.outliers.len();
        candidates.extend(result.cluster_map.into_values().map(sorted));
    }
    candidates
}

/// Every round searches the same pool, so a genome only one k finds is proposed alongside the
/// blob that swallows it, and the bar is asked which of the two to keep rather than which came first.
fn best_of_every_round(
    pot: &Pot,
    pool: &mut HashSet<usize>,
    settings: DissolveSettings,
    top: usize,
    ledger: &mut DissolveLedger,
    partition: impl Fn(&HashSet<usize>, RoundParams) -> Result<HDBSCANResult>,
) -> Vec<Vec<usize>> {
    if pool.len() < settings.min_contigs {
        return Vec::new();
    }
    let candidates = proposals(pool, settings, ledger, partition);
    ledger.proposed += candidates.len();
    let order = ranked(pot.features, pot.quality, candidates);
    let scored = pot.quality.is_some();

    let mut promoted = Vec::new();
    let mut claimed = HashSet::new();
    loop {
        let bar = rung(settings, top, ledger.rung, scored);
        let mut taken = Taken::default();
        for (contigs, _) in &order {
            let left = remaining(contigs, &claimed);
            if left.len() < 2 {
                continue;
            }
            match judge(pot.features, pot.quality, &left, bar) {
                Verdict::Adopt if settings.improve && !pot.improves(&left) => {
                    taken.refused_worse += 1
                }
                Verdict::Adopt => {
                    claimed.extend(left.iter().copied());
                    taken.clusters.push(left);
                }
                Verdict::TooSmall => taken.refused_small += 1,
                Verdict::Duplicated => taken.refused_duplicated += 1,
            }
        }
        ledger.refused_small += taken.refused_small;
        ledger.refused_duplicated += taken.refused_duplicated;
        ledger.refused_worse += taken.refused_worse;
        let empty = taken.clusters.is_empty();
        promoted.extend(taken.clusters);
        if !empty || !settings.ladder || ledger.rung + 1 >= LADDER.len() {
            break;
        }
        ledger.rung += 1;
    }
    for contig in &claimed {
        pool.remove(contig);
    }
    promoted
}

/// Bins the sketch says hold their own sequence twice, and bins under the genome floor, go back
/// in the pot with the unbinned and are searched again without the bins that already left, which
/// is the one thing re-cutting inside a bin cannot do.
pub fn dissolve(
    features: &ContigFeatures,
    quality: Option<&ContigQuality>,
    bins: &mut BTreeMap<usize, Vec<usize>>,
    unbinned: &mut Vec<usize>,
    settings: DissolveSettings,
    partition: impl Fn(&HashSet<usize>, RoundParams) -> Result<HDBSCANResult>,
) -> DissolveLedger {
    let mut ledger = DissolveLedger::default();
    let top = floor_for(settings);
    let dissolved = dissolving(features, bins, top, settings, &mut ledger);

    let mut pool = unbinned.iter().copied().collect::<HashSet<_>>();
    for (_, contigs) in &dissolved {
        pool.extend(contigs.iter().copied());
    }
    ledger.pool_contigs = pool.len();
    ledger.pool_bp = bases(features, &pool);

    if pool.len() < settings.min_contigs {
        return ledger;
    }
    let handed = pool.clone();
    let pot = Pot {
        features,
        quality,
        origin: dissolved
            .iter()
            .flat_map(|(bin_id, contigs)| contigs.iter().map(|contig| (*contig, *bin_id)))
            .collect(),
        held: match (settings.improve, quality) {
            (true, Some(quality)) => dissolved
                .iter()
                .map(|(bin_id, contigs)| {
                    (
                        *bin_id,
                        (quality.score(contigs).score(), features.bin_size(contigs)),
                    )
                })
                .collect(),
            _ => HashMap::new(),
        },
    };
    let mut promoted = search(&pot, &mut pool, settings, top, &mut ledger, partition);
    if promoted.is_empty() {
        return ledger;
    }
    // Hash order would give the same partition different bin names on every run.
    promoted.sort_unstable_by_key(|contigs| contigs.first().copied().unwrap_or(usize::MAX));

    if let Err(error) = placed_once(promoted.iter().flatten().copied(), &handed) {
        warn!("The pool did not come back partitioned, leaving the bins alone: {error}");
        return ledger;
    }

    let next_bin_id = bins.keys().max().map_or(0, |id| id + 1);
    let mut claimed = HashSet::new();
    ledger.promoted = promoted.len();
    for (offset, contigs) in promoted.into_iter().enumerate() {
        ledger.adopted_contigs += contigs.len();
        ledger.adopted_bp += features.bin_size(&contigs);
        claimed.extend(contigs.iter().copied());
        bins.insert(next_bin_id + offset, contigs);
    }

    // Whatever the pool did not claim goes back to the bin it came from, or a long contig the
    // writer would stand alone leaves as a singleton for no gain.
    for (bin_id, contigs) in dissolved {
        let kept = contigs
            .iter()
            .copied()
            .filter(|contig| !claimed.contains(contig))
            .collect::<Vec<_>>();
        if kept.is_empty() {
            ledger.emptied += 1;
            bins.remove(&bin_id);
            continue;
        }
        ledger.returned_contigs += kept.len();
        ledger.returned_bp += features.bin_size(&kept);
        for contig in &kept {
            pool.remove(contig);
        }
        bins.insert(bin_id, kept);
    }

    ledger.left_contigs = pool.len();
    ledger.left_bp = bases(features, &pool);
    *unbinned = sorted(pool);
    ledger
}
