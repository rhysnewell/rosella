//! Counting communities over a 50 per cent bar rewards splitting: a genome cut in two scores
//! twice where it scored once, so the statistic peaks past the true genome count.

use std::collections::{HashMap, HashSet};

use rosella::clustering::clusterer::Partitioning;
use rosella::quality::{Quality, Scorer, Worth};
use rosella::recover::ladder::{Judge, RungStatistic, pick_rung};

const GENOME: usize = 4;
const CONTIGS: usize = 8;

struct Planted;

impl Scorer for Planted {
    fn score(&self, contigs: &[usize]) -> Quality {
        let mut per_genome = [0usize; 2];
        for contig in contigs {
            per_genome[contig / GENOME] += 1;
        }
        let whole = per_genome[0].max(per_genome[1]);
        let foreign = per_genome[0] + per_genome[1] - whole;
        Quality {
            completeness: 100.0 * whole as f64 / GENOME as f64,
            contamination: 100.0 * foreign as f64 / GENOME as f64,
        }
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        contigs.iter().map(|contig| *contig as u32).collect()
    }
}

fn partitioning(clusters: &[&[usize]]) -> Partitioning {
    Partitioning {
        cluster_map: clusters
            .iter()
            .enumerate()
            .map(|(id, contigs)| (id, contigs.iter().copied().collect::<HashSet<_>>()))
            .collect::<HashMap<_, _>>(),
        outliers: HashSet::new(),
        score: 0.0,
    }
}

fn ladder() -> Vec<Partitioning> {
    vec![
        partitioning(&[&[0, 1, 2, 3], &[4, 5, 6, 7]]),
        partitioning(&[&[0, 1], &[2, 3], &[4, 5], &[6, 7]]),
    ]
}

fn chosen(statistic: RungStatistic) -> usize {
    let contigs = (0..CONTIGS).collect::<Vec<_>>();
    let scorer = Planted;
    let judge = Judge {
        quality: &scorer,
        contigs: &contigs,
        worth: Worth { contamination: 2.0, allowance: 0.0 },
        rung_statistic: statistic,
    };
    pick_rung(ladder(), &judge).cluster_map.len()
}

#[test]
fn pass50_takes_the_over_split_rung_and_pass90_takes_the_whole_genomes() {
    assert_eq!(chosen(RungStatistic::Pass50), 4);
    assert_eq!(chosen(RungStatistic::Pass90), 2);
    assert_eq!(chosen(RungStatistic::Pass80), 2);
}