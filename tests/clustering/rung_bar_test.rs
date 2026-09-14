//! Counting communities at 50 per cent counts a halved genome twice, which is what took the judge
//! to the 875 community rung on cami_i_high where the truth is 1075 genomes and the coarse rung
//! scores 35 more bins.

use std::collections::{HashMap, HashSet};

use rosella::clustering::clusterer::Partitioning;
use rosella::clustering::graph_partition::Partition;
use rosella::quality::{Quality, Scorer, Worth};
use rosella::recover::ladder::{Judge, pick_rung};

const GENOME: usize = 4;
const CONTIGS: usize = 8;

struct Planted {
    offset: f64,
}

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

    fn completeness_bar(&self, requested: f64) -> f64 {
        (requested - self.offset).max(0.0)
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
        score: None,
        arm: Partition::Leiden,
    }
}

fn chosen(offset: f64) -> usize {
    let contigs = (0..CONTIGS).collect::<Vec<_>>();
    let scorer = Planted { offset };
    let judge = Judge {
        quality: &scorer,
        contigs: &contigs,
        worth: Worth {
            contamination: 2.0,
            allowance: 0.0,
        },
        completeness: 90.0,
    };
    let ladder = vec![
        partitioning(&[&[0, 1, 2, 3], &[4, 5, 6, 7]]),
        partitioning(&[&[0, 1], &[2, 3], &[4, 5], &[6, 7]]),
    ];
    pick_rung(ladder, &judge).cluster_map.len()
}

#[test]
fn the_judge_counts_at_the_scorers_own_accept_bar() {
    assert_eq!(chosen(10.0), 2, "an 80 per cent bar refuses the halves");
    assert_eq!(chosen(40.0), 4, "a 50 per cent bar counts each half as a genome");
}
