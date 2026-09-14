//! Both of these bins pass a 5 per cent accept bar, so ranking them on contamination charged from
//! zero is what lets the cleaner piece take the whole genome's contigs.

use std::collections::{HashMap, HashSet};

use rosella::clustering::clusterer::Partitioning;
use rosella::clustering::graph_partition::Partition;
use rosella::quality::{Quality, Scorer, Worth};
use rosella::recover::ladder::{Judge, combine};
use rosella::refine::rung::Bars;

#[path = "../support/bars.rs"]
mod bars;

struct Table;

impl Scorer for Table {
    fn score(&self, contigs: &[usize]) -> Quality {
        let (completeness, contamination) = match contigs {
            [0, 1, 2, 3] => (100.0, 4.0),
            [0, 1] => (95.0, 0.0),
            _ => (10.0, 0.0),
        };
        Quality {
            completeness,
            contamination,
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
        score: None,
        arm: Partition::Leiden,
        seed: 0,
    }
}

fn combined(allowance: f64) -> usize {
    let contigs = (0..4).collect::<Vec<_>>();
    let scorer = Table;
    let judge = Judge {
        quality: &scorer,
        contigs: &contigs,
        bars: Bars {
            worth: Worth {
                contamination: 2.0,
                allowance,
            },
            ..bars::bars(80.0)
        },
        rungs: 0,
        size_tie: false,
    };
    let ladder = vec![
        partitioning(&[&[0, 1, 2, 3]]),
        partitioning(&[&[0, 1], &[2, 3]]),
    ];
    combine(ladder, &judge).cluster_map.len()
}

#[test]
fn an_allowance_keeps_the_whole_genome_the_clean_piece_would_have_cut() {
    assert_eq!(
        combined(0.0),
        2,
        "charged from zero, the 95/0 piece outranks the 100/4 genome"
    );
    assert_eq!(combined(5.0), 1);
}
