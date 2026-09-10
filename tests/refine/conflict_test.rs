use std::collections::BTreeMap;

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::markers::{ContigMarkers, Hit, MarkerRules, MarkerSet};
use rosella::refine::conflict::{Bars, eject_conflicts};

const TABLE: &str = "model_name\tdomain\n\
                     alpha\tbac120\n\
                     beta\tbac120\n\
                     gamma\tbac120\n\
                     delta\tbac120\n";

const BARS: Bars = Bars {
    completeness: 80.0,
    contamination: 5.0,
};
const FLOOR: usize = 200_000;

struct Bin {
    lengths: Vec<usize>,
    coverage: Array2<f64>,
    tnf: Array2<f64>,
    markers: ContigMarkers,
}

fn bin(lengths: Vec<usize>, hits: Vec<Vec<u16>>) -> Bin {
    let rows = lengths.len();
    let per_contig = hits
        .into_iter()
        .map(|held| {
            held.into_iter()
                .map(|marker| Hit {
                    marker,
                    partial: false,
                })
                .collect()
        })
        .collect();
    Bin {
        coverage: Array2::from_elem((rows, 2), 10.0),
        tnf: Array2::from_elem((rows, 4), 0.25),
        markers: ContigMarkers::new(per_contig, MarkerSet::parse(TABLE), MarkerRules::default()),
        lengths,
    }
}

impl Bin {
    fn features(&self) -> ContigFeatures<'_> {
        ContigFeatures::new(&self.coverage, &self.tnf, &self.lengths)
    }

    fn run(&self, contigs: Vec<usize>) -> (Vec<usize>, Vec<usize>) {
        let mut bins = BTreeMap::from([(0, contigs)]);
        let (ejected, _) = eject_conflicts(&self.features(), &self.markers, &mut bins, BARS, FLOOR);
        (ejected, bins.remove(&0).unwrap_or_default())
    }
}

#[test]
fn the_rider_holding_a_second_copy_leaves_and_the_whole_genome_stays() {
    let held = bin(vec![1_000_000, 100_000], vec![vec![0, 1, 2, 3], vec![0]]);
    assert_eq!(held.run(vec![0, 1]), (vec![1], vec![0]));
}

#[test]
fn a_rider_carrying_the_only_copy_of_a_marker_stays_however_small() {
    let held = bin(vec![1_000_000, 100_000], vec![vec![0, 1, 2], vec![0, 3]]);
    assert_eq!(held.run(vec![0, 1]), (Vec::new(), vec![0, 1]));
}

#[test]
fn the_bin_floor_stops_the_peel_before_the_bar_is_reached() {
    let held = bin(vec![150_000, 100_000], vec![vec![0, 1, 2, 3], vec![0]]);
    assert_eq!(held.run(vec![0, 1]), (Vec::new(), vec![0, 1]));
}

#[test]
fn a_bin_the_scorer_calls_short_is_left_alone_however_dirty() {
    let held = bin(vec![1_000_000, 100_000], vec![vec![0, 1], vec![0]]);
    assert_eq!(held.run(vec![0, 1]), (Vec::new(), vec![0, 1]));
}
