use std::collections::{BTreeMap, HashSet};

use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::quality::{Quality, Scorer};
use rosella::recover::coherence::Coherence;
use rosella::recover::strains::{Depths, split_strains};

const MARKERS: u32 = 10;

struct Markers(Vec<Option<u32>>);

impl Markers {
    fn counts(&self, contigs: &[usize]) -> BTreeMap<u32, usize> {
        let mut counts = BTreeMap::new();
        for marker in contigs.iter().filter_map(|contig| self.0[*contig]) {
            *counts.entry(marker).or_insert(0) += 1;
        }
        counts
    }

    fn carriers(&self, contigs: &[usize]) -> Vec<Vec<u32>> {
        let mut families = BTreeMap::<u32, Vec<u32>>::new();
        for (at, contig) in contigs.iter().enumerate() {
            if let Some(marker) = self.0[*contig] {
                families.entry(marker).or_default().push(at as u32);
            }
        }
        families
            .into_values()
            .filter(|family| family.len() > 1)
            .collect()
    }
}

impl Scorer for Markers {
    fn score(&self, contigs: &[usize]) -> Quality {
        let counts = self.counts(contigs);
        let extra = counts.values().map(|count| count - 1).sum::<usize>();
        Quality {
            completeness: 100.0 * counts.len() as f64 / f64::from(MARKERS),
            contamination: 100.0 * extra as f64 / f64::from(MARKERS),
            ..Default::default()
        }
    }

    fn features(&self, contigs: &[usize]) -> HashSet<u32> {
        self.counts(contigs).into_keys().collect()
    }
}

struct World {
    depths: Array2<f64>,
    markers: Vec<Option<u32>>,
}

const JITTER: [f64; 7] = [1.0, 1.08, 0.93, 1.04, 0.96, 1.02, 0.98];

impl World {
    fn new() -> Self {
        Self {
            depths: Array2::zeros((0, 2)),
            markers: Vec::new(),
        }
    }

    fn genome(&mut self, profile: [f64; 2], contigs: usize, markers: &[u32]) -> Vec<usize> {
        let start = self.depths.nrows();
        for at in 0..contigs {
            let scale = JITTER[at % JITTER.len()];
            self.depths
                .push_row(ndarray::aview1(&[profile[0] * scale, profile[1] * scale]))
                .expect("two samples");
            self.markers.push(markers.get(at).copied());
        }
        (start..start + contigs).collect()
    }
}

fn run(world: &World, hosts: &[Vec<usize>], bin: Vec<usize>) -> BTreeMap<usize, Vec<usize>> {
    let lengths = vec![20_000; world.depths.nrows()];
    let hosts = hosts.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let coherence = Coherence::of(&hosts, world.depths.view(), &lengths)
        .expect("three hosts measure the spread");
    let markers = Markers(world.markers.clone());
    let depths = Depths {
        depths: world.depths.view(),
        lengths: &lengths,
        read_lengths: Some(&[4_000.0, 4_000.0]),
    };
    let rows = world.depths.nrows();
    let mut coverage = Array2::zeros((rows, 4));
    for row in 0..rows {
        coverage[[row, 0]] = world.depths[[row, 0]];
        coverage[[row, 2]] = world.depths[[row, 1]];
    }
    let composition = Array2::from_elem((rows, 4), 0.25);
    let features = ContigFeatures::new(&coverage, &composition, &lengths);
    let mut bins = BTreeMap::from([(0, bin)]);
    split_strains(
        &mut bins,
        &markers,
        |members| markers.carriers(members),
        &depths,
        &coherence,
        &features,
    );
    bins
}

fn hosts(world: &mut World) -> Vec<Vec<usize>> {
    let all = (0..MARKERS).collect::<Vec<_>>();
    [[3.0, 30.0], [40.0, 2.0], [15.0, 15.0]]
        .iter()
        .map(|profile| world.genome(*profile, 12, &all))
        .collect()
}

#[test]
fn two_strains_part_where_their_depths_differ() {
    let mut world = World::new();
    let hosts = hosts(&mut world);
    let all = (0..MARKERS).collect::<Vec<_>>();
    let high = world.genome([12.0, 6.0], 20, &all);
    let low = world.genome([4.0, 6.0], 20, &all);
    let bins = run(&world, &hosts, [high.clone(), low.clone()].concat());
    let mut parts = bins.into_values().collect::<Vec<_>>();
    parts.sort();
    assert_eq!(parts, vec![high, low]);
}

#[test]
fn one_genome_with_a_doubled_marker_stays_whole() {
    let mut world = World::new();
    let hosts = hosts(&mut world);
    let doubled = [&(0..MARKERS).collect::<Vec<_>>()[..], &[0]].concat();
    let genome = world.genome([8.0, 6.0], 30, &doubled);
    let bins = run(&world, &hosts, genome.clone());
    assert_eq!(bins.into_values().collect::<Vec<_>>(), vec![genome]);
}
