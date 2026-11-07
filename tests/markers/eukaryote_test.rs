use rosella::markers::eukaryote::{Reading, eukaryotic};
use rosella::markers::replicon::Shape;
use rosella::markers::{ContigMarkers, Hit, MarkerSet};
use rosella::quality::{Bars, Quality, Scorer};

const BARS: Bars = Bars {
    completeness: 90.0,
    contamination: 5.0,
};

fn bin(density: f64, completeness: f64) -> Reading {
    Reading {
        density: Some(density),
        quality: Quality {
            completeness,
            contamination: 0.0,
            set: 0,
        },
    }
}

fn finished(count: usize) -> Vec<Reading> {
    (0..count)
        .map(|at| bin(0.86 + 0.01 * at as f64, 95.0))
        .collect()
}

#[test]
fn only_a_sparse_bin_short_of_medium_quality_is_eukaryotic() {
    let mut readings = finished(6);
    readings.extend([
        bin(0.55, 20.0),
        bin(0.55, 60.0),
        bin(0.80, 10.0),
        Reading {
            density: None,
            quality: Quality::default(),
        },
    ]);
    let called = eukaryotic(&readings, BARS);
    assert_eq!(
        called[6..],
        [true, false, false, false],
        "fence sits under the finished bins and a medium quality bin is a prokaryote"
    );
    assert!(!called[..6].iter().any(|called| *called));
}

#[test]
fn too_few_finished_bins_name_nothing() {
    let mut readings = finished(3);
    readings.push(bin(0.20, 0.0));
    assert!(!eukaryotic(&readings, BARS).iter().any(|called| *called));
}

#[test]
fn a_sparse_carrier_falls_under_the_fence_and_an_unread_one_stays() {
    let mut shapes = (0..6)
        .map(|at| {
            Some(Shape {
                coding_bases: 860 + 10 * at,
                genes: 1,
            })
        })
        .collect::<Vec<_>>();
    shapes.extend([
        Some(Shape {
            coding_bases: 550,
            genes: 1,
        }),
        None,
    ]);
    let markers = ContigMarkers::new(vec![Vec::new(); 8], MarkerSet::joined(None))
        .with_lengths(vec![1000; 8])
        .with_shapes(shapes);
    let kept = markers.prokaryotic(&(0..8).collect::<Vec<_>>());
    assert_eq!(kept, [true, true, true, true, true, true, false, true]);
}

#[test]
fn a_marker_on_a_sparse_carrier_counts_as_present_but_never_as_a_copy() {
    let hit = |marker| Hit {
        marker,
        ..Default::default()
    };
    let dense = |at: u64| Shape {
        coding_bases: 860 + 10 * at,
        genes: 1,
    };
    let mut shapes = (0..6).map(|at| Some(dense(at))).collect::<Vec<_>>();
    shapes.push(Some(Shape {
        coding_bases: 550,
        genes: 1,
    }));
    let mut per_contig = vec![vec![hit(0)]; 6];
    per_contig.push(vec![hit(0), hit(1)]);
    let table = "model_name\tsets\nalpha\tbac\nbeta\tbac\n";
    let mut markers = ContigMarkers::new(per_contig, MarkerSet::parse_with(table, None))
        .with_lengths(vec![1000; 7])
        .with_shapes(shapes);
    assert_eq!(markers.score(&[0, 6]).contamination, 50.0);
    markers.count_sparse_as_presence();
    let scored = markers.score(&[0, 6]);
    assert_eq!((scored.completeness, scored.contamination), (100.0, 0.0));
    assert_eq!(markers.score(&[0, 1]).contamination, 50.0);
}

#[test]
fn pieces_of_one_genome_split_their_families_and_two_genomes_share_them() {
    let hits = |markers: &[u16]| {
        markers
            .iter()
            .map(|marker| Hit {
                marker: *marker,
                ..Default::default()
            })
            .collect::<Vec<_>>()
    };
    let table = "model_name\tsets\nalpha\tbac\nbeta\tbac\ngamma\tbac\ndelta\tbac\n";
    let markers = ContigMarkers::new(
        vec![hits(&[0, 1]), hits(&[2, 3]), hits(&[0, 1, 2]), Vec::new()],
        MarkerSet::parse_with(table, None),
    )
    .with_lengths(vec![1000; 4]);
    assert_eq!(markers.shared(&[0], &[1]), Some(0.0));
    assert_eq!(markers.shared(&[0], &[2]), Some(1.0));
    assert_eq!(
        markers.shared(&[0], &[3]),
        None,
        "a piece with no markers says nothing"
    );
}
