use rosella::markers::eukaryote::{Reading, eukaryotic};
use rosella::quality::{Bars, Quality};

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
