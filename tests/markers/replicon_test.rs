use ndarray::Array2;
use rosella::markers::replicon::{Shape, small_replicons};

const CHROMOSOME_GENE: usize = 1000;
const ELEMENT_GENE: usize = 500;
const ANCHORS: usize = 40;
const FOREIGN: f64 = 5.0;
const OWN: f64 = 0.001;

fn contig(genes: u32, gene_bases: usize) -> Shape {
    let mut shape = Shape::default();
    for _ in 0..genes {
        shape.add(gene_bases);
    }
    shape
}

struct Candidate {
    shape: Shape,
    carrier: bool,
    length: usize,
    offset: f64,
}

fn element() -> Candidate {
    Candidate {
        shape: contig(60, ELEMENT_GENE),
        carrier: false,
        length: 32_000,
        offset: FOREIGN,
    }
}

fn run(candidate: Candidate, anchors: usize) -> Vec<usize> {
    let mut shapes = (0..anchors)
        .map(|at| contig(400, CHROMOSOME_GENE - at * 5))
        .collect::<Vec<_>>();
    let mut carries = vec![true; anchors];
    let mut lengths = vec![450_000; anchors];
    shapes.push(candidate.shape);
    carries.push(candidate.carrier);
    lengths.push(candidate.length);
    let composition = Array2::from_shape_fn((anchors + 1, 2), |(row, column)| match row {
        row if row == anchors => candidate.offset,
        row => (row + column) as f64 * OWN,
    });
    let bin = (0..=anchors).collect::<Vec<_>>();
    small_replicons(
        &shapes,
        &carries,
        &lengths,
        [bin.as_slice()],
        composition.view(),
    )
}

#[test]
fn a_short_gene_contig_foreign_to_its_bin_is_a_replicon() {
    assert_eq!(run(element(), ANCHORS), vec![ANCHORS]);
}

#[test]
fn a_contig_that_fits_its_bin_is_a_fragment() {
    let fits = Candidate {
        offset: OWN,
        ..element()
    };
    assert!(run(fits, ANCHORS).is_empty());
}

#[test]
fn only_a_marker_free_dense_short_gene_contig_over_the_floor_qualifies() {
    let refused = [
        Candidate {
            carrier: true,
            ..element()
        },
        Candidate {
            shape: contig(30, ELEMENT_GENE),
            length: 90_000,
            ..element()
        },
        Candidate {
            shape: contig(30, CHROMOSOME_GENE),
            ..element()
        },
        Candidate {
            length: 9_000,
            ..element()
        },
        Candidate {
            shape: contig(4, ELEMENT_GENE),
            length: 12_000,
            ..element()
        },
    ];
    for candidate in refused {
        assert!(run(candidate, ANCHORS).is_empty());
    }
}

#[test]
fn an_assembly_with_too_few_anchors_flags_nothing() {
    assert!(run(element(), 2).is_empty());
}
