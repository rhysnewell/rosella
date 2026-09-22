use rosella::markers::replicon::{Shape, small_replicons};

const CHROMOSOME_GENE: usize = 1000;
const ELEMENT_GENE: usize = 500;

fn contig(genes: u32, gene_bases: usize) -> Shape {
    let mut shape = Shape::default();
    for _ in 0..genes {
        shape.add(gene_bases);
    }
    shape
}

fn assembly(extra: &[(Shape, bool, usize)]) -> (Vec<Shape>, Vec<bool>, Vec<usize>) {
    let mut shapes = (0..40)
        .map(|at| contig(400, CHROMOSOME_GENE - at * 5))
        .collect::<Vec<_>>();
    let mut carries = vec![true; 40];
    let mut lengths = vec![450_000; 40];
    for (shape, carrier, length) in extra {
        shapes.push(*shape);
        carries.push(*carrier);
        lengths.push(*length);
    }
    (shapes, carries, lengths)
}

#[test]
fn a_dense_short_gene_contig_with_no_marker_is_a_replicon() {
    let (shapes, carries, lengths) = assembly(&[(contig(60, ELEMENT_GENE), false, 32_000)]);
    assert_eq!(small_replicons(&shapes, &carries, &lengths), vec![40]);
}

#[test]
fn a_marker_carrier_is_never_a_replicon() {
    let (shapes, carries, lengths) = assembly(&[(contig(60, ELEMENT_GENE), true, 32_000)]);
    assert!(small_replicons(&shapes, &carries, &lengths).is_empty());
}

#[test]
fn a_sparse_contig_is_not_a_replicon() {
    let (shapes, carries, lengths) = assembly(&[(contig(30, ELEMENT_GENE), false, 90_000)]);
    assert!(small_replicons(&shapes, &carries, &lengths).is_empty());
}

#[test]
fn chromosome_length_genes_are_not_a_replicon() {
    let (shapes, carries, lengths) = assembly(&[(contig(30, CHROMOSOME_GENE), false, 32_000)]);
    assert!(small_replicons(&shapes, &carries, &lengths).is_empty());
}

#[test]
fn a_contig_under_the_size_floor_is_not_a_replicon() {
    let (shapes, carries, lengths) = assembly(&[(contig(60, ELEMENT_GENE), false, 9_000)]);
    assert!(small_replicons(&shapes, &carries, &lengths).is_empty());
}

#[test]
fn too_few_genes_to_read_a_shape_is_not_a_replicon() {
    let (shapes, carries, lengths) = assembly(&[(contig(4, ELEMENT_GENE), false, 12_000)]);
    assert!(small_replicons(&shapes, &carries, &lengths).is_empty());
}

#[test]
fn an_assembly_with_too_few_anchors_flags_nothing() {
    let shapes = [contig(400, CHROMOSOME_GENE), contig(60, ELEMENT_GENE)].to_vec();
    let carries = [true, false].to_vec();
    let lengths = [450_000, 32_000].to_vec();
    assert!(small_replicons(&shapes, &carries, &lengths).is_empty());
}
