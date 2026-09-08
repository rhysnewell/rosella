use ndarray::Array2;
use rosella::embedding::features::ContigFeatures;
use rosella::embedding::umap::EmbedOverrides;
use rosella::homology::{Homology, HomologySettings, Pair};
use rosella::seeds::Seeds;

const CONTIGS: usize = 24;
const LENGTH: usize = 400_000;
const NEIGHBOURS: usize = 5;

fn seeds() -> Seeds {
    Seeds {
        knn: 42,
        init: 42,
        layout: 42,
        sample: 42,
        partition: 42,
    }
}

fn two_blobs() -> (Array2<f64>, Array2<f64>, Vec<usize>) {
    let mut coverage = Array2::zeros((CONTIGS, 2));
    let mut tnf = Array2::zeros((CONTIGS, 2));
    for contig in 0..CONTIGS {
        let blob = if contig < CONTIGS / 2 { 0.0 } else { 10.0 };
        let position = blob + contig as f64 * 0.01;
        coverage[[contig, 0]] = position;
        coverage[[contig, 1]] = 1.0;
        tnf[[contig, 0]] = position;
        tnf[[contig, 1]] = 0.5;
    }
    (coverage, tnf, vec![LENGTH; CONTIGS])
}

#[test]
fn a_severed_edge_leaves_the_graph_rather_than_holding_a_zero() {
    let (coverage, tnf, lengths) = two_blobs();
    let indices = (0..CONTIGS).collect::<Vec<_>>();
    let overrides = EmbedOverrides::default();

    let open = ContigFeatures::new(&coverage, &tnf, &lengths)
        .graph_of(&indices, NEIGHBOURS, seeds(), &overrides);
    assert!(open.get(0, 1).is_some() && open.get(1, 0).is_some());

    let homology = Homology::from_pairs(
        [Pair {
            one: 0,
            other: 1,
            identity: 99.0,
            aligned_one: 90.0,
            aligned_other: 90.0,
        }],
        &lengths,
        HomologySettings::default(),
    );
    let severed = ContigFeatures::new(&coverage, &tnf, &lengths)
        .with_homology(Some(&homology))
        .graph_of(&indices, NEIGHBOURS, seeds(), &overrides);

    assert!(severed.get(0, 1).is_none() && severed.get(1, 0).is_none());
    assert_eq!(severed.nnz(), open.nnz() - 2);
}
