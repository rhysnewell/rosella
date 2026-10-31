use rosella::markers::{ContigMarkers, Hit, MarkerSet};
use rosella::quality::Scorer;

const TABLE: &str = "model_name\tsets\n\
                     alpha\tbac\n\
                     beta\tbac\n\
                     gamma\tbac\n\
                     delta\tbac\n";

const SETS: &str = "set\tmedian_genome_bp\tmax_genome_bp\n\
                    bac\t300000\t0\n";

fn hit(marker: u16) -> Hit {
    Hit {
        marker,
        ..Default::default()
    }
}

#[test]
fn a_marker_doubled_alone_in_many_bins_stops_counting_as_contamination() {
    let mut per_contig = (0..60u16)
        .map(|bin| {
            let mut hits = vec![hit(0), hit(1), hit(2), hit(3)];
            if bin % 2 == 0 {
                hits.push(hit(0));
            }
            hits
        })
        .collect::<Vec<_>>();
    per_contig[1].push(hit(1));
    let lengths = vec![1_000_000; per_contig.len()];
    let markers = ContigMarkers::new(per_contig, MarkerSet::parse(TABLE).with_bounds(SETS))
        .with_lengths(lengths);
    let bins = (0..60).map(|contig| vec![contig]).collect::<Vec<_>>();

    markers.learn_multicopy(bins.iter().map(Vec::as_slice));

    assert_eq!(markers.score(&[0]).contamination, 0.0);
    assert!(markers.score(&[1]).contamination > 0.0);
}
