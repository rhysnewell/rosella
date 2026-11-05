use rosella::markers::checkm::{Copies, Panel};
use rosella::markers::replicon::Shape;
use rosella::markers::{ContigMarkers, Hit, MarkerSet};
use rosella::quality::{Bars, Scorer};

const BARS: Bars = Bars {
    completeness: 90.0,
    contamination: 5.0,
};

struct Bac {
    gtdb: u16,
    copy: Copies,
}

fn bac() -> Bac {
    let set = MarkerSet::embedded();
    let panel = Panel::embedded(|model| set.id(model));
    let at = panel.lineage("bac").unwrap();
    Bac {
        gtdb: set.id("ATPsyn_F1gamma").unwrap(),
        copy: Copies {
            set: at as u8,
            model: panel.models(at)[0],
            score: 100.0,
            partial: false,
        },
    }
}

fn scored(copy: Copies, score: f32) -> Copies {
    Copies { score, ..copy }
}

fn hit(marker: u16) -> Hit {
    Hit {
        marker,
        ..Default::default()
    }
}

fn whole() -> Vec<Hit> {
    (0..MarkerSet::embedded().len() as u16).map(hit).collect()
}

fn markers(hits: Vec<Vec<Hit>>, checkm: Vec<Vec<Copies>>) -> ContigMarkers {
    let lengths = vec![900_000; hits.len()];
    ContigMarkers::new(hits, MarkerSet::embedded())
        .with_lengths(lengths)
        .with_checkm(checkm)
}

fn coding(held: ContigMarkers, coding_bases: &[u64]) -> ContigMarkers {
    held.with_shapes(
        coding_bases
            .iter()
            .map(|coding_bases| {
                Some(Shape {
                    coding_bases: *coding_bases,
                    genes: 1,
                })
            })
            .collect(),
    )
}

#[test]
fn checkm_reads_once_learned_and_the_higher_panel_sets_the_contamination() {
    let Bac { gtdb, copy } = bac();
    let checkm_doubled = markers(
        vec![vec![hit(gtdb)], Vec::new()],
        vec![vec![copy], vec![copy]],
    );
    let gtdb_doubled = markers(
        vec![vec![hit(gtdb)], vec![hit(gtdb)]],
        vec![vec![copy], vec![]],
    );

    assert_eq!(checkm_doubled.score(&[0, 1]).contamination, 0.0);
    assert!(checkm_doubled.redundant(&[0, 1]).is_empty());
    let gtdb_alone = gtdb_doubled.score(&[0, 1]).contamination;
    for held in [&checkm_doubled, &gtdb_doubled] {
        held.learn_checkm([[0, 1].as_slice()], BARS);
    }

    assert!(checkm_doubled.score(&[0, 1]).contamination > 0.0);
    assert_eq!(checkm_doubled.redundant(&[0, 1]), vec![1]);
    assert_eq!(gtdb_doubled.score(&[0, 1]).contamination, gtdb_alone);
}

#[test]
fn a_bin_coding_like_a_eukaryote_keeps_the_gtdb_contamination() {
    let Bac { gtdb, copy } = bac();
    let mut hits = vec![whole(); 4];
    hits.extend([vec![hit(gtdb)], Vec::new(), vec![hit(gtdb)], Vec::new()]);
    let mut checkm = vec![Vec::new(); 4];
    checkm.extend([vec![copy], vec![copy], vec![copy], vec![copy]]);
    let dense = 810_000;
    let held = coding(
        markers(hits, checkm),
        &[dense, dense, dense, dense, 495_000, 495_000, dense, dense],
    );
    held.learn_checkm([[0], [1], [2], [3]].iter().map(|bin| bin.as_slice()), BARS);

    assert_eq!(held.score(&[4, 5]).contamination, 0.0);
    assert!(held.score(&[6, 7]).contamination > 0.0);
}

#[test]
fn copies_of_a_model_on_one_contig_count_once() {
    let Bac { gtdb, copy } = bac();
    let held = markers(
        vec![vec![hit(gtdb)], Vec::new()],
        vec![vec![copy, scored(copy, 90.0)], Vec::new()],
    );
    held.learn_checkm([[0, 1].as_slice()], BARS);

    let (chosen, sets) = held.readings(&[0, 1]).unwrap();
    assert!(sets[chosen].checkm.unwrap().contamination > 0.0);
    assert_eq!(held.score(&[0, 1]).contamination, 0.0);
}

#[test]
fn a_copy_far_under_the_bins_best_copy_is_a_paralog() {
    let Bac { gtdb, copy } = bac();
    let ratios = [0.15, 0.17, 0.95, 0.98];
    let mut hits = Vec::new();
    let mut checkm = Vec::new();
    for ratio in ratios {
        hits.extend([whole(), Vec::new()]);
        checkm.extend([vec![copy], vec![scored(copy, 100.0 * ratio)]]);
    }
    for second in [15.0, 96.0] {
        hits.extend([vec![hit(gtdb)], Vec::new()]);
        checkm.extend([vec![copy], vec![scored(copy, second)]]);
    }
    let held = markers(hits, checkm);
    let learned = (0..ratios.len())
        .map(|bin| [2 * bin, 2 * bin + 1])
        .collect::<Vec<_>>();
    held.learn_checkm(learned.iter().map(|bin| bin.as_slice()), BARS);

    assert_eq!(held.score(&[8, 9]).contamination, 0.0);
    assert!(held.score(&[10, 11]).contamination > 0.0);
}

#[test]
fn a_checkm_model_doubled_in_every_bin_counts_once() {
    let Bac { copy, .. } = bac();
    let set = MarkerSet::embedded();
    let panel = Panel::embedded(|model| set.id(model));
    let every = panel
        .models(usize::from(copy.set))
        .iter()
        .map(|model| Copies {
            model: *model,
            ..copy
        })
        .collect::<Vec<_>>();
    let bins = 8;
    let mut hits = Vec::new();
    let mut checkm = Vec::new();
    for _ in 0..bins {
        hits.extend([whole(), Vec::new()]);
        checkm.extend([every.clone(), vec![copy]]);
    }
    let held = markers(hits, checkm);
    let pairs = (0..bins)
        .map(|bin| [2 * bin, 2 * bin + 1])
        .collect::<Vec<_>>();
    held.learn_checkm(pairs.iter().map(|bin| bin.as_slice()), BARS);
    let gtdb = held.gtdb_score(&[0, 1]).contamination;

    assert!(held.score(&[0, 1]).contamination > gtdb);
    held.learn_multicopy(pairs.iter().map(|bin| bin.as_slice()));
    assert_eq!(held.score(&[0, 1]).contamination, gtdb);
    assert!(held.redundant(&[0, 1]).is_empty());
}
