use rosella::markers::checkm::{Cutoffs, Panel};

const SETS: &str = "set\tpanel\tmodel\tgroup\n\
                    bac\tcheckm\tPF00001.1\t0\n\
                    bac\tcheckm\tPF00002.1\t0\n\
                    bac\tcheckm\tTIGR00003\t1\n";
const CLANS: &str = "pfam\tclan\tnested\nPF00001\tCL0001\t\nPF00002\tCL0001\t\n";
const HMM: &str = "ACC   PF00001.1\nGA    20.0 20.0;\n//\n\
                   ACC   PF00002.1\nGA    20.0 20.0;\n//\n\
                   ACC   TIGR00003\nNC    20.0 20.0;\n//\n";

fn row(protein: usize, model: &str, e_value: f64, score: f64, from: u32, to: u32) -> String {
    format!(
        "{protein} - 300 name {model} 100 {e_value} {score} 0.0 1 1 {e_value} {e_value} {score} 0.0 \
         1 90 {from} {to} {from} {to} 0.9 -\n"
    )
}

fn copies(table: &str) -> Vec<Vec<(u16, u16)>> {
    let panel = Panel::parse(SETS, CLANS, |_| None);
    let contig_of = |protein: usize| Some(usize::from(protein >= 10));
    panel
        .tally(table, &Cutoffs::parse(HMM), contig_of, 2)
        .into_iter()
        .map(|held| {
            held.into_iter()
                .map(|entry| (entry.model, entry.copies))
                .collect()
        })
        .collect()
}

#[test]
fn a_gene_split_across_neighbouring_calls_is_one_copy_and_a_distant_one_is_two() {
    let tigr = |protein| row(protein, "TIGR00003", 1e-30, 80.0, 1, 90);
    let split = format!("{}{}", tigr(3), tigr(4));
    let distant = format!("{}{}", tigr(3), tigr(6));
    let three = format!("{}{}{}", tigr(3), tigr(4), tigr(5));

    assert_eq!(copies(&split)[0], vec![(2, 1)]);
    assert_eq!(copies(&distant)[0], vec![(2, 2)]);
    assert_eq!(copies(&three)[0], vec![(2, 2)]);
}

#[test]
fn two_pfams_of_one_clan_on_one_stretch_keep_the_better() {
    let clash = format!(
        "{}{}",
        row(3, "PF00001.1", 1e-20, 60.0, 10, 80),
        row(3, "PF00002.1", 1e-40, 90.0, 40, 120),
    );
    let apart = format!(
        "{}{}",
        row(3, "PF00001.1", 1e-20, 60.0, 10, 80),
        row(3, "PF00002.1", 1e-40, 90.0, 80, 160),
    );

    assert_eq!(copies(&clash)[0], vec![(1, 1)]);
    assert_eq!(copies(&apart)[0], vec![(0, 1), (1, 1)]);
}

#[test]
fn a_hit_under_its_model_cutoff_or_aligned_over_too_little_is_no_hit() {
    let weak = row(3, "PF00001.1", 1e-3, 15.0, 10, 80);
    let sliver = row(3, "PF00001.1", 1e-30, 60.0, 10, 30);

    assert!(copies(&weak)[0].is_empty());
    assert!(copies(&sliver)[0].is_empty());
}

#[test]
fn a_group_counts_once_however_many_markers_it_holds() {
    let panel = Panel::parse(SETS, CLANS, |_| None);
    let held = [2u32, 0, 1];
    let (completeness, contamination) = panel.score(0, |model| held[model as usize]);

    assert!((completeness - 75.0).abs() < 1e-9);
    assert!((contamination - 25.0).abs() < 1e-9);
}
