use rosella::markers::genome_tree::Tree;

const NEWICK: &str = "(((IMG_a:0.1,IMG_b:0.1)'UID5|g__Gen|1.000':0.1,(IMG_c:0.1,IMG_d:0.1)UID6||1.000:0.1)'UID2|k__Bacteria|1.000':0.2,(IMG_e:0.1,IMG_f:0.1)'UID3|k__Archaea|1.000':0.2)UID1||;";

const NODES: &str = "uid\tgenomes\ttaxonomy\tremoved\tgroups
UID1\t10\t\t\tR1,R2;R3
UID2\t6\tk__Bacteria\t\tB1,PF00002.7;B3,TIGR00398
UID3\t4\tk__Archaea\t\tA1
UID5\t3\tk__Bacteria;g__Gen\tPF00002\tPF00001.1;G2
UID6\t1\t\t\t
";

const SELECTED: &str = "uid\tselected
UID1\tUID1
UID2\tUID2
UID3\tUID3
UID4\tUID2
UID5\tUID4
UID6\tUID6
";

#[test]
fn walk_picks_the_set_lineage_wf_picks_and_searches_every_recorded_set() {
    let tree = Tree::new(NEWICK, NODES, SELECTED).unwrap();
    let cases = [
        (
            Some("IMG_a"),
            false,
            "k__Bacteria (UID2)",
            vec![vec!["B1"], vec!["B3"]],
            vec!["B1", "B3", "G2", "PF00001.1", "R1", "R2", "R3"],
        ),
        (
            Some("IMG_c"),
            false,
            "k__Bacteria (UID2)",
            vec![vec!["B1", "PF00002.7"], vec!["B3"]],
            vec!["B1", "B3", "PF00002.7", "R1", "R2", "R3"],
        ),
        (
            Some("IMG_a"),
            true,
            "k__Bacteria (UID2)",
            vec![vec!["B1"], vec!["B3"]],
            vec!["B1", "B3", "R1", "R2", "R3"],
        ),
        (
            Some("IMG_e"),
            true,
            "k__Archaea (UID3)",
            vec![vec!["A1"]],
            vec!["A1", "R1", "R2", "R3"],
        ),
        (
            None,
            false,
            "root (UID1)",
            vec![vec!["R1", "R2"], vec!["R3"]],
            vec!["R1", "R2", "R3"],
        ),
    ];
    for (leaf, forced, name, groups, searched) in cases {
        let chosen = tree
            .choose(leaf.map(|leaf| tree.leaf(leaf).unwrap()), forced)
            .unwrap();
        assert_eq!(chosen.name, name, "{leaf:?}");
        assert_eq!(chosen.groups, groups, "{leaf:?}");
        assert!(chosen.searched.iter().eq(searched.iter()), "{leaf:?}");
    }
}
