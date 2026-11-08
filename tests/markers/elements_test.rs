use rosella::markers::elements::Families;

const TABLE: &str = "family\tplasmid\tchromosome\tvirus\tplasmid_chromosome
Terminase\t0.0001\t0.0001\t0.2\t0.0001
Rep\t0.001\t0.00001\t0.0001\t0.001
Mob\t0.001\t0.00001\t0.0001\t0.001
Tra\t0.001\t0.00001\t0.0001\t0.001
Host\t0.001\t0.011\t0.0001\t0.01
";

#[test]
fn families_count_once_and_a_host_family_vetoes_a_plasmid() {
    let families = Families::parse(TABLE).unwrap();
    let cases = [
        (vec!["Terminase"], true),
        (vec!["Rep", "Rep"], false),
        (vec!["Rep", "Mob", "Tra"], true),
        (vec!["Rep", "Mob", "Tra", "Host"], false),
        (vec!["Unknown"], false),
    ];
    for (names, element) in cases {
        assert_eq!(families.call(names.clone()).element(), element, "{names:?}");
    }
}
