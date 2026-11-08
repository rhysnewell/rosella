use rosella::markers::elements::Families;

const TABLE: &str = "family\tplasmid\tchromosome\tvirus\tplasmid_chromosome
Terminase\t0.0001\t0.0001\t0.2\t0.0001
Rep\t0.001\t0.00001\t0.0001\t0.001
Mob\t0.001\t0.00001\t0.0001\t0.001
";

#[test]
fn a_family_counts_once_toward_the_bar() {
    let families = Families::parse(TABLE).unwrap();
    let cases = [
        (vec!["Terminase"], true),
        (vec!["Rep", "Rep"], false),
        (vec!["Rep", "Mob"], true),
        (vec!["Unknown"], false),
    ];
    for (names, element) in cases {
        assert_eq!(families.call(names.clone()).element(), element, "{names:?}");
    }
}
