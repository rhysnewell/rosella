use rosella::markers::viruses::{Families, viral};

const TABLE: &str = "family\tplasmid\tchromosome\tvirus\tplasmid_chromosome
Terminase\t0.0001\t0.0001\t0.2\t0.0001
Tail\t0.0001\t0.0001\t0.01\t0.0001
Portal\t0.0001\t0.0001\t0.01\t0.0001
Rep\t0.05\t0.00001\t0.0001\t0.002
";

#[test]
fn virus_families_count_once_and_plasmid_genes_count_against() {
    let families = Families::parse(TABLE).unwrap();
    let cases = [
        (vec!["Terminase"], true),
        (vec!["Tail", "Tail"], false),
        (vec!["Tail", "Portal"], true),
        (vec!["Rep"], false),
        (vec!["Rep", "Terminase"], false),
    ];
    for (names, called) in cases {
        assert_eq!(viral(families.virus(names.clone())), called, "{names:?}");
    }
}
