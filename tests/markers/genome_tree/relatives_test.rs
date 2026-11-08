use rosella::markers::genome_tree::Relatives;

const MODELS: &str = "model\tsearched\tsequence\tdomain
PF00376.18\tPF00376.18\t20\t20
PF13411.1\tPF13411.1\t20\t20
PF00583.19\tPF00583.19\t20\t20
TIGR00001\tTIGR00001\t20\t20
";

const CLANS: &str = "pfam\tclan\tnested
PF00376\tCL0123\t
PF13411\tCL0123\t
PF00583\tCL0257\t
";

#[test]
fn a_chosen_model_brings_only_its_own_clan() {
    let relatives = Relatives::parse(MODELS, CLANS);
    let chosen = ["PF00376.18".to_string(), "TIGR00001".to_string()];
    let widened = relatives.widened(chosen.iter());
    assert_eq!(
        widened.into_iter().collect::<Vec<_>>(),
        ["PF00376.18", "PF13411.1", "TIGR00001"]
    );
}
