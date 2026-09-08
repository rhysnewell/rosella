//! The feature vector the trained models were fitted on, built from gene counts.

use rosella::quality::tables::{METADATA, Tables};

/// Pathways and categories count a gene once however many copies are seen; modules sum the
/// raw counts and divide by the definition length before any missing gene is dropped. Getting
/// either wrong shifts 1220 columns and the models read them all.
#[test]
fn the_derived_columns_follow_their_own_arithmetic() {
    let tables = Tables::load().expect("embedded tables");
    assert_eq!(tables.gene_count, 19999);
    assert_eq!(
        tables.width(),
        METADATA + 19999 + tables.pathways.len() + tables.modules.len() + tables.categories.len()
    );

    let mut genes = vec![0.0; tables.gene_count];
    let module = tables
        .modules
        .iter()
        .find(|group| group.columns.len() > 2)
        .expect("a module with members");
    let length = module.length;
    let columns = module.columns.clone();
    for column in &columns {
        genes[*column as usize] = 3.0;
    }

    let mut vector = Vec::new();
    tables.fill(&[0.0; METADATA], &genes, &mut vector);
    assert_eq!(vector.len(), tables.width());

    let modules_at = METADATA + tables.gene_count + tables.pathways.len();
    let at = tables
        .modules
        .iter()
        .position(|group| group.columns == columns)
        .expect("the module back again");
    let expected = 3.0 * columns.len() as f64 / length as f64;
    assert!(
        (vector[modules_at + at] - expected).abs() < 1e-12,
        "modules sum raw counts: {} against {expected}",
        vector[modules_at + at]
    );

    let pathway = tables
        .pathways
        .iter()
        .enumerate()
        .find(|(_, group)| !group.columns.is_empty())
        .map(|(at, group)| (at, group.columns.clone(), group.length));
    let Some((at, columns, length)) = pathway else {
        panic!("no pathway holds a gene");
    };
    let mut genes = vec![0.0; tables.gene_count];
    for column in &columns {
        genes[*column as usize] = 7.0;
    }
    tables.fill(&[0.0; METADATA], &genes, &mut vector);
    let expected = columns.len() as f64 / length as f64;
    let value = vector[METADATA + tables.gene_count + at];
    assert!(
        (value - expected).abs() < 1e-12,
        "pathways clip to presence: {value} against {expected}"
    );
}

/// Metadata leads the vector and the gene counts follow it, because the models index columns
/// by position and nothing downstream would notice the two being swapped.
#[test]
fn metadata_leads_the_vector() {
    let tables = Tables::load().expect("embedded tables");
    let mut metadata = [0.0; METADATA];
    metadata[21] = 42.0;
    let mut genes = vec![0.0; tables.gene_count];
    genes[0] = 5.0;

    let mut vector = Vec::new();
    tables.fill(&metadata, &genes, &mut vector);
    assert_eq!(vector[21], 42.0);
    assert_eq!(vector[METADATA], 5.0);
}
