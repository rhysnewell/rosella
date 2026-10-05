use rosella::kmers::embedding::Embedding;

fn write(dir: &std::path::Path, file: &str, values: &[f32]) {
    let bytes = values
        .iter()
        .flat_map(|value| value.to_le_bytes())
        .collect::<Vec<_>>();
    std::fs::write(dir.join(file), bytes).unwrap();
}

#[test]
fn rows_and_halves_come_back_by_name_in_the_order_asked() {
    let home = tempfile::tempdir().unwrap();
    std::fs::write(home.path().join("whole.names"), "a\nb\n").unwrap();
    write(home.path(), "whole.f32", &[1.0, 2.0, 3.0, 4.0]);
    std::fs::write(home.path().join("halves.names"), "b\n").unwrap();
    write(home.path(), "first.f32", &[5.0, 6.0]);
    write(home.path(), "second.f32", &[7.0, 8.0]);

    let embedding = Embedding::read(home.path()).unwrap();

    let rows = embedding.rows(["b", "a"].into_iter()).unwrap();
    assert_eq!(rows.row(0).to_vec(), vec![3.0, 4.0]);
    assert_eq!(rows.row(1).to_vec(), vec![1.0, 2.0]);
    let [first, second] = embedding.halves(&["b"]).unwrap();
    assert_eq!(first.row(0).to_vec(), vec![5.0, 6.0]);
    assert_eq!(second.row(0).to_vec(), vec![7.0, 8.0]);
    assert!(embedding.halves(&["a"]).is_err());
}
