use rosella::external::hmmer_engine::HmmerEngine;

fn dealt(shards: usize, records: usize) -> Vec<String> {
    let directory = tempfile::tempdir().unwrap();
    let engine = HmmerEngine::new(shards * 2, Some(shards));
    let mut sink = engine.protein_shards(directory.path()).unwrap();
    for id in 0..records {
        sink.write(id, "MKRLGSEVCDAIREF").unwrap();
    }
    sink.finish()
        .unwrap()
        .iter()
        .map(|path| std::fs::read_to_string(path).unwrap())
        .collect()
}

/// The reader that this replaced advanced its counter before writing, so record zero belongs in
/// shard one. Anything else reshuffles the search and stops the ledger comparing.
#[test]
fn a_record_lands_one_shard_past_its_index() {
    for shards in [1usize, 3, 4] {
        let pieces = dealt(shards, 10);
        assert_eq!(pieces.len(), shards);
        for id in 0..10 {
            let wanted = (id + 1) % shards;
            assert!(
                pieces[wanted].contains(&format!(">{id}\nMKRLGSEVCDAIREF\n")),
                "{shards} shards, record {id} is not in shard {wanted}"
            );
        }
        assert_eq!(
            pieces
                .iter()
                .map(|piece| piece.lines().count())
                .sum::<usize>(),
            20
        );
    }
}
