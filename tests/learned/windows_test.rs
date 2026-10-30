use rosella::coverage::depth_file::{self, Depth};
use rosella::learned::pieces::{VIEWS, view_key, views};
use rosella::learned::windows::Windows;

const WINDOW: usize = 100;

fn half(value: u16) -> u16 {
    if value == 0 {
        return 0;
    }
    let exponent = 15 - value.leading_zeros() as u16;
    let mantissa = (u32::from(value) << 10 >> exponent) as u16 & 0x3ff;
    ((exponent + 15) << 10) | mantissa
}

#[test]
fn each_view_reads_the_mean_depth_of_its_own_span() {
    let directory = tempfile::tempdir().unwrap();
    let lengths: [usize; 2] = [1234, 2000];
    let names = ["first".to_string(), "second".to_string()];
    let windows = lengths
        .iter()
        .map(|length| length.div_ceil(WINDOW))
        .collect::<Vec<_>>();
    let depth = |contig: usize, window: usize| (3 * contig + window % 7 + 1) as u16;
    let rows = names
        .iter()
        .zip(lengths)
        .map(|(name, length)| format!("{name}\t{length}\n"))
        .collect::<String>();
    std::fs::write(
        directory.path().join("contigs.tsv"),
        format!("contig\tlength\n{rows}"),
    )
    .unwrap();
    let bytes = windows
        .iter()
        .enumerate()
        .flat_map(|(contig, count)| (0..*count).map(move |window| half(depth(contig, window))))
        .flat_map(u16::to_le_bytes)
        .collect::<Vec<_>>();
    let none = vec![0.0; names.len()];
    let table = Depth {
        means: none.clone(),
        variances: none,
        reads: vec![0; names.len()],
    };
    depth_file::write(&depth_file::path(directory.path(), "s0"), &table, &bytes).unwrap();
    let spans = Windows::read(directory.path(), &["s0".to_string()])
        .unwrap()
        .spans(&names, &lengths, 42)
        .unwrap();
    for (contig, (name, length)) in names.iter().zip(lengths).enumerate() {
        let views = views(view_key(42, name), length);
        assert_eq!(views.len(), VIEWS);
        for (view, (start, end)) in views.into_iter().enumerate() {
            let sum = (start..end)
                .map(|base| f64::from(depth(contig, base / WINDOW)))
                .sum::<f64>();
            let expected = sum / (end - start) as f64;
            let given = f64::from(spans.view(contig, view)[0]);
            assert!(
                (given - expected).abs() < 1e-4 * expected,
                "{name} view {view}: {given} against {expected}"
            );
        }
    }
}
