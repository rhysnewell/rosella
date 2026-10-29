use rosella::coverage::depth_file::{self, Depth};
use rosella::learned::pieces::{VIEWS, view_key, views};
use rosella::learned::windows::{SpanDepth, Windows};

const WINDOW: usize = 100;

fn half(value: u16) -> u16 {
    if value == 0 {
        return 0;
    }
    let exponent = 15 - value.leading_zeros() as u16;
    let mantissa = (u32::from(value) << 10 >> exponent) as u16 & 0x3ff;
    ((exponent + 15) << 10) | mantissa
}

const LENGTHS: [usize; 2] = [1234, 2000];

fn depth(contig: usize, window: usize) -> u16 {
    (3 * contig + window % 7 + 1) as u16
}

fn span_mean(contig: usize, (start, end): (usize, usize)) -> f64 {
    let sum = (start..end)
        .map(|base| f64::from(depth(contig, base / WINDOW)))
        .sum::<f64>();
    sum / (end - start) as f64
}

fn written() -> (tempfile::TempDir, Vec<String>, Windows) {
    let directory = tempfile::tempdir().unwrap();
    let names = vec!["first".to_string(), "second".to_string()];
    let counts = LENGTHS
        .iter()
        .map(|length| length.div_ceil(WINDOW))
        .collect::<Vec<_>>();
    let rows = names
        .iter()
        .zip(LENGTHS)
        .map(|(name, length)| format!("{name}\t{length}\n"))
        .collect::<String>();
    std::fs::write(
        directory.path().join("contigs.tsv"),
        format!("contig\tlength\n{rows}"),
    )
    .unwrap();
    let bytes = counts
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
    let windows = Windows::read(directory.path(), &["s0".to_string()]).unwrap();
    (directory, names, windows)
}

fn check(spans: &SpanDepth, cuts: &[Vec<(usize, usize)>]) {
    for (contig, cuts) in cuts.iter().enumerate() {
        for (view, span) in cuts.iter().enumerate() {
            let expected = span_mean(contig, *span);
            let given = f64::from(spans.view(contig, view)[0]);
            assert!(
                (given - expected).abs() < 1e-4 * expected,
                "contig {contig} view {view}: {given} against {expected}"
            );
        }
    }
}

#[test]
fn each_view_reads_the_mean_depth_of_its_own_span() {
    let (_directory, names, windows) = written();
    let spans = windows.spans(&names, &LENGTHS, 42).unwrap();
    let cuts = names
        .iter()
        .zip(LENGTHS)
        .map(|(name, length)| views(view_key(42, name), length))
        .collect::<Vec<_>>();
    assert!(cuts.iter().all(|cuts| cuts.len() == VIEWS));
    check(&spans, &cuts);
}

#[test]
fn a_contig_finds_its_spans_after_a_contig_cut_fewer_times() {
    let (_directory, names, windows) = written();
    let cuts = [
        vec![
            (0, LENGTHS[0]),
            (0, LENGTHS[0] / 2),
            (LENGTHS[0] / 2, LENGTHS[0]),
        ],
        (0..8)
            .map(|piece| (250 * piece, 250 * (piece + 1)))
            .collect(),
    ];
    let spans = windows.cut(&names, |at| cuts[at].clone()).unwrap();
    check(&spans, &cuts);
}
