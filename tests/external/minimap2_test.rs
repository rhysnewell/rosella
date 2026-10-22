use std::ffi::OsString;
use std::io::Read;
use std::path::Path;

use rosella::external::minimap2::minimap2;

const DATA: &str = "tests/data/minimap2";

fn words(words: &[&str]) -> Vec<OsString> {
    words.iter().map(OsString::from).collect()
}

fn reads() -> [OsString; 2] {
    ["reads_1.fq.gz", "reads_2.fq.gz"].map(|name| Path::new(DATA).join(name).into())
}

fn sam(args: &[OsString]) -> String {
    minimap2(args, |sam| {
        let mut text = String::new();
        sam.read_to_string(&mut text)?;
        Ok(text)
    })
    .unwrap()
}

#[test]
fn the_library_writes_the_command_lines_sam_for_a_whole_and_a_split_index() {
    let home = tempfile::tempdir().unwrap();
    for (case, split) in [("whole", None), ("split", Some("5k"))] {
        let index = home.path().join(format!("{case}.mmi"));
        let mut build = words(&["-x", "sr", "-t", "2"]);
        if let Some(split) = split {
            build.extend(words(&["-I", split]));
        }
        build.extend(["-d".into(), index.clone().into()]);
        build.push(Path::new(DATA).join("assembly.fa").into());
        assert_eq!(sam(&build), "", "{case}");

        let mut map = words(&["-a", "--secondary=yes", "-x", "sr", "-t", "2"]);
        if split.is_some() {
            map.extend(["--split-prefix".into(), home.path().join(case).into()]);
        }
        map.push(index.into());
        map.extend(reads());
        let body = sam(&map);
        let got: Vec<&str> = body
            .lines()
            .filter(|line| !line.starts_with("@PG"))
            .collect();
        let expected = std::fs::read_to_string(Path::new(DATA).join(format!("{case}.sam")));
        let expected = expected.unwrap();
        let want: Vec<&str> = expected.lines().collect();
        for (row, (got, want)) in got.iter().zip(&want).enumerate() {
            assert_eq!(got, want, "{case} row {row}");
        }
        assert_eq!(got.len(), want.len(), "{case}");
    }
}

#[test]
fn a_failed_run_or_a_stopped_reader_is_an_error_and_not_an_exit() {
    let home = tempfile::tempdir().unwrap();
    let mut missing = words(&["-a", "-x", "sr"]);
    missing.push(home.path().join("absent.mmi").into());
    missing.extend(reads());
    let failed = minimap2(&missing, |_| Ok(())).unwrap_err();
    assert!(failed.to_string().contains("status 1"), "{failed}");

    let mut map = words(&["-a", "--secondary=yes", "-x", "sr"]);
    map.push(Path::new(DATA).join("assembly.fa").into());
    map.extend(reads());
    let stopped = minimap2(&map, |_| -> anyhow::Result<()> {
        anyhow::bail!("full disk")
    });
    assert_eq!(stopped.unwrap_err().to_string(), "full disk");
}
