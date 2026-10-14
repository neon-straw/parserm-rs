//! Runs the binary on corpus files and compares with outputs produced by
//! parseRM.pl (tests/golden/caseN, regenerate with the oracle if the corpus
//! changes). Line order is ignored (Perl prints hash keys in random order),
//! average columns use relative-only 1e-9 tolerance; other fields are exact. Perl's
//! `REF(0x...)` addresses are normalised.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn manifest() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn normalise(text: &str) -> Vec<(String, Vec<String>)> {
    let re_free = |l: &str| -> String {
        match (l.find("REF(0x"), l.find(')')) {
            (Some(s), Some(e)) if e > s => format!("{}REF(ADDR){}", &l[..s], &l[e + 1..]),
            _ => l.to_string(),
        }
    };
    let mut header = Vec::new();
    let mut v = Vec::new();
    for line in text.lines() {
        let is_header = line.starts_with('#');
        if is_header && line.contains('\t') {
            header = line.trim_start_matches('#').split('\t').map(String::from).collect();
        }
        v.push((re_free(line), if is_header { Vec::new() } else { header.clone() }));
    }
    v.sort();
    v
}

fn fields_match(a: &str, b: &str, column: &str) -> bool {
    if a == b {
        return true;
    }
    if !matches!(column, "AVG_%DIV" | "AVG_%DEL" | "AVG_%INS" | "AVG_LEN_MASKED" | "%_GENOME" | "%_masked") {
        return false;
    }
    match (a.parse::<f64>(), b.parse::<f64>()) {
        (Ok(x), Ok(y)) => x == y || (x - y).abs() / x.abs().max(y.abs()) <= 1e-9,
        _ => false,
    }
}

fn check_case(dir: &Path, input: &str, args: &str) {
    let tmp = std::env::temp_dir().join(format!(
        "parserm-golden-{}-{}",
        std::process::id(),
        dir.file_name().unwrap().to_string_lossy()
    ));
    let _ = fs::remove_dir_all(&tmp);
    fs::create_dir_all(&tmp).unwrap();
    let corpus = manifest().join("tests/corpus");
    for e in fs::read_dir(&corpus).unwrap() {
        let e = e.unwrap();
        fs::copy(e.path(), tmp.join(e.file_name())).unwrap();
    }
    let status = Command::new(env!("CARGO_BIN_EXE_parserm-rs"))
        .current_dir(&tmp)
        .arg("-i")
        .arg(input)
        .args(args.split_whitespace())
        .stderr(std::process::Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "{input} {args}: exit {status}");
    for e in fs::read_dir(dir).unwrap() {
        let e = e.unwrap();
        let name = e.file_name();
        let want = normalise(&fs::read_to_string(e.path()).unwrap());
        let got = normalise(
            &fs::read_to_string(tmp.join(&name))
                .unwrap_or_else(|_| panic!("{input} {args}: missing {}", name.to_string_lossy())),
        );
        assert_eq!(want.len(), got.len(), "{input} {args}: {} line count", name.to_string_lossy());
        for ((w, columns), (g, _)) in want.iter().zip(&got) {
            let ok = w.split('\t').count() == g.split('\t').count()
                && w.split('\t').zip(g.split('\t')).enumerate().all(|(i, (x, y))| fields_match(x, y, columns.get(i).map_or("", String::as_str)));
            assert!(ok, "{input} {args}: {}\n want: {w}\n  got: {g}", name.to_string_lossy());
        }
    }
    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn golden_outputs_match_oracle() {
    let golden = manifest().join("tests/golden");
    let cases = fs::read_to_string(golden.join("cases.txt")).unwrap();
    for (i, line) in cases.lines().filter(|l| !l.trim().is_empty()).enumerate() {
        let (input, args) = line.split_once('|').unwrap();
        check_case(&golden.join(format!("case{}", i + 1)), input, args);
    }
}

#[test]
fn edit_option_dies_like_oracle() {
    let tmp = std::env::temp_dir().join(format!("parserm-edit-{}", std::process::id()));
    fs::create_dir_all(&tmp).unwrap();
    let corpus = manifest().join("tests/corpus");
    for f in ["nesting.out", "edit.txt"] {
        fs::copy(corpus.join(f), tmp.join(f)).unwrap();
    }
    let out = Command::new(env!("CARGO_BIN_EXE_parserm-rs"))
        .current_dir(&tmp)
        .args(["-i", "nesting.out", "-p", "-edit", "edit.txt"])
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(!tmp.join("nesting.out.parseRM.summary.tab").exists());
    let fixed = Command::new(env!("CARGO_BIN_EXE_parserm-rs"))
        .current_dir(&tmp)
        .args(["-i", "nesting.out", "-p", "-edit", "edit.txt", "--fix-quirks"])
        .output()
        .unwrap();
    assert!(fixed.status.success());
    assert!(tmp.join("nesting.out.parseRM.summary.tab").exists());
    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn directory_named_like_input_is_skipped() {
    // Perl opens the directory, readline returns undef: no outputs, exit 0.
    let tmp = std::env::temp_dir().join(format!("parserm-dirin-{}", std::process::id()));
    fs::create_dir_all(tmp.join("x.out")).unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_parserm-rs"))
        .current_dir(&tmp)
        .args(["-i", "x.out", "-p"])
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(fs::read_dir(&tmp).unwrap().count(), 1);
    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn median_numeric_ties_preserve_position_order() {
    let tmp = std::env::temp_dir().join(format!("parserm-median-{}", std::process::id()));
    fs::create_dir_all(&tmp).unwrap();
    fs::copy(
        manifest().join("tests/corpus/partial.out"),
        tmp.join("partial.out"),
    )
    .unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_parserm-rs"))
        .current_dir(&tmp)
        .args(["-i", "partial.out", "-p"])
        .status()
        .unwrap();
    assert!(status.success());
    let output = fs::read_to_string(tmp.join("partial.out.parseRM.all-repeats.tab")).unwrap();
    for name in ["MedianTie", "MedianLex"] {
        let row = output
            .lines()
            .find(|row| row.starts_with(&format!("{name}\t")))
            .unwrap();
        let fields: Vec<_> = row.split('\t').collect();
        // This binary path loads hits, merges histogram weights, then applies
        // tied_median; MedianTie has input strings [5.0, 5, 5.0].
        assert_eq!(fields[6], if name == "MedianTie" { "3" } else { "25" });
        for column in [8, 10, 12] {
            assert_eq!(fields[column], "5", "{name} median column {column}");
        }
    }
    fs::remove_dir_all(tmp).unwrap();
}

#[test]
fn average_fields_use_relative_tolerance_without_floor() {
    assert!(!fields_match("0", "9e-13", "AVG_%DIV"));
    assert!(!fields_match("1e-15", "9e-13", "AVG_%DIV"));
    assert!(fields_match("1e-15", "1.0000000001e-15", "AVG_%DIV"));
    assert!(fields_match("0", "0", "AVG_%DIV"));
}

#[test]
fn median_fields_require_exact_strings() {
    assert!(!fields_match("5", "5.0", "MED_%DIV"));
    assert!(fields_match("5", "5", "MED_%DIV"));
}

#[test]
fn non_average_fields_require_exact_strings() {
    assert!(!fields_match("5", "5.0", "LEN_MASKED_NR"));
}
