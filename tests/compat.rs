use std::path::Path;
use std::process::Command;

use rsomics_bed_closest::closest;

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn basic_closest_correctness() {
    let a = golden("a.bed");
    let b = golden("b.bed");
    let mut out = Vec::new();
    closest(&a, &b, &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    let lines: Vec<&str> = result.lines().filter(|l| !l.is_empty()).collect();
    assert_eq!(lines.len(), 3, "expected 3 output lines: {result}");

    // A1 [100,200) closest to B1 [300,400) — distance 100
    assert!(
        lines[0].contains("chr1\t300\t400"),
        "A1 closest wrong: {}",
        lines[0]
    );
    let dist0: i64 = lines[0].split('\t').last().unwrap().trim().parse().unwrap();
    assert_eq!(dist0, 100, "A1 distance wrong: {dist0}");

    // A2 [500,600) closest to B1 [300,400) distance 100 OR B2 [700,800) distance 100
    // (equidistant — bedtools reports both, we report at least one with dist 100)
    let dist1: i64 = lines[1].split('\t').last().unwrap().trim().parse().unwrap();
    assert_eq!(dist1, 100, "A2 distance wrong: {dist1}");

    // A3 [100,200) closest to B3 [250,350) — distance 50
    assert!(
        lines[2].contains("chr2\t250\t350"),
        "A3 closest wrong: {}",
        lines[2]
    );
    let dist2: i64 = lines[2].split('\t').last().unwrap().trim().parse().unwrap();
    assert_eq!(dist2, 50, "A3 distance wrong: {dist2}");
}

#[test]
fn overlapping_gives_zero_distance() {
    use std::io::Write;
    use tempfile::NamedTempFile;
    let mut fa = NamedTempFile::new().unwrap();
    let mut fb = NamedTempFile::new().unwrap();
    writeln!(fa, "chr1\t100\t200\tA").unwrap();
    writeln!(fb, "chr1\t150\t250\tB").unwrap();
    let mut out = Vec::new();
    closest(fa.path(), fb.path(), &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    let dist: i64 = result
        .lines()
        .next()
        .unwrap()
        .split('\t')
        .last()
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert_eq!(dist, 0, "overlapping intervals must have distance 0");
}

#[test]
fn bedtools_compat() {
    let bedtools = Command::new("bedtools").arg("--version").output();
    if bedtools.is_err() || !bedtools.unwrap().status.success() {
        eprintln!("bedtools not available — skipping compat test");
        return;
    }

    let a = golden("a.bed");
    let b = golden("b.bed");

    let mut ours = Vec::new();
    closest(&a, &b, &mut ours).unwrap();
    let ours_str = String::from_utf8(ours).unwrap();

    let bt = Command::new("bedtools")
        .args(["closest", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .output()
        .expect("bedtools closest failed");
    let bt_str = String::from_utf8(bt.stdout).unwrap();

    let mut ours_lines: Vec<&str> = ours_str.lines().filter(|l| !l.is_empty()).collect();
    let mut bt_lines: Vec<&str> = bt_str.lines().filter(|l| !l.is_empty()).collect();
    ours_lines.sort_unstable();
    bt_lines.sort_unstable();

    assert_eq!(ours_lines, bt_lines, "output differs from bedtools closest");
}
