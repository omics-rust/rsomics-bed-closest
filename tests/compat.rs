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
    // A1→1 line, A2→2 lines (equidistant tie: B1 and B2 both 100 bp away), A3→1 line = 4.
    assert_eq!(
        lines.len(),
        4,
        "expected 4 output lines (A2 ties): {result}"
    );

    // A1 [100,200) closest to B1 [300,400) — distance 100.
    assert!(
        lines[0].contains("chr1\t300\t400"),
        "A1 closest wrong: {}",
        lines[0]
    );
    let dist0: i64 = lines[0]
        .split('\t')
        .next_back()
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert_eq!(dist0, 100, "A1 distance wrong: {dist0}");

    // A2 [500,600): equidistant to B1 [300,400) and B2 [700,800) — both at distance 100.
    for idx in 1..=2 {
        let dist: i64 = lines[idx]
            .split('\t')
            .next_back()
            .unwrap()
            .trim()
            .parse()
            .unwrap();
        assert_eq!(dist, 100, "A2 line[{idx}] distance wrong: {dist}");
        assert!(
            lines[idx].starts_with("chr1\t500\t600"),
            "A2 line[{idx}] A columns wrong: {}",
            lines[idx]
        );
    }

    // A3 chr2:[100,200) closest to B3 chr2:[250,350) — distance 50.
    assert!(
        lines[3].contains("chr2\t250\t350"),
        "A3 closest wrong: {}",
        lines[3]
    );
    let dist3: i64 = lines[3]
        .split('\t')
        .next_back()
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert_eq!(dist3, 50, "A3 distance wrong: {dist3}");
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
        .next_back()
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

    // bedtools closest (without -d) emits A+B columns only; we always emit
    // the distance as a trailing column. Strip our distance column before
    // comparing — the distance value itself is verified in
    // basic_closest_correctness.
    fn strip_last_col(s: &str) -> &str {
        s.rfind('\t').map(|i| &s[..i]).unwrap_or(s)
    }

    let bt = Command::new("bedtools")
        .args(["closest", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .output()
        .expect("bedtools closest failed");
    let bt_str = String::from_utf8(bt.stdout).unwrap();

    let mut ours_lines: Vec<&str> = ours_str
        .lines()
        .filter(|l| !l.is_empty())
        .map(strip_last_col)
        .collect();
    let bt_lines: Vec<&str> = bt_str.lines().filter(|l| !l.is_empty()).collect();

    // Sort both for order-independent comparison.
    let mut bt_sorted = bt_lines.clone();
    ours_lines.sort_unstable();
    bt_sorted.sort_unstable();

    assert_eq!(
        ours_lines, bt_sorted,
        "A+B columns differ from bedtools closest"
    );
}
