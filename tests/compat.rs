use std::io::Write;
use std::path::Path;
use std::process::Command;

use rsomics_bed_closest::closest;
use tempfile::NamedTempFile;

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn bedtools_present() -> bool {
    Command::new("bedtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
}

fn ours(args: &[&std::path::Path], extra: &[&str]) -> String {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rsomics-bed-closest"));
    cmd.arg(args[0]).arg("-b").arg(args[1]).args(extra);
    let out = cmd.output().expect("run rsomics-bed-closest");
    assert!(out.status.success());
    String::from_utf8(out.stdout).unwrap()
}

fn bedtools(a: &Path, b: &Path, extra: &[&str]) -> String {
    let out = Command::new("bedtools")
        .arg("closest")
        .args(extra)
        .arg("-a")
        .arg(a)
        .arg("-b")
        .arg(b)
        .output()
        .expect("run bedtools");
    assert!(out.status.success());
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn matches_bedtools_default() {
    if !bedtools_present() {
        eprintln!("SKIP: bedtools not on PATH");
        return;
    }
    let (a, b) = (golden("a.bed"), golden("b.bed"));
    assert_eq!(
        bedtools(&a, &b, &[]),
        ours(&[&a, &b], &[]),
        "default (no -d)"
    );
}

#[test]
fn matches_bedtools_with_distance() {
    if !bedtools_present() {
        eprintln!("SKIP: bedtools not on PATH");
        return;
    }
    let (a, b) = (golden("a.bed"), golden("b.bed"));
    assert_eq!(
        bedtools(&a, &b, &["-d"]),
        ours(&[&a, &b], &["-d"]),
        "-d mode"
    );
}

#[test]
fn overlap_is_zero_book_ended_is_one() {
    let mut fb = NamedTempFile::new().unwrap();
    writeln!(fb, "chr1\t150\t250").unwrap(); // overlaps [100,200)
    writeln!(fb, "chr1\t200\t300").unwrap(); // book-ended with [100,200)
    let mut fa = NamedTempFile::new().unwrap();
    writeln!(fa, "chr1\t100\t200").unwrap();

    let mut out = Vec::new();
    closest(fa.path(), fb.path(), true, &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    // Only the overlap (distance 0) is closest; the book-ended one is distance 1.
    let lines: Vec<&str> = result.lines().collect();
    assert_eq!(lines.len(), 1, "one closest row: {result}");
    let dist: i64 = lines[0].rsplit('\t').next().unwrap().parse().unwrap();
    assert_eq!(dist, 0, "overlap distance must be 0");
    assert!(
        lines[0].contains("\t150\t250"),
        "overlap B chosen: {}",
        lines[0]
    );
}

#[test]
fn no_b_on_chromosome_emits_placeholder() {
    let fb = NamedTempFile::new().unwrap();
    std::fs::write(fb.path(), "chr2\t10\t20\n").unwrap();
    let fa = NamedTempFile::new().unwrap();
    std::fs::write(fa.path(), "chr1\t100\t200\n").unwrap();

    let mut out = Vec::new();
    closest(fa.path(), fb.path(), false, &mut out).unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "chr1\t100\t200\t.\t-1\t-1\n"
    );
}

// Equidistant B records with the same start must come out in B-file order, like
// bedtools — a stable sort preserves it; sort_unstable could permute them.
#[test]
fn equal_start_ties_keep_file_order() {
    let mut fb = NamedTempFile::new().unwrap();
    writeln!(fb, "chr1\t30\t40\tX").unwrap();
    writeln!(fb, "chr1\t30\t40\tY").unwrap();
    let mut fa = NamedTempFile::new().unwrap();
    writeln!(fa, "chr1\t10\t20").unwrap();

    let mut out = Vec::new();
    closest(fa.path(), fb.path(), false, &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    let names: Vec<&str> = result
        .lines()
        .map(|l| l.rsplit('\t').next().unwrap())
        .collect();
    assert_eq!(names, ["X", "Y"], "ties must keep B-file order");
}
