use criterion::{Criterion, criterion_group, criterion_main};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;

const N_RECORDS: usize = 50_000;
const CHROM_SIZE: u64 = 100_000_000;
const SEED: u64 = 0x00BE_DC10_5E51;

fn xorshift(x: &mut u64) -> u64 {
    *x ^= *x << 13;
    *x ^= *x >> 7;
    *x ^= *x << 17;
    *x
}

fn synth_bed(path: &PathBuf, n: usize, seed: u64) {
    let f = File::create(path).expect("create bench fixture");
    let mut w = BufWriter::new(f);
    let mut rng = seed;
    let chroms = ["chr1", "chr2", "chr3", "chr4", "chr5"];
    for i in 0..n {
        let chrom = chroms[(xorshift(&mut rng) % chroms.len() as u64) as usize];
        let start = xorshift(&mut rng) % (CHROM_SIZE - 1000);
        let end = start + 100 + (xorshift(&mut rng) % 900);
        writeln!(w, "{chrom}\t{start}\t{end}\tfeature{i}").unwrap();
    }
}

fn sort_bed(path: &PathBuf, sorted: &PathBuf) {
    let out = Command::new("sort")
        .args(["-k1,1", "-k2,2n"])
        .arg(path)
        .output()
        .expect("sort");
    std::fs::write(sorted, &out.stdout).expect("write sorted");
}

fn ensure_fixtures() -> (PathBuf, PathBuf) {
    let mut a = std::env::temp_dir();
    a.push(format!("rsomics-bed-closest-bench-a-{N_RECORDS}.bed"));
    let mut a_sorted = std::env::temp_dir();
    a_sorted.push(format!(
        "rsomics-bed-closest-bench-a-{N_RECORDS}-sorted.bed"
    ));
    let mut b_sorted = std::env::temp_dir();
    b_sorted.push(format!(
        "rsomics-bed-closest-bench-b-{N_RECORDS}-sorted.bed"
    ));
    let mut b = std::env::temp_dir();
    b.push(format!("rsomics-bed-closest-bench-b-{N_RECORDS}.bed"));
    if !a_sorted.exists() {
        synth_bed(&a, N_RECORDS, SEED);
        sort_bed(&a, &a_sorted);
    }
    if !b_sorted.exists() {
        synth_bed(&b, N_RECORDS, SEED ^ 0xDEAD_BEEF);
        sort_bed(&b, &b_sorted);
    }
    (a_sorted, b_sorted)
}

fn bench(c: &mut Criterion) {
    let (a, b) = ensure_fixtures();
    let ours = env!("CARGO_BIN_EXE_rsomics-bed-closest");
    let mut group = c.benchmark_group(format!("bed_closest/{N_RECORDS}"));
    group.sample_size(10);

    group.bench_function("rsomics-bed-closest", |bm| {
        bm.iter(|| {
            let out = Command::new(ours)
                .arg(&a)
                .arg("-b")
                .arg(&b)
                .output()
                .expect("ours run");
            assert!(
                out.status.success(),
                "{}",
                String::from_utf8_lossy(&out.stderr)
            );
        });
    });

    if Command::new("bedtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
    {
        group.bench_function("bedtools-closest", |bm| {
            bm.iter(|| {
                let out = Command::new("bedtools")
                    .args(["closest", "-a"])
                    .arg(&a)
                    .arg("-b")
                    .arg(&b)
                    .output()
                    .expect("bedtools run");
                assert!(
                    out.status.success(),
                    "{}",
                    String::from_utf8_lossy(&out.stderr)
                );
            });
        });
    } else {
        eprintln!("bedtools not on PATH — skipping upstream comparison");
    }

    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
