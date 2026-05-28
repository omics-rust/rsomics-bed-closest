use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bed_closest(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bed-closest");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let a = manifest.join("tests/golden/a.bed");
    let b = manifest.join("tests/golden/b.bed");
    c.bench_function("rsomics-bed-closest golden", |b_| {
        b_.iter(|| {
            let out = Command::new(black_box(bin))
                .args(["-a", a.to_str().unwrap(), "-b", b.to_str().unwrap()])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_bed_closest);
criterion_main!(benches);
