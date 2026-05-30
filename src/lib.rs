//! Find the closest feature in B for each A interval — `bedtools closest`.
//!
//! Each A row is paired with the nearest B interval(s) on the same chromosome;
//! all B at the minimum distance are emitted, in B-file order (the bedtools tie
//! rule). Output is A columns then B columns; `report_distance` adds a trailing
//! signed-distance column (`bedtools closest -d`). Overlap is distance 0,
//! book-ended is 1. When no B shares the chromosome, B columns are `.`/-1/-1.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Clone)]
struct BRecord {
    start: u64,
    end: u64,
    rest: String,
}

/// Load B BED file; raw lines are stored for output, intervals grouped by chrom.
fn load_b(path: &Path) -> Result<HashMap<String, Vec<BRecord>>> {
    let raw = std::fs::read_to_string(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut map: HashMap<String, Vec<BRecord>> = HashMap::new();
    for line in raw.lines() {
        let bytes = line.as_bytes();
        if bytes.is_empty()
            || bytes[0] == b'#'
            || bytes.starts_with(b"track")
            || bytes.starts_with(b"browser")
        {
            continue;
        }
        let mut fields = line.splitn(4, '\t');
        let chrom = fields.next().unwrap_or("").to_owned();
        let start_str = fields.next().unwrap_or("0");
        let end_str = fields.next().unwrap_or("0");
        let rest = fields.next().unwrap_or("").to_owned();
        let start: u64 = start_str.parse().unwrap_or(0);
        let end: u64 = end_str.parse().unwrap_or(0);
        map.entry(chrom)
            .or_default()
            .push(BRecord { start, end, rest });
    }
    // Stable sort by start so equal-start records keep file order — bedtools
    // breaks distance ties by B-file order.
    for v in map.values_mut() {
        v.sort_by_key(|r| r.start);
    }
    Ok(map)
}

/// Distance between [a_start, a_end) and [b_start, b_end): 0 if they overlap,
/// otherwise the gap plus one, matching `bedtools closest -d` (book-ended
/// features are distance 1, not 0, so they never tie with a true overlap).
#[allow(clippy::manual_range_contains)]
fn distance(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> i64 {
    if a_start < b_end && b_start < a_end {
        0
    } else if a_end <= b_start {
        (b_start - a_end) as i64 + 1
    } else {
        (a_start - b_end) as i64 + 1
    }
}

/// All B records at minimum distance to [a_start, a_end), in B-file (start) order
/// — the tie semantics of `bedtools closest`. B is sorted by start; two passes
/// (min, then collect) keep ties complete and ordered, which an early-exit scan
/// gets wrong when several B share the minimum.
fn find_closest(b: &[BRecord], a_start: u64, a_end: u64) -> Vec<(&BRecord, i64)> {
    let mut best = i64::MAX;
    for r in b {
        let d = distance(a_start, a_end, r.start, r.end);
        if d < best {
            best = d;
        }
    }
    if best == i64::MAX {
        return vec![];
    }
    b.iter()
        .map(|r| (r, distance(a_start, a_end, r.start, r.end)))
        .filter(|(_, d)| *d == best)
        .collect()
}

/// Run closest on file A vs file B, writing to `output`. `report_distance` appends
/// a signed-distance column (`bedtools closest -d`); off by default, matching bedtools.
pub fn closest(
    a_path: &Path,
    b_path: &Path,
    report_distance: bool,
    output: &mut dyn Write,
) -> Result<()> {
    let b_map = load_b(b_path)?;
    let file = File::open(a_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", a_path.display())))?;
    closest_reader(BufReader::new(file), &b_map, report_distance, output)
}

/// Same as [`closest`] but reads A from stdin.
pub fn closest_stdin(b_path: &Path, report_distance: bool, output: &mut dyn Write) -> Result<()> {
    let b_map = load_b(b_path)?;
    closest_reader(BufReader::new(io::stdin()), &b_map, report_distance, output)
}

fn closest_reader<R: io::Read>(
    reader: BufReader<R>,
    b_map: &HashMap<String, Vec<BRecord>>,
    report_distance: bool,
    output: &mut dyn Write,
) -> Result<()> {
    let mut out = BufWriter::new(output);

    for (lineno_0, line) in reader.lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        let bytes = line.as_bytes();

        if bytes.is_empty()
            || bytes[0] == b'#'
            || bytes.starts_with(b"track")
            || bytes.starts_with(b"browser")
        {
            continue;
        }
        let lineno = lineno_0 + 1;
        let mut fields = line.splitn(4, '\t');
        let chrom = fields.next().unwrap_or("");
        let start_str = fields.next().unwrap_or("");
        let end_str = fields.next().unwrap_or("");
        let a_rest = fields.next().unwrap_or("");

        let start: u64 = start_str.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!("line {lineno}: bad start {start_str:?}"))
        })?;
        let end: u64 = end_str.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!("line {lineno}: bad end {end_str:?}"))
        })?;

        let b_ivs = b_map.get(chrom).map(Vec::as_slice).unwrap_or(&[]);
        let hits = find_closest(b_ivs, start, end);

        write!(out, "{chrom}\t{start}\t{end}").map_err(RsomicsError::Io)?;
        if !a_rest.is_empty() {
            write!(out, "\t{a_rest}").map_err(RsomicsError::Io)?;
        }
        if hits.is_empty() {
            // No B on this chromosome: bedtools emits `.`/-1/-1 (then -1 distance under -d).
            out.write_all(b"\t.\t-1\t-1").map_err(RsomicsError::Io)?;
            if report_distance {
                out.write_all(b"\t-1").map_err(RsomicsError::Io)?;
            }
            out.write_all(b"\n").map_err(RsomicsError::Io)?;
        } else {
            let a_prefix = if a_rest.is_empty() {
                format!("{chrom}\t{start}\t{end}")
            } else {
                format!("{chrom}\t{start}\t{end}\t{a_rest}")
            };
            for (idx, (b_rec, dist)) in hits.iter().enumerate() {
                if idx > 0 {
                    write!(out, "{a_prefix}").map_err(RsomicsError::Io)?;
                }
                write!(out, "\t{chrom}\t{}\t{}", b_rec.start, b_rec.end)
                    .map_err(RsomicsError::Io)?;
                if !b_rec.rest.is_empty() {
                    write!(out, "\t{}", b_rec.rest).map_err(RsomicsError::Io)?;
                }
                if report_distance {
                    write!(out, "\t{dist}").map_err(RsomicsError::Io)?;
                }
                out.write_all(b"\n").map_err(RsomicsError::Io)?;
            }
        }
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}
