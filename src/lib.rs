//! Find the closest feature in B for each A interval — bedtools closest equivalent.
//!
//! For each A interval, reports the nearest B interval(s). If A overlaps B,
//! distance is 0 and the overlapping B is reported. Both files must be sorted
//! by (chrom, start) — the same requirement as `bedtools closest`.
//!
//! Output format: A columns, then B columns, then distance (bp).
//! When multiple B intervals are equidistant, all are reported (one per line),
//! matching `bedtools closest` default.
//!
//! When there is no B interval on the same chromosome, the A interval is
//! emitted with `.` for all B columns and `-1` for distance, matching
//! bedtools `-io` absent-handling (default mode).
//!
//! Algorithm: B intervals are grouped by chromosome into sorted Vec; for each
//! A record, a binary search finds the candidate and then neighbours are
//! checked for equidistant ties. O(N log M) per chromosome.

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
    // Sort each chrom's B records by start.
    for v in map.values_mut() {
        v.sort_unstable_by_key(|r| r.start);
    }
    Ok(map)
}

/// Distance between interval [as, ae) and [bs, be) — 0 if they overlap.
fn distance(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> i64 {
    // Overlap condition: neither interval is entirely to the left/right of the other.
    // clippy::manual_range_contains does not apply here — this is a 2-interval
    // overlap test, not a single-point containment test.
    #[allow(clippy::manual_range_contains)]
    if a_start < b_end && b_start < a_end {
        // Overlapping.
        0
    } else if a_end <= b_start {
        // A is left of B.
        (b_start - a_end) as i64
    } else {
        // B is left of A.
        (a_start - b_end) as i64
    }
}

/// Find all B records at minimum distance to [a_start, a_end).
fn find_closest(b: &[BRecord], a_start: u64, a_end: u64) -> Vec<(&BRecord, i64)> {
    if b.is_empty() {
        return vec![];
    }

    // Binary search for first B record with start >= a_start.
    let pos = b.partition_point(|r| r.start < a_start);

    let mut best_dist = i64::MAX;
    let mut best: Vec<(&BRecord, i64)> = Vec::new();

    // Expand outward from pos to find the closest record(s).
    // Scan left from pos-1 and right from pos.
    let check = |idx: usize| -> i64 { distance(a_start, a_end, b[idx].start, b[idx].end) };

    // Scan right from pos.
    let mut r = pos;
    while r < b.len() {
        let d = check(r);
        if d < best_dist {
            best_dist = d;
            best.clear();
        }
        if d == best_dist {
            best.push((&b[r], d));
        }
        // Once the B record starts beyond a_end + best_dist, no closer record can follow.
        // Guard against overflow when best_dist is i64::MAX (no hit found yet).
        if (0..i64::MAX).contains(&best_dist) && b[r].start > a_end + best_dist as u64 {
            break;
        }
        r += 1;
    }

    // Scan left from pos-1.
    // Records are sorted by start, so as l decreases, b[l].start also decreases.
    // The minimum possible distance between A and any record at index ≤ l is
    //   a_start.saturating_sub(b[l].end)
    // because b[l].end ≥ b[l].start and b[j].start ≤ b[l].start for j ≤ l.
    // Once that lower-bound exceeds best_dist, no improvement is possible to the left.
    if pos > 0 {
        let mut l = pos - 1;
        loop {
            // Early exit: all records to the left have end ≤ b[l].end ≤ a_start (or less),
            // so their distance ≥ a_start - b[l].end ≥ best_dist already found.
            if best_dist >= 0 && a_start.saturating_sub(b[l].end) as i64 > best_dist {
                break;
            }
            let d = check(l);
            if d < best_dist {
                best_dist = d;
                // Drop right-side records that are no longer optimal.
                best.retain(|(_, bd)| *bd == best_dist);
            }
            if d == best_dist {
                best.push((&b[l], d));
            }
            if l == 0 {
                break;
            }
            l -= 1;
        }
    }

    // Filter to only best_dist.
    best.retain(|(_, d)| *d == best_dist);
    best
}

/// Run closest on file A vs file B, writing to `output`.
///
/// Output columns: all A columns, all B columns, distance.
pub fn closest(a_path: &Path, b_path: &Path, output: &mut dyn Write) -> Result<()> {
    let b_map = load_b(b_path)?;
    let file = File::open(a_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", a_path.display())))?;
    closest_reader(BufReader::new(file), &b_map, output)
}

/// Same as [`closest`] but reads A from stdin.
pub fn closest_stdin(b_path: &Path, output: &mut dyn Write) -> Result<()> {
    let b_map = load_b(b_path)?;
    closest_reader(BufReader::new(io::stdin()), &b_map, output)
}

fn closest_reader<R: io::Read>(
    reader: BufReader<R>,
    b_map: &HashMap<String, Vec<BRecord>>,
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

        if hits.is_empty() {
            // No B on this chrom — emit with null B and distance -1.
            write!(out, "{chrom}\t{start}\t{end}").map_err(RsomicsError::Io)?;
            if !a_rest.is_empty() {
                write!(out, "\t{a_rest}").map_err(RsomicsError::Io)?;
            }
            out.write_all(b"\t.\t.\t.\t-1\n")
                .map_err(RsomicsError::Io)?;
        } else {
            for (b_rec, dist) in hits {
                write!(out, "{chrom}\t{start}\t{end}").map_err(RsomicsError::Io)?;
                if !a_rest.is_empty() {
                    write!(out, "\t{a_rest}").map_err(RsomicsError::Io)?;
                }
                write!(out, "\t{chrom}\t{}\t{}", b_rec.start, b_rec.end)
                    .map_err(RsomicsError::Io)?;
                if !b_rec.rest.is_empty() {
                    write!(out, "\t{}", b_rec.rest).map_err(RsomicsError::Io)?;
                }
                writeln!(out, "\t{dist}").map_err(RsomicsError::Io)?;
            }
        }
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}
