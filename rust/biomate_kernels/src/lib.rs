//! Native kernels for BioMate.
//!
//! `IndexPattern` is a parity-checked native replacement for the `regex`-module
//! fuzzy patterns built by `biomate.index.index.compile_index_regex`.
//!
//! # Parity model (verified against the `regex` C engine, mrab-regex `_regex.c`)
//!
//! Patterns have the shape `((PART1){err}[.*PART2{err}...]*){err}`.
//!
//! 1. Leftmost successful match start wins, even with more errors than a
//!    later start.
//! 2. At each mismatch a depth-first search tries SUB -> INS -> DEL; the first
//!    completing path wins. Counts are path-dependent, NOT minimum-edit.
//! 3. One global (sub, ins, del) counter is shared by all fuzzy nodes. A
//!    commit needs `count[t] < caps[t]` and `total < total_cap`.
//!    - `e<=d`: caps (d,d,d), total_cap = d.
//!    - `i<=a,d<=b,s<=c`: caps (c,a,b), total_cap = a+b+c (SUM; probed).
//! 4. INS is forbidden only at the very start of the whole match (text_pos ==
//!    match start); always allowed elsewhere.
//! 5. A string part is scanned exactly; on mismatch one error is committed and
//!    the scan resumes. The engine's `fuzzy_ext_match` lookahead is a no-op for
//!    our node shapes (next node is group machinery).
//! 6. Trailing lazy insertions and multi-copy repetitions never change the
//!    first-success outcome (probed); omitted.
//! 7. The `.*` gap is greedy: rightmost part start first, backtracking left;
//!    at a fixed gap length ALL error paths of the following part are explored
//!    before the gap shrinks.
//! 8. `fuzzy_counts` tuple order is (subs, ins, dels).

use flate2::read::MultiGzDecoder;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::{BTreeSet, HashMap};
use std::io::BufRead;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

const MAX_TEXT: usize = 1_000_000; // sanity guard per record

#[derive(Clone, Copy)]
struct Budget {
    /// caps by error type: 0 = sub, 1 = ins, 2 = del
    caps: [u32; 3],
    total_cap: u32,
}

struct Ctx<'t> {
    text: &'t [u8],
    /// Search origin (0 for `pat.search(seq)`); INS is forbidden only here.
    anchor: usize,
    start: usize,
    budget: Budget,
}

/// Full DFS over parts from part `k`, part cursor `pos`, text cursor `tp`.
/// Returns the final text position of the first (engine-ordered) successful
/// completion; `counts` restored on failure. `spans` holds per-part
/// (start, end); `spans[k].0` is set on entry, `.1` on completion.
fn dfs(
    ctx: &Ctx,
    parts: &[&[u8]],
    k: usize,
    pos: usize,
    tp: usize,
    counts: &mut [u32; 3],
    spans: &mut Vec<(usize, usize)>,
) -> Option<usize> {
    let text = ctx.text;
    let text_len = text.len();
    let part = parts[k];
    if pos == part.len() {
        // Part complete.
        spans[k].1 = tp;
        if k + 1 == parts.len() {
            return Some(tp);
        }
        // Greedy `.*` gap before part k+1: rightmost first.
        for gtp in (tp..=text_len).rev() {
            spans.push((gtp, 0));
            if let Some(r) = dfs(ctx, parts, k + 1, 0, gtp, counts, spans) {
                return Some(r);
            }
            spans.pop();
        }
        return None;
    }
    if tp >= text_len {
        // Only deletions of the whole remaining part can still work.
        let rest = (part.len() - pos) as u32;
        let total: u32 = counts.iter().sum::<u32>();
        if rest <= ctx.budget.total_cap - total && counts[2] + rest <= ctx.budget.caps[2] {
            counts[2] += rest;
            return dfs(ctx, parts, k, part.len(), tp, counts, spans);
        }
        return None;
    }
    // Fast-fail: even deleting every remaining part char cannot bridge the
    // length deficit.
    let remaining = ctx.budget.total_cap - counts.iter().sum::<u32>();
    if (part.len() - pos) as u32 > (text_len - tp) as u32 + remaining {
        return None;
    }
    if part[pos] == text[tp] {
        return dfs(ctx, parts, k, pos + 1, tp + 1, counts, spans);
    }
    // Mismatch: engine order SUB -> INS -> DEL.
    let total: u32 = counts.iter().sum::<u32>();
    // SUB
    if total < ctx.budget.total_cap && counts[0] < ctx.budget.caps[0] {
        counts[0] += 1;
        if let Some(r) = dfs(ctx, parts, k, pos + 1, tp + 1, counts, spans) {
            return Some(r);
        }
        counts[0] -= 1;
    }
    // INS (forbidden only at the search origin, not at each match start)
    if tp != ctx.anchor && total < ctx.budget.total_cap && counts[1] < ctx.budget.caps[1] {
        counts[1] += 1;
        if let Some(r) = dfs(ctx, parts, k, pos, tp + 1, counts, spans) {
            return Some(r);
        }
        counts[1] -= 1;
    }
    // DEL
    if total < ctx.budget.total_cap && counts[2] < ctx.budget.caps[2] {
        counts[2] += 1;
        if let Some(r) = dfs(ctx, parts, k, pos + 1, tp, counts, spans) {
            return Some(r);
        }
        counts[2] -= 1;
    }
    None
}

/// Leftmost-start driver: returns (start, end, counts, spans) of the first
/// successful match, or None.
fn find_match(
    text: &[u8],
    parts: &[&[u8]],
    budget: Budget,
    starts: &mut dyn Iterator<Item = usize>,
) -> Option<(usize, usize, [u32; 3], Vec<(usize, usize)>)> {
    while let Some(start) = starts.next() {
        let ctx = Ctx {
            text,
            anchor: 0,
            start,
            budget,
        };
        let mut counts = [0u32; 3];
        let mut spans = vec![(start, 0)];
        if let Some(end) = dfs(&ctx, parts, 0, 0, start, &mut counts, &mut spans) {
            return Some((start, end, counts, spans));
        }
    }
    None
}

/// Seed-and-verify start candidates: split part 0 into total_cap+1 blocks;
/// any match within budget leaves >=1 block intact at offset off_j shifted by
/// <= total_cap insertions. Yields candidate starts in ascending order.
fn seeded_starts(text: &[u8], part0: &[u8], d: usize) -> BTreeSet<usize> {
    let mut cands = BTreeSet::new();
    let nb = d + 1;
    let block_len = part0.len() / nb;
    if block_len == 0 {
        // Part shorter than the budget: brute-force every start.
        for s in 0..=text.len() {
            cands.insert(s);
        }
        return cands;
    }
    for j in 0..nb {
        let off = j * block_len;
        let block = &part0[off..off + block_len];
        let mut p = 0usize;
        while p + block.len() <= text.len() {
            let mut found = None;
            let mut rel = 0usize;
            while rel + block.len() <= text.len() - p {
                if &text[p + rel..p + rel + block.len()] == block {
                    found = Some(rel);
                    break;
                }
                rel += 1;
            }
            match found {
                Some(rel) => {
                    let pos = p + rel;
                    // Block may be shifted by -d..=+d (dels move it left,
                    // ins right).
                    for kshift in -(d as i64)..=(d as i64) {
                        let bpos = off as i64 + kshift;
                        if bpos >= 0 && (bpos as usize) <= pos {
                            cands.insert(pos - bpos as usize);
                        }
                    }
                    p = pos + 1;
                }
                None => break,
            }
        }
    }
    cands
}

fn all_starts(text: &[u8]) -> std::ops::RangeInclusive<usize> {
    0..=text.len()
}

#[pyclass]
#[derive(Clone)]
pub struct IndexPattern {
    parts: Vec<Vec<u8>>,
    budget: Budget,
}

#[pymethods]
impl IndexPattern {
    /// Budget semantics mirroring `compile_index_regex`.
    #[new]
    fn new(motif: &str, distance: u32, error_type: &str) -> PyResult<Self> {
        let motif = motif.to_uppercase();
        let parts: Vec<Vec<u8>> = motif
            .split('+')
            .map(|p| {
                if p.is_empty() {
                    Err(PyValueError::new_err("empty motif part"))
                } else {
                    Ok(p.as_bytes().to_vec())
                }
            })
            .collect::<Result<_, _>>()?;
        let budget = if error_type == "e" {
            Budget {
                caps: [distance, distance, distance],
                total_cap: distance,
            }
        } else {
            let cap = |c: char| if error_type.contains(c) { distance } else { 0 };
            let (c_sub, c_ins, c_del) = (cap('s'), cap('i'), cap('d'));
            Budget {
                caps: [c_sub, c_ins, c_del],
                total_cap: c_sub + c_ins + c_del,
            }
        };
        Ok(Self { parts, budget })
    }

    /// Single-record search: `(start, end, (subs, ins, dels), [part texts])`
    /// or None. Returns None for non-ASCII / oversized input (caller falls
    /// back to the `regex` engine).
    fn search(&self, seq: &str) -> Option<(usize, usize, (u32, u32, u32), Vec<String>)> {
        if seq.len() > MAX_TEXT || seq.bytes().any(|b| b > 127) {
            return None;
        }
        let text: Vec<u8> = seq.bytes().collect();
        let text = text.as_slice();
        let parts: Vec<&[u8]> = self.parts.iter().map(|p| p.as_slice()).collect();
        let d = self.budget.total_cap as usize;
        let part0 = parts[0];
        let block_len = part0.len() / (d + 1);
        let cands: Option<BTreeSet<usize>> = if block_len >= 1 {
            Some(seeded_starts(text, part0, d))
        } else {
            None
        };
        let res = match cands {
            Some(c) => {
                let starts = c.iter().copied();
                find_match(text, &parts, self.budget, &mut Box::new(starts))
            }
            None => {
                let s = all_starts(text);
                find_match(text, &parts, self.budget, &mut Box::new(s))
            }
        }?;
        let (start, end, counts, spans) = res;
        let part_texts: Vec<String> = spans
            .iter()
            .map(|(a, b)| String::from_utf8_lossy(&text[*a..*b]).into_owned())
            .collect();
        let _ = start;
        let _ = end;
        Some((start, end, (counts[0], counts[1], counts[2]), part_texts))
    }

    /// Batch scan aggregating `results[fuzziness][groups] = count` where
    /// fuzziness = str(subs+ins+dels) and groups = part texts space-joined
    /// (mirrors `" ".join(match.groups()[1:])` used by index.py).
    /// Returns `(results, fallback_indices)`.
    fn scan_batch(
        &self,
        seqs: Vec<String>,
    ) -> (HashMap<String, HashMap<String, u64>>, Vec<usize>) {
        let mut results: HashMap<String, HashMap<String, u64>> = HashMap::new();
        let mut fallback: Vec<usize> = Vec::new();
        let parts: Vec<&[u8]> = self.parts.iter().map(|p| p.as_slice()).collect();
        let d = self.budget.total_cap as usize;
        let part0 = parts[0];
        let block_len = part0.len() / (d + 1);
        for (i, seq) in seqs.iter().enumerate() {
            if seq.len() > MAX_TEXT || seq.bytes().any(|b| b > 127) {
                fallback.push(i);
                continue;
            }
            let text: Vec<u8> = seq.bytes().collect();
            let text = text.as_slice();
            let cands: Option<BTreeSet<usize>> = if block_len >= 1 {
                Some(seeded_starts(text, part0, d))
            } else {
                None
            };
            let res = match cands {
                Some(c) => {
                    let starts = c.iter().copied();
                    find_match(text, &parts, self.budget, &mut Box::new(starts))
                }
                None => {
                    let s = all_starts(text);
                    find_match(text, &parts, self.budget, &mut Box::new(s))
                }
            };
            if let Some((_, _, counts, spans)) = res {
                let fuzz = (counts[0] + counts[1] + counts[2]).to_string();
                let groups: Vec<String> = spans
                    .iter()
                    .map(|(a, b)| String::from_utf8_lossy(&text[*a..*b]).into_owned())
                    .collect();
                let key = groups.join(" ");
                *results.entry(fuzz).or_default().entry(key).or_insert(0) += 1;
            }
        }
        (results, fallback)
    }
}

/// Per-tile accumulator for one FASTQ file (nspector).
#[derive(Default)]
struct TileAcc {
    total: u64,
    reads_with_n: u64,
    n_bases: u64,
    /// reads with at least one called (non-N) base below `min_qual`
    lowq_reads: u64,
    /// reads with an N at cycle `i + 1`
    cycles: Vec<u64>,
    /// reads with a called base below `min_qual` at cycle `i + 1`
    lowq_cycles: Vec<u64>,
}

/// Cycles at each end of the read used to classify where the first N falls.
const EDGE_CYCLES: usize = 10;

/// Region of the first N: 0 = first `EDGE_CYCLES` cycles, 2 = last
/// `EDGE_CYCLES` cycles, 1 = middle (start wins on short reads).
fn first_n_region(first: usize, len: usize) -> u8 {
    if first < EDGE_CYCLES {
        0
    } else if first + EDGE_CYCLES >= len {
        2
    } else {
        1
    }
}

/// Per-file N-content scan result (nspector), aggregated while streaming.
///
/// Memory depends on lanes × tiles × cycles (plus occupied spatial bins),
/// never on the number of reads. `format_error` is the 0-based index of the
/// first record whose name fails coordinate parsing; the file is then
/// abandoned, exactly like the Python reference.
struct FileScan {
    format_error: Option<usize>,
    max_len: u32,
    tiles: HashMap<(i32, i32), TileAcc>,
    n_per_read: HashMap<(i32, u32), u64>,
    spatial: HashMap<(i32, i32, i32, i32), u64>,
    /// (lane, longest N run, first-N region, all-N flag) -> reads
    n_shape: HashMap<(i32, u32, u8, u8), u64>,
}

/// Python tuple returned per file by `scan_n_content`:
/// `(format_error, max_len, tiles, tile_cycle, n_per_read, spatial, n_shape)` with
/// `tiles = [(lane, tile, total, reads_with_n, n_bases, lowq_reads)]`,
/// `tile_cycle = [(lane, tile, cycle, reads_with_n, reads_with_lowq)]`
/// (1-based, rows where either count is non-zero),
/// `n_per_read = [(lane, n_count, reads)]`,
/// `spatial = [(lane, tile, x // bin, y // bin, reads_with_n)]`,
/// `n_shape = [(lane, longest_run, first_region, all_n, reads)]`, all sorted.
type ScanTuple = (
    Option<usize>,
    u32,
    Vec<(i32, i32, u64, u64, u64, u64)>,
    Vec<(i32, i32, u32, u64, u64)>,
    Vec<(i32, u32, u64)>,
    Vec<(i32, i32, i32, i32, u64)>,
    Vec<(i32, u32, u8, u8, u64)>,
);

/// Port of `biomate.nspector.scan.extract_tile_coordinates`.
///
/// Python `int()` trims ASCII whitespace and allows a leading sign; `parse`
/// after `trim()` matches that for the realistic (digit-only) coordinate
/// fields. Returns `(lane, tile, x, y)`, or None for <7 colon-fields or
/// non-numeric fields, exactly like the reference.
fn parse_tile_coordinates(name: &str) -> Option<(i32, i32, i32, i32)> {
    let base = name.split(' ').next()?;
    let parts = base.split(':').collect::<Vec<_>>();
    if parts.len() < 7 {
        return None;
    }
    let lane = parts[3].trim().parse::<i32>().ok()?;
    let tile = parts[4].trim().parse::<i32>().ok()?;
    let x = parts[5].trim().parse::<i32>().ok()?;
    let y = parts[6].trim().parse::<i32>().ok()?;
    Some((lane, tile, x, y))
}

/// `Read` adapter fed by a background thread that decompresses ahead, so
/// gzip inflation and FASTQ parsing run on two cores for a single file.
struct ChannelReader {
    rx: std::sync::mpsc::Receiver<std::io::Result<Vec<u8>>>,
    buf: Vec<u8>,
    pos: usize,
}

impl ChannelReader {
    fn spawn<R: std::io::Read + Send + 'static>(mut inner: R) -> Self {
        let (tx, rx) = std::sync::mpsc::sync_channel(4);
        std::thread::spawn(move || loop {
            let mut chunk = vec![0u8; 1 << 20];
            let mut filled = 0;
            while filled < chunk.len() {
                match inner.read(&mut chunk[filled..]) {
                    Ok(0) => break,
                    Ok(n) => filled += n,
                    Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                    Err(e) => {
                        let _ = tx.send(Err(e));
                        return;
                    }
                }
            }
            chunk.truncate(filled);
            let done = filled == 0;
            if tx.send(Ok(chunk)).is_err() || done {
                return; // receiver gone (early stop) or EOF sent
            }
        });
        ChannelReader { rx, buf: Vec::new(), pos: 0 }
    }
}

impl std::io::Read for ChannelReader {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        if self.pos == self.buf.len() {
            match self.rx.recv() {
                Ok(Ok(chunk)) => {
                    self.buf = chunk;
                    self.pos = 0;
                }
                Ok(Err(e)) => return Err(e),
                Err(_) => return Ok(0),
            }
        }
        let n = out.len().min(self.buf.len() - self.pos);
        out[..n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

fn trim_eol(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && (line[end - 1] == b'\n' || line[end - 1] == b'\r') {
        end -= 1;
    }
    &line[..end]
}

/// Scan one FASTQ (`.gz` suffix => gzip) and aggregate N counters.
///
/// Line protocol: `@name` / sequence / `+` / qualities. A name line not
/// starting with `@`, or a truncated final record, is a hard error (the
/// reference raises a `dnaio` format error in the same situations).
fn scan_fastq_file(path: &str, bin_size: i32, min_qual: u8) -> Result<FileScan, String> {
    // Phred+33: a quality byte below this is a low-quality call
    let lowq_byte = min_qual.saturating_add(33);
    let file = std::fs::File::open(path).map_err(|e| e.to_string())?;
    let mut reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(std::io::BufReader::with_capacity(
            1 << 20,
            ChannelReader::spawn(MultiGzDecoder::new(std::io::BufReader::new(file))),
        ))
    } else {
        Box::new(std::io::BufReader::with_capacity(1 << 20, file))
    };

    let mut out = FileScan {
        format_error: None,
        max_len: 0,
        tiles: HashMap::new(),
        n_per_read: HashMap::new(),
        spatial: HashMap::new(),
        n_shape: HashMap::new(),
    };
    let mut name_line: Vec<u8> = Vec::new();
    let mut seq_line: Vec<u8> = Vec::new();
    let mut line: Vec<u8> = Vec::new();
    let mut record_idx = 0usize;

    loop {
        name_line.clear();
        let n = reader
            .read_until(b'\n', &mut name_line)
            .map_err(|e| e.to_string())?;
        if n == 0 {
            break; // clean EOF
        }
        if !name_line.starts_with(b"@") {
            return Err(format!(
                "malformed FASTQ at record {record_idx}: name line missing '@'"
            ));
        }
        seq_line.clear();
        reader
            .read_until(b'\n', &mut seq_line)
            .map_err(|e| e.to_string())?;
        // '+' separator
        line.clear();
        reader.read_until(b'\n', &mut line).map_err(|e| e.to_string())?;
        // qualities
        line.clear();
        let nq = reader.read_until(b'\n', &mut line).map_err(|e| e.to_string())?;
        if nq == 0 {
            return Err(format!("truncated FASTQ at record {record_idx}"));
        }

        // dnaio strips the leading '@' from record.name
        let rec_name = std::str::from_utf8(&trim_eol(&name_line)[1..]).unwrap_or("");
        let Some((lane, tile, x, y)) = parse_tile_coordinates(rec_name) else {
            out.format_error = Some(record_idx);
            break; // reference stops reading this file
        };

        let seq = trim_eol(&seq_line);
        out.max_len = out.max_len.max(seq.len() as u32);
        let acc = out.tiles.entry((lane, tile)).or_default();
        acc.total += 1;

        let qual = trim_eol(&line);
        let mut n_count = 0u32;
        let mut any_lowq = false;
        let (mut run, mut longest, mut first) = (0u32, 0u32, usize::MAX);
        for (i, b) in seq.iter().enumerate() {
            if *b == b'N' {
                if acc.cycles.len() <= i {
                    acc.cycles.resize(i + 1, 0);
                }
                acc.cycles[i] += 1;
                n_count += 1;
                run += 1;
                longest = longest.max(run);
                first = first.min(i);
            } else {
                run = 0;
                if qual.get(i).is_some_and(|q| *q < lowq_byte) {
                    if acc.lowq_cycles.len() <= i {
                        acc.lowq_cycles.resize(i + 1, 0);
                    }
                    acc.lowq_cycles[i] += 1;
                    any_lowq = true;
                }
            }
        }
        if any_lowq {
            acc.lowq_reads += 1;
        }
        if n_count > 0 {
            acc.reads_with_n += 1;
            acc.n_bases += n_count as u64;
            *out.n_per_read.entry((lane, n_count)).or_insert(0) += 1;
            let all_n = (n_count as usize == seq.len()) as u8;
            *out
                .n_shape
                .entry((lane, longest, first_n_region(first, seq.len()), all_n))
                .or_insert(0) += 1;
            *out
                .spatial
                .entry((lane, tile, x / bin_size, y / bin_size))
                .or_insert(0) += 1;
        }
        record_idx += 1;
    }
    Ok(out)
}

fn to_tuple(s: FileScan) -> ScanTuple {
    let mut tiles = Vec::with_capacity(s.tiles.len());
    let mut tile_cycle = Vec::new();
    for ((lane, tile), acc) in s.tiles {
        tiles.push((
            lane,
            tile,
            acc.total,
            acc.reads_with_n,
            acc.n_bases,
            acc.lowq_reads,
        ));
        for i in 0..acc.cycles.len().max(acc.lowq_cycles.len()) {
            let n = acc.cycles.get(i).copied().unwrap_or(0);
            let lowq = acc.lowq_cycles.get(i).copied().unwrap_or(0);
            if n > 0 || lowq > 0 {
                tile_cycle.push((lane, tile, (i + 1) as u32, n, lowq));
            }
        }
    }
    tiles.sort_unstable();
    tile_cycle.sort_unstable();
    let mut n_per_read: Vec<_> = s.n_per_read.into_iter().map(|((l, n), c)| (l, n, c)).collect();
    n_per_read.sort_unstable();
    let mut spatial: Vec<_> = s
        .spatial
        .into_iter()
        .map(|((l, t, xb, yb), c)| (l, t, xb, yb, c))
        .collect();
    spatial.sort_unstable();
    let mut n_shape: Vec<_> = s
        .n_shape
        .into_iter()
        .map(|((l, run, region, all_n), c)| (l, run, region, all_n, c))
        .collect();
    n_shape.sort_unstable();
    (s.format_error, s.max_len, tiles, tile_cycle, n_per_read, spatial, n_shape)
}

/// Scan N content of one or more FASTQ files.
///
/// `paths` are processed by up to `threads` OS threads (work-stealing over
/// the list); results come back in input order, one `ScanTuple` per file.
/// `bin_size` is the spatial bin width in read-name coordinate units;
/// `min_qual` is the Phred (+33) threshold for low-quality called bases.
#[pyfunction]
#[pyo3(signature = (paths, threads, bin_size=1000, min_qual=20))]
fn scan_n_content(
    py: Python<'_>,
    paths: Vec<String>,
    threads: u32,
    bin_size: i32,
    min_qual: u8,
) -> PyResult<Vec<ScanTuple>> {
    if bin_size < 1 {
        return Err(PyValueError::new_err("bin_size must be >= 1"));
    }
    let n = paths.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let nthreads = (threads.max(1) as usize).min(n);

    let scanned: Vec<Result<FileScan, String>> = py.detach(|| {
        if nthreads <= 1 {
            return paths
                .iter()
                .map(|p| scan_fastq_file(p, bin_size, min_qual))
                .collect();
        }
        let paths = Arc::new(paths);
        let next = Arc::new(AtomicUsize::new(0));
        let (tx, rx) = std::sync::mpsc::channel();
        let mut handles = Vec::new();
        for _ in 0..nthreads {
            let paths = Arc::clone(&paths);
            let next = Arc::clone(&next);
            let tx = tx.clone();
            handles.push(std::thread::spawn(move || loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= n {
                    break;
                }
                tx.send((i, scan_fastq_file(&paths[i], bin_size, min_qual))).ok();
            }));
        }
        drop(tx);
        for h in handles {
            h.join().ok();
        }
        let mut out: Vec<Option<Result<FileScan, String>>> = (0..n).map(|_| None).collect();
        for (i, r) in rx {
            out[i] = Some(r);
        }
        out.into_iter().map(|o| o.expect("result for every index")).collect()
    });

    scanned
        .into_iter()
        .map(|r| r.map(to_tuple).map_err(PyValueError::new_err))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_and_fuzzy() {
        let p = IndexPattern {
            parts: vec!["ATCGATCGAT".as_bytes().to_vec()],
            budget: Budget {
                caps: [2, 2, 2],
                total_cap: 2,
            },
        };
        // leftmost start wins: 2 deletions at start=1 beat the exact start=3
        let m = p.search("NNNATCGATCGATNNN").unwrap();
        assert_eq!(m.0, 1);
        assert_eq!(m.1, 13);
        assert_eq!(m.2, (0, 2, 0));
        let m = p.search("XXXATCGATCGAT").unwrap();
        assert_eq!(m.0, 1); // leftmost wins over the 0-error start at 3
        assert_eq!(m.2, (0, 2, 0));
    }

    #[test]
    fn first_n_region_edges() {
        assert_eq!(first_n_region(0, 150), 0);
        assert_eq!(first_n_region(9, 150), 0);
        assert_eq!(first_n_region(10, 150), 1);
        assert_eq!(first_n_region(139, 150), 1);
        assert_eq!(first_n_region(140, 150), 2);
        assert_eq!(first_n_region(5, 12), 0); // start wins on short reads
        assert_eq!(first_n_region(10, 12), 2);
    }
}

#[pymodule]
fn biomate_kernels(m: &pyo3::Bound<'_, pyo3::types::PyModule>) -> PyResult<()> {
    m.add_class::<IndexPattern>()?;
    m.add_function(wrap_pyfunction!(scan_n_content, m)?)?;
    Ok(())
}
