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

use flate2::read::GzDecoder;
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

/// Per-file N-content scan result (nspector).
///
/// `format_error` is the 0-based index of the first N-containing record whose
/// name fails coordinate parsing; when set, the reference implementation stops
/// reading the file and skips it entirely (parity: data after the error is
/// neither produced nor needed).
///
/// `tiles`/`xs`/`ys`/`n_positions` are parallel: entry `i` describes the
/// i-th N-containing record (records without Ns are omitted).
struct FileScan {
    format_error: Option<usize>,
    tiles: Vec<i32>,
    xs: Vec<i32>,
    ys: Vec<i32>,
    n_positions: Vec<Vec<u16>>,
}

/// Port of `biomate.nspector.nspector.extract_tile_coordinates`.
///
/// Python `int()` trims ASCII whitespace and allows a leading sign; `parse`
/// after `trim()` matches that for the realistic (digit-only) coordinate
/// fields. Returns None for <7 colon-fields or non-numeric fields, exactly
/// like the reference.
fn parse_tile_coordinates(name: &str) -> Option<(i32, i32, i32)> {
    let base = name.split(' ').next()?;
    let parts = base.split(':').collect::<Vec<_>>();
    if parts.len() < 7 {
        return None;
    }
    let tile = parts[4].trim().parse::<i32>().ok()?;
    let x = parts[5].trim().parse::<i32>().ok()?;
    let y = parts[6].trim().parse::<i32>().ok()?;
    Some((tile, x, y))
}

/// Scan one FASTQ (`.gz` suffix => gzip) for N-containing records.
///
/// Line protocol: `@name` / sequence / `+` / qualities. A name line not
/// starting with `@`, or a truncated final record, is a hard error (the
/// reference raises a `dnaio` format error in the same situations).
fn scan_fastq_file(path: &str) -> Result<FileScan, String> {
    let file = std::fs::File::open(path).map_err(|e| e.to_string())?;
    let mut reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(std::io::BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(std::io::BufReader::new(file))
    };

    let mut out = FileScan {
        format_error: None,
        tiles: Vec::new(),
        xs: Vec::new(),
        ys: Vec::new(),
        n_positions: Vec::new(),
    };
    let mut name_line = String::new();
    let mut line = String::new();
    let mut record_idx = 0usize;

    loop {
        name_line.clear();
        let n = reader
            .read_line(&mut name_line)
            .map_err(|e| e.to_string())?;
        if n == 0 {
            break; // clean EOF
        }
        if !name_line.starts_with('@') {
            return Err(format!(
                "malformed FASTQ at record {record_idx}: name line missing '@'"
            ));
        }
        // sequence: scan for N positions immediately (line is reused below)
        line.clear();
        reader.read_line(&mut line).map_err(|e| e.to_string())?;
        let positions: Vec<u16> = line
            .trim_end_matches(['\n', '\r'])
            .bytes()
            .enumerate()
            .filter(|(_, b)| *b == b'N')
            .map(|(i, _)| i as u16)
            .collect();
        // '+' separator
        line.clear();
        reader.read_line(&mut line).map_err(|e| e.to_string())?;
        // qualities
        line.clear();
        let nq = reader.read_line(&mut line).map_err(|e| e.to_string())?;
        if nq == 0 {
            return Err(format!("truncated FASTQ at record {record_idx}"));
        }

        // dnaio strips the leading '@' from record.name
        let rec_name = &name_line.trim_end_matches(['\n', '\r'])[1..];

        if !positions.is_empty() {
            match parse_tile_coordinates(rec_name) {
                Some((tile, x, y)) => {
                    out.tiles.push(tile);
                    out.xs.push(x);
                    out.ys.push(y);
                    out.n_positions.push(positions);
                }
                None => {
                    out.format_error = Some(record_idx);
                    break; // reference stops reading this file
                }
            }
        }
        record_idx += 1;
    }
    Ok(out)
}

/// Scan N content of one or more FASTQ files.
///
/// `paths` are processed by up to `threads` OS threads (work-stealing over
/// the list); results come back in input order. Returns, per file:
/// `(format_error, tiles, xs, ys, n_positions)` where `format_error` is the
/// index of the first bad record or None, and the four arrays are parallel
/// over the N-containing records (see `FileScan`).
#[pyfunction]
fn scan_n_content(
    paths: Vec<String>,
    threads: u32,
) -> PyResult<Vec<(Option<usize>, Vec<i32>, Vec<i32>, Vec<i32>, Vec<Vec<u16>>)>> {
    let n = paths.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let nthreads = (threads.max(1) as usize).min(n);

    let scanned: Vec<Result<FileScan, String>> = if nthreads <= 1 {
        paths.iter().map(|p| scan_fastq_file(p)).collect()
    } else {
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
                tx.send((i, scan_fastq_file(&paths[i]))).ok();
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
    };

    scanned
        .into_iter()
        .map(|r| match r {
            Ok(s) => Ok((s.format_error, s.tiles, s.xs, s.ys, s.n_positions)),
            Err(e) => Err(PyValueError::new_err(e)),
        })
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
}

#[pymodule]
fn biomate_kernels(m: &pyo3::Bound<'_, pyo3::types::PyModule>) -> PyResult<()> {
    m.add_class::<IndexPattern>()?;
    m.add_function(wrap_pyfunction!(scan_n_content, m)?)?;
    Ok(())
}
