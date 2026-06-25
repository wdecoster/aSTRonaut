//! Kmer counting, rotational canonicalization and selection.
//!
//! Ports `count_kmers`, `get_rotations`, `prune_counts`, `select_kmers` and
//! `get_kmers` from the Python implementation.

use std::collections::HashMap;

use crate::cli::{Args, KmerArg};
use crate::model::RepeatRecord;

/// Bounds and fallback for `-k auto`, ported from trout's `features.rs`.
/// k=1 is too crude; motifs longer than 6 are rare and dilute the signal.
pub const AUTO_K_MIN: usize = 2;
pub const AUTO_K_MAX: usize = 6;
/// Fallback when no clean period is recoverable (trinucleotide repeats are the
/// most common pathogenic motif length).
pub const AUTO_K_DEFAULT: usize = 3;
/// Composition-corrected self-shift threshold for `detect_period`. The raw
/// self-shift rate is inflated by base-composition bias (an A-rich motif matches
/// itself often just because most bases are A), which made short periods clear a
/// raw threshold spuriously. We instead score the rate *above chance*,
/// `(observed - chance) / (1 - chance)`, so composition cancels out: a pure
/// tandem repeat still scores ~1.0, GC-rich hexamers (C9orf72 GGCCCC) clear it,
/// and A-rich motifs no longer match at short periods. 0.30 is low enough to
/// keep recovering impure/biased repeats while rejecting chance-level periods.
const PERIOD_SCORE_THRESHOLD: f64 = 0.30;

/// Probability that two independently drawn bases of `seq` are identical
/// (`Σ frequencyₐ²`) — the self-match rate expected by chance given the
/// sequence's composition.
fn chance_match_rate(seq: &[u8]) -> f64 {
    let mut counts = [0usize; 256];
    for &b in seq {
        counts[b as usize] += 1;
    }
    let n = seq.len() as f64;
    if n == 0.0 {
        return 0.0;
    }
    counts.iter().map(|&c| (c as f64 / n).powi(2)).sum()
}

/// Detect the dominant tandem-repeat period in `seq`. For each candidate period
/// it measures the self-shift match rate (`seq[i] == seq[i + p]`), corrects it
/// for the sequence's base composition, and returns the smallest period whose
/// corrected score clears [`PERIOD_SCORE_THRESHOLD`] — the fundamental period,
/// since multiples of the true period also score highly. Adapted from trout's
/// `detect_period`, with the composition correction added.
pub fn detect_period(seq: &[u8], k_max: usize) -> Option<usize> {
    if seq.len() < 2 * k_max {
        return None;
    }
    let chance = chance_match_rate(seq);
    let denom = 1.0 - chance;
    if denom <= f64::EPSILON {
        return None; // homopolymer-like: no meaningful period
    }
    for p in AUTO_K_MIN..=k_max {
        let total = seq.len() - p;
        let matches = (0..total).filter(|&i| seq[i] == seq[i + p]).count();
        let observed = matches as f64 / total as f64;
        let score = (observed - chance) / denom;
        if score >= PERIOD_SCORE_THRESHOLD {
            return Some(p);
        }
    }
    None
}

/// Infer k for one repeat locus by consensus.
///
/// Adapted from trout's `detect_locus_k`. trout cross-checked the reference
/// allele against a single median allele, but a single atypical allele could
/// flip the result, making it sample-dependent. Instead we detect the period of
/// the reference allele (VCF REF) *and* every allele, and take the most-voted
/// period — robust to a few noisy alleles. Ties prefer the reference's period,
/// then the smaller period. Falls back to the default when nothing yields a
/// period (e.g. `--table` input with only impure sequences). Returns the k and a
/// provenance label.
pub fn detect_locus_k(records: &[RepeatRecord]) -> (usize, &'static str) {
    // REF is identical for every record of a locus.
    let ref_seq = records
        .iter()
        .map(|r| r.reference.as_str())
        .find(|s| !s.is_empty())
        .unwrap_or("");
    let p_ref = detect_period(ref_seq.as_bytes(), AUTO_K_MAX);

    let mut votes: HashMap<usize, usize> = HashMap::new();
    if let Some(p) = p_ref {
        *votes.entry(p).or_default() += 1;
    }
    for r in records {
        if r.sequence.is_empty() {
            continue;
        }
        if let Some(p) = detect_period(r.sequence.as_bytes(), AUTO_K_MAX) {
            *votes.entry(p).or_default() += 1;
        }
    }

    let Some(&max_votes) = votes.values().max() else {
        return (AUTO_K_DEFAULT, "fallback");
    };
    let mut top: Vec<usize> = votes
        .iter()
        .filter(|(_, &v)| v == max_votes)
        .map(|(&p, _)| p)
        .collect();
    top.sort_unstable();
    // Tie-break: prefer the reference's period, otherwise the smallest.
    let k = match p_ref {
        Some(pr) if top.contains(&pr) => pr,
        _ => top[0],
    };
    (k, "detected")
}

/// Return the lexicographically-smallest rotation of `kmer` together with all
/// of its rotations. Ports `get_rotations`.
pub fn get_rotations(kmer: &str) -> (String, Vec<String>) {
    let n = kmer.len();
    let mut rotations = Vec::with_capacity(n);
    for i in 0..n {
        let mut s = String::with_capacity(n);
        s.push_str(&kmer[i..]);
        s.push_str(&kmer[..i]);
        rotations.push(s);
    }
    let canonical = rotations.iter().min().cloned().unwrap_or_default();
    (canonical, rotations)
}

/// Collapse rotation-equivalent kmers onto their canonical form and normalise
/// counts to fractions of the total. Ports `prune_counts`.
pub fn prune_counts(counts: &HashMap<String, usize>) -> HashMap<String, f64> {
    let mut pruned: HashMap<String, usize> = HashMap::new();
    for key in counts.keys() {
        let (canonical, rotations) = get_rotations(key);
        if pruned.contains_key(&canonical) {
            continue;
        }
        let sum: usize = rotations
            .iter()
            .map(|r| counts.get(r).copied().unwrap_or(0))
            .sum();
        pruned.insert(canonical, sum);
    }
    let total: usize = pruned.values().sum();
    if total == 0 {
        return HashMap::new();
    }
    pruned
        .into_iter()
        .map(|(k, v)| (k, v as f64 / total as f64))
        .collect()
}

/// Count kmers of length `k` in `seq`, returning canonical fractions. Ports
/// `count_kmers`.
pub fn count_kmers(seq: &str, k: usize) -> HashMap<String, f64> {
    if k == 0 || seq.len() < k {
        return HashMap::new();
    }
    let mut counts: HashMap<String, usize> = HashMap::new();
    for i in 0..=(seq.len() - k) {
        *counts.entry(seq[i..i + k].to_string()).or_insert(0) += 1;
    }
    prune_counts(&counts)
}

/// Select the `number` most frequent kmers across all records. Ports
/// `select_kmers`.
///
/// For each record the canonical kmer fractions are computed, rounded to two
/// decimals (as the Python code does via `.round(2)`), and summed across all
/// records. Kmers are then ranked by total descending. Ties are broken by
/// first-appearance order (a stable sort), which can differ from pandas'
/// unstable quicksort for exact ties but is deterministic.
pub fn select_kmers(records: &[RepeatRecord], k: usize, number: usize) -> Vec<String> {
    let mut order: Vec<String> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    let mut totals: Vec<f64> = Vec::new();

    for rec in records {
        if rec.sequence.is_empty() {
            continue;
        }
        let kmers = count_kmers(&rec.sequence, k);
        if kmers.is_empty() {
            continue;
        }
        // Iterate in sorted key order for deterministic first-appearance order.
        let mut keys: Vec<&String> = kmers.keys().collect();
        keys.sort();
        for key in keys {
            let rounded = (kmers[key] * 100.0).round() / 100.0;
            match index.get(key) {
                Some(&i) => totals[i] += rounded,
                None => {
                    index.insert(key.clone(), order.len());
                    order.push(key.clone());
                    totals.push(rounded);
                }
            }
        }
    }

    let mut indices: Vec<usize> = (0..order.len()).collect();
    indices.sort_by(|&a, &b| {
        totals[b]
            .partial_cmp(&totals[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    indices
        .into_iter()
        .take(number)
        .map(|i| order[i].clone())
        .collect()
}

/// Pick the kmers to plot: user-specified motifs (sorted longest first to
/// accommodate overlap) or the most frequent kmers. Ports `get_kmers`.
///
/// With `-k auto` the kmer length is detected per locus from the data via
/// [`detect_locus_k`]; `coords` is used only for the log message.
pub fn get_kmers(records: &[RepeatRecord], args: &Args, coords: &str) -> Vec<String> {
    if let Some(motifs) = &args.motifs {
        let mut v: Vec<String> = motifs.split(',').map(|s| s.to_string()).collect();
        // Stable sort, longest first.
        v.sort_by_key(|m| std::cmp::Reverse(m.len()));
        v
    } else {
        let k = match args.kmer {
            KmerArg::Fixed(k) => k,
            KmerArg::Auto => {
                let (k, source) = detect_locus_k(records);
                let label = if coords.is_empty() { "table" } else { coords };
                eprintln!("Repeat {label}: auto-detected k={k} ({source})");
                k
            }
        };
        select_kmers(records, k, args.number)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotations_canonical() {
        let (canon, rots) = get_rotations("TGGCGC");
        assert_eq!(rots.len(), 6);
        assert_eq!(canon, "CGCTGG");
    }

    #[test]
    fn count_kmers_normalises() {
        // "AATAAT": kmers of length 3: AAT, ATA, TAA, AAT -> AAT x2, ATA, TAA
        // canonicalised: AAT(=AAT,ATA,TAA rotations)...
        let counts = count_kmers("AATAATAAT", 3);
        let total: f64 = counts.values().sum();
        assert!((total - 1.0).abs() < 1e-9);
    }

    #[test]
    fn detect_period_trinucleotide() {
        let seq = "CAGCAGCAGCAGCAGCAG".as_bytes();
        assert_eq!(detect_period(seq, 6), Some(3));
    }

    #[test]
    fn detect_period_pentamer() {
        let seq = "AAAATAAAATAAAATAAAATAAAAT".as_bytes();
        assert_eq!(detect_period(seq, 6), Some(5));
    }

    #[test]
    fn detect_period_arich_pentamer() {
        // AAAAG is ~80% A; without composition correction the abundance of A
        // lets short periods match by chance. Correction must still recover 5.
        let seq = "AAAAGAAAAGAAAAGAAAAGAAAAGAAAAG".as_bytes();
        assert_eq!(detect_period(seq, 6), Some(5));
    }

    #[test]
    fn detect_period_too_short() {
        assert_eq!(detect_period("CAG".as_bytes(), 6), None);
    }

    #[test]
    fn detect_period_tolerates_substitutions() {
        // A CAG repeat with a couple of point substitutions still scores above
        // the 0.65 threshold.
        let seq = "CAGCAGCTGCAGCAGCAGCGGCAG".as_bytes();
        assert_eq!(detect_period(seq, 6), Some(3));
    }

    fn rec(seq: &str, reference: &str) -> RepeatRecord {
        RepeatRecord {
            coords: "chr1:1".into(),
            sample: "s".into(),
            allele: "Allele1".into(),
            sequence: seq.into(),
            reference: reference.into(),
            case: false,
        }
    }

    #[test]
    fn locus_k_ref_and_median_agree() {
        let recs = vec![
            rec("CAGCAGCAGCAGCAGCAGCAG", "CAGCAGCAGCAGCAG"),
            rec("CAGCAGCAGCAGCAGCAGCAGCAGCAG", "CAGCAGCAGCAGCAG"),
        ];
        assert_eq!(detect_locus_k(&recs), (3, "detected"));
    }

    #[test]
    fn locus_k_no_ref_uses_median() {
        // No REF (e.g. --table input): detection rests on the median allele.
        let recs = vec![
            rec("AAAATAAAATAAAATAAAAT", ""),
            rec("AAAATAAAATAAAATAAAATAAAAT", ""),
        ];
        assert_eq!(detect_locus_k(&recs), (5, "detected"));
    }

    #[test]
    fn locus_k_empty_falls_back() {
        assert_eq!(detect_locus_k(&[]), (AUTO_K_DEFAULT, "fallback"));
    }
}
