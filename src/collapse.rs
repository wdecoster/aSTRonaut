//! Collapse identical (or near-identical) repeat sequences into clusters with
//! allele counts, for the `--collapse` view.

use std::collections::HashMap;

use crate::model::RepeatRecord;

/// A group of collapsed alleles sharing one (representative) sequence.
#[derive(Debug, Clone)]
pub struct Cluster {
    /// Sequence drawn for this row (the most abundant exact member).
    pub representative: String,
    /// Total alleles in the cluster.
    pub count: usize,
    /// Alleles from `case` samples.
    pub case_count: usize,
    /// Alleles from non-case samples.
    pub control_count: usize,
}

/// Collapse records into clusters. With `threshold <= 0` only exact duplicates
/// are merged; otherwise sequences within `threshold` normalized edit distance
/// are merged greedily, most-abundant-first (the abundant sequence becomes the
/// representative). Result is sorted by count descending.
pub fn collapse(records: &[&RepeatRecord], threshold: f64) -> Vec<Cluster> {
    // Exact collapse: sequence -> (count, case_count).
    let mut map: HashMap<&str, (usize, usize)> = HashMap::new();
    for r in records {
        if r.sequence.is_empty() {
            continue;
        }
        let e = map.entry(r.sequence.as_str()).or_insert((0, 0));
        e.0 += 1;
        if r.case {
            e.1 += 1;
        }
    }
    let mut clusters: Vec<Cluster> = map
        .into_iter()
        .map(|(seq, (count, case))| Cluster {
            representative: seq.to_string(),
            count,
            case_count: case,
            control_count: count - case,
        })
        .collect();
    sort_clusters(&mut clusters);

    if threshold <= 0.0 {
        return clusters;
    }

    // Fuzzy greedy merge over the (already abundance-sorted) unique sequences.
    let mut merged: Vec<Cluster> = Vec::new();
    for c in clusters {
        let cb = c.representative.as_bytes();
        let mut joined = false;
        for m in merged.iter_mut() {
            let mb = m.representative.as_bytes();
            let maxlen = cb.len().max(mb.len());
            let k = (threshold * maxlen as f64).floor() as usize;
            if cb.len().abs_diff(mb.len()) > k {
                continue; // length difference alone exceeds the budget
            }
            if within_distance(cb, mb, k) {
                m.count += c.count;
                m.case_count += c.case_count;
                m.control_count += c.control_count;
                joined = true;
                break;
            }
        }
        if !joined {
            merged.push(c);
        }
    }
    sort_clusters(&mut merged);
    merged
}

/// Sort by count descending, breaking ties on the representative for
/// determinism.
fn sort_clusters(clusters: &mut [Cluster]) {
    clusters.sort_by(|a, b| {
        b.count
            .cmp(&a.count)
            .then_with(|| a.representative.cmp(&b.representative))
    });
}

/// Banded Levenshtein with cutoff: returns `true` iff the edit distance between
/// `a` and `b` is at most `k`. Runs in O(n·k) by only filling the diagonal band
/// of width `2k+1` and bailing out once an entire row exceeds `k`.
fn within_distance(a: &[u8], b: &[u8], k: usize) -> bool {
    let (n, m) = (a.len(), b.len());
    if n.abs_diff(m) > k {
        return false;
    }
    const INF: usize = usize::MAX / 2;
    let mut prev = vec![INF; m + 1];
    for (j, slot) in prev.iter_mut().enumerate().take(k.min(m) + 1) {
        *slot = j;
    }
    for i in 1..=n {
        let mut cur = vec![INF; m + 1];
        cur[0] = if i <= k { i } else { INF };
        let lo = i.saturating_sub(k).max(1);
        let hi = (i + k).min(m);
        let mut row_min = INF;
        for j in lo..=hi {
            let cost = usize::from(a[i - 1] != b[j - 1]);
            let v = (prev[j].saturating_add(1))
                .min(cur[j - 1].saturating_add(1))
                .min(prev[j - 1].saturating_add(cost));
            cur[j] = v;
            row_min = row_min.min(v);
        }
        if row_min > k {
            return false;
        }
        prev = cur;
    }
    prev[m] <= k
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(seq: &str, case: bool) -> RepeatRecord {
        RepeatRecord {
            coords: "c".into(),
            sample: "s".into(),
            allele: "Allele1".into(),
            sequence: seq.into(),
            reference: String::new(),
            case,
        }
    }

    #[test]
    fn exact_collapse_counts() {
        let r = [
            rec("AAAAT", false),
            rec("AAAAT", true),
            rec("AAAAT", false),
            rec("CCCCG", false),
        ];
        let refs: Vec<&RepeatRecord> = r.iter().collect();
        let clusters = collapse(&refs, 0.0);
        assert_eq!(clusters.len(), 2);
        // Most abundant first.
        assert_eq!(clusters[0].representative, "AAAAT");
        assert_eq!(clusters[0].count, 3);
        assert_eq!(clusters[0].case_count, 1);
        assert_eq!(clusters[0].control_count, 2);
        assert_eq!(clusters[1].count, 1);
    }

    #[test]
    fn within_distance_basics() {
        assert!(within_distance(b"AAAAA", b"AAAAA", 0));
        assert!(within_distance(b"AAAAA", b"AAATA", 1));
        assert!(!within_distance(b"AAAAA", b"AATTA", 1));
        assert!(within_distance(b"AAAAA", b"AAAA", 1)); // one deletion
        assert!(!within_distance(b"AAAAA", b"AAA", 1)); // two deletions
    }

    #[test]
    fn fuzzy_merges_near_duplicates() {
        // Two abundant exact sequences plus a single-error variant of the first.
        let mut r = vec![rec("AAAAAAAAAAAAAAAAAAAA", false); 5]; // 20 A's, count 5
        r.push(rec("AAAAAAAAAATAAAAAAAAA", false)); // one substitution, count 1
        r.extend(vec![rec("CCCCCCCCCCCCCCCCCCCC", false); 3]); // distinct, count 3
        let refs: Vec<&RepeatRecord> = r.iter().collect();

        // Exact: 3 clusters.
        assert_eq!(collapse(&refs, 0.0).len(), 3);

        // 5% of 20 = 1 edit allowed -> the variant merges into the A-cluster.
        let fuzzy = collapse(&refs, 0.05);
        assert_eq!(fuzzy.len(), 2);
        assert_eq!(fuzzy[0].representative, "AAAAAAAAAAAAAAAAAAAA");
        assert_eq!(fuzzy[0].count, 6);
    }
}
