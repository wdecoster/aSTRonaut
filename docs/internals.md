# aSTRonaut internals

Implementation notes for contributors and curious users. None of this is needed
to *use* aSTRonaut — see the [README](../README.md) for that.

## Motif-length detection (`-k auto`)

Adapted from [trout](https://github.com/wdecoster/trout). For a single sequence,
`detect_period` measures the self-shift match rate (fraction of positions where
`seq[i] == seq[i + p]`) for each candidate period 2–6 and returns the smallest
one whose **composition-corrected** score clears a threshold — that is the
fundamental period, since multiples of the true period also score highly.

The composition correction matters: the raw match rate is inflated by base bias
(an A-rich motif like RFC1's `AAAAG` matches itself often just because most bases
are A), which let short periods pass spuriously. So instead of the raw rate we
score it above chance:

```
score = (observed − chance) / (1 − chance),   chance = Σ frequencyₐ²
```

`chance` is the probability two random bases of the sequence are equal, so
composition cancels out: a pure repeat still scores ~1.0, GC-rich hexamers
(C9orf72 `GGCCCC`) still clear it, and A-rich motifs no longer match at short
periods.

For a whole locus, `detect_locus_k` takes a **consensus** rather than trusting a
single allele: it detects the period of the reference allele (VCF `REF`) and of
every allele, and returns the most-voted period (ties prefer the reference's
period, then the smaller one). This is robust to a few atypical alleles, so the
result no longer depends on which cohort subset you happen to plot. It falls back
to 3 (the most common pathogenic motif length) only when nothing yields a period.
With `--table` input there is no `REF`, so the vote is over the alleles alone. The
chosen `k` is reported per locus on stderr.

It is still a heuristic and can occasionally be wrong; pass an explicit `-k` to
override.

## Collapsing (`--collapse`)

Exact collapse hashes sequences into unique haplotypes with allele (and, with
`--sampleinfo`, case/control) counts. The fuzzy variant (`--collapse=0.05`) then
merges near-identical sequences greedily, most-abundant-first, using a length
pre-filter plus a banded (Ukkonen) edit distance with a cutoff at
`threshold × length`. Working on the *unique* set keeps it fast even for a full
cohort. Normalised distance is `edits / max(len_a, len_b)`, so length differences
count — which is why a small threshold mostly absorbs sequencing error rather than
merging genuine length variants.

## Run-length rendering

Each row is a horizontal stacked bar where consecutive nucleotides of the same
motif collapse into one segment, instead of one marker per nucleotide. This keeps
SVG/HTML files small even for kilobase expansions (a 150-allele RFC1 plot drops
from ~20 MB to under 1 MB) and renders faster. The trade-off: the legend is a
static colour key rather than a clickable filter — use the search box for
filtering instead.

## Interactive HTML

`.html` output embeds a self-contained interactive SVG produced by
[kuva](https://psy-fer.github.io/kuva/) (no external JavaScript). Hovering a
segment shows `sample · motif · position`; the search box matches that text, so
typing a sample name locates an individual and typing a motif highlights it.
Because we set an explicit, row-scaled canvas height, kuva skips its automatic
bottom-margin reservation for the search strip, so the generated SVG is
post-processed to add that space and to make the legend non-interactive.

## Relationship to the Python tool

This is a port of the original Python `aSTRonaut.py` (part of pathSTR). The Rust
version is a single static binary with no Python/Plotly/kaleido runtime. The main
visible differences: the default `.html` output is an interactive SVG embedded in
HTML rather than a Plotly page, and `jpeg`/`webp` outputs are dropped (use `png`
or `pdf`).
