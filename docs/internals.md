# aSTRonaut internals

Implementation notes for contributors and curious users. None of this is needed
to *use* aSTRonaut — see the [README](../README.md) for that.

## Motif-length detection (`-k auto`)

Ported from [trout](https://github.com/wdecoster/trout). For each locus it picks
the smallest period in 2–6 whose self-shift match rate (fraction of positions
where `seq[i] == seq[i + p]`) clears 0.65 — that is the fundamental period, since
multiples of the true period also score highly. It cross-checks the reference
allele (VCF `REF`, short but accurate) against the median-length allele (longer,
more informative): if both yield a period they must agree, if only one has signal
that one is trusted, otherwise it falls back to 3 (trinucleotide, the most common
pathogenic motif length). With `--table` input there is no `REF`, so detection
rests on the median allele alone. The chosen `k` is reported per locus on stderr.

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
