# aSTRonaut

Rust implementation of the aSTRonaut tandem repeat visualization tool.

aSTRonaut creates a stand-alone "barcode" plot of the kmer composition of short
tandem repeats (STRs), similar to the pathSTR sequence composition visualization.
It reads STR genotyper VCFs (e.g. [STRdust](https://github.com/wdecoster/STRdust),
LongTR) or a TSV table, picks the most frequent kmers (collapsing rotationally
equivalent motifs), and draws one row per sample/allele with every nucleotide
coloured by the kmer that covers it.

Plotting is done with [kuva](https://psy-fer.github.io/kuva/); VCF parsing uses
[noodles](https://github.com/zaeleus/noodles). The result is a single static
binary with no Python/Plotly/kaleido runtime.

## Build

```bash
cargo build --release
# binary at target/release/aSTRonaut
```

## Usage

```bash
# One repeat from several VCFs, interactive HTML (default)
aSTRonaut sample1.vcf.gz sample2.vcf.gz --repeat chr1:57367043 -o out.html

# Let aSTRonaut detect the motif length per repeat
aSTRonaut *.vcf.gz -k auto -o plot.png

# Static image of a single repeat
aSTRonaut *.vcf.gz --repeat chr1:57367043 -o plot.png

# All repeats found in the VCFs (one file per repeat for image formats)
aSTRonaut sample.vcf.gz -o plot.svg

# From a TSV table with name, sequence and optional group columns
aSTRonaut --table sequences.tsv -o plot.png

# Publication styling, with case/control annotation
aSTRonaut *.vcf.gz --repeat chr1:57367043 \
    --sampleinfo info.tsv --publication --title "STR composition" -o fig.pdf
```

### Key options

- `-k/--kmer`, `-n/--number`: kmer length and how many kmers to colour.
  Pass `-k auto` to **detect the repeat period per locus** from the data
  (ported from [trout](https://github.com/wdecoster/trout)): for each repeat it
  picks the smallest period 2–6 whose self-shift match rate clears 0.65,
  cross-checking the reference allele (VCF REF) against the median-length allele
  and falling back to 3 when no clean period is found (with `--table` input,
  which has no REF, detection rests on the median allele alone). The chosen k is
  reported per locus on stderr.
- `--motifs`: colour specific motifs instead of the most frequent ones.
- `-m/--minlen`: minimum allele length to plot.
- `--repeat chr:pos`: restrict to one repeat coordinate.
- `--sort`: row ordering — `length` (default), `alphabetic`, or `motif`.
  **Motif clusters rows by composition**, so alleles dominated by the same motif
  group together (e.g. RFC1 `AAGGG` expansions cluster apart from the benign
  `AAAAG` alleles).
- `--collapse`: collapse identical sequences into one row per unique haplotype,
  with a marginal **allele-count histogram** beside the barcode (stacked
  case/control when `--sampleinfo` is given). Bare `--collapse` merges only exact
  duplicates; `--collapse=0.05` also merges sequences within 5% normalized edit
  distance (greedy, abundance-first, banded edit distance) to absorb sequencing
  error. Rows are ordered by count, most abundant at the top.
- `--somatic`: read per-read sequences from the `SEQS` INFO field.
- `--longest_only`: keep only the longest allele per individual.
- `--sampleinfo`: TSV (`name`, `group`) marking `case` samples with `>>>`.
- `--publication` / `--minimal`: publication-ready styling.
- `--hide-labels`, `--hide_allele_label`, `--label_size`,
  `--size`, `--legend_corner`, `--width`, `--title`.
- `--height`: plot height in pixels. If omitted, it **auto-scales with the
  number of rows** (~22 px/row, minimum 800) so labels stay legible for large
  cohorts; pass a value to override.

## Output formats

| Extension | Output |
|-----------|--------|
| `.html`   | Interactive SVG(s) embedded in one HTML page; all repeats in one file |
| `.svg`    | Static SVG; one file per repeat |
| `.png`    | Static PNG; one file per repeat |
| `.pdf`    | Static PDF; one file per repeat |

Each row is drawn as a **run-length-encoded bar** (consecutive same-motif
nucleotides collapse into one segment), so files stay small even for long
expansions — a 150-allele RFC1 plot drops from ~20 MB to <1 MB.

> **Note on interactivity:** the original Python tool produced an interactive
> Plotly HTML by default. kuva renders a self-contained interactive *SVG*
> instead, embedded in an HTML wrapper for the `.html` output (no JavaScript
> CDN). **Hover** a segment for `sample · motif · position`; use the **search
> box** (top-left) to filter — type a sample name to locate that individual's
> row, or a motif to highlight where it occurs. The legend is a static colour
> key. `jpeg`/`webp` outputs from the Python version are not supported; use
> `png` or `pdf`.

## Input

VCF (single sample, plain or BGZF-compressed) with `CHROM`, `POS`, `ALT` and a
`GT` genotype. Somatic mode additionally reads the `SEQS` (and optional
`OUTLIERS`) INFO fields. TSV input requires `name` and `sequence` columns and an
optional `group` column (where the value `case` drives the annotation).
