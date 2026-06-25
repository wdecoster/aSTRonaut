//! Command-line interface, mirroring `get_args()` in the Python tool.

use std::str::FromStr;

use clap::{Parser, ValueEnum};

/// Value of `-k`: either a fixed kmer length or `auto` (detect the repeat period
/// per locus from the data, ported from trout).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KmerArg {
    Auto,
    Fixed(usize),
}

impl FromStr for KmerArg {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.eq_ignore_ascii_case("auto") {
            Ok(KmerArg::Auto)
        } else {
            match s.parse::<usize>() {
                Ok(k) if k >= 1 => Ok(KmerArg::Fixed(k)),
                _ => Err(format!("expected a positive integer or 'auto', got '{s}'")),
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum LegendCorner {
    #[value(name = "topright")]
    TopRight,
    #[value(name = "bottomright")]
    BottomRight,
    #[value(name = "outside")]
    Outside,
}

/// How to order the sample/allele rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum SortMode {
    /// By sequence length (default).
    Length,
    /// Alphabetically by sample name.
    Alphabetic,
    /// Cluster by motif composition, so alleles dominated by the same motif
    /// group together (e.g. RFC1 AAGGG carriers).
    Motif,
}

#[derive(Debug, Parser)]
#[command(
    name = "aSTRonaut",
    version = env!("CARGO_PKG_VERSION"),
    about = "Create a repeat sequence plot similar to the pathSTR sequence composition visualization, but stand-alone"
)]
pub struct Args {
    /// VCF files to analyze
    #[arg(value_name = "vcf")]
    pub vcf: Vec<String>,

    /// Sample names to use, comma-separated
    #[arg(long)]
    pub names: Option<String>,

    /// Kmer length, or 'auto' to detect the repeat period per locus from the
    /// data. Ignored with --motifs.
    #[arg(short = 'k', long, default_value = "3")]
    pub kmer: KmerArg,

    /// Number of kmers to plot, ignored with --motifs
    #[arg(short = 'n', long, default_value_t = 10)]
    pub number: usize,

    /// Manually specify the motifs to plot, comma separated
    #[arg(long)]
    pub motifs: Option<String>,

    /// Chromosome and POS of repeat to plot from VCF e.g. chr1:12345, default: all repeats
    #[arg(long)]
    pub repeat: Option<String>,

    /// Output file name. Supported formats: html, svg, png, pdf
    #[arg(short = 'o', long, default_value = "astronaut.html")]
    pub out: String,

    /// Minimal allele length to plot
    #[arg(short = 'm', long, default_value_t = 20)]
    pub minlen: usize,

    /// Hide sample labels
    #[arg(long = "hide-labels", default_value_t = false)]
    pub hide_labels: bool,

    /// Size of sample labels
    #[arg(long, default_value_t = 8)]
    pub label_size: u32,

    /// Row ordering: length (default), alphabetic, or motif (cluster by composition)
    #[arg(long, value_enum, default_value = "length")]
    pub sort: SortMode,

    /// Collapse identical sequences into one row with an allele-count panel.
    /// Bare `--collapse` collapses only exact duplicates; `--collapse=0.05`
    /// also merges sequences within 5% normalized edit distance (to absorb
    /// sequencing error).
    #[arg(long, num_args = 0..=1, require_equals = true, default_missing_value = "0")]
    pub collapse: Option<f64>,

    /// Hide 'Allele' from labels
    #[arg(long, default_value_t = false)]
    pub hide_allele_label: bool,

    /// Create a plot suitable for publication
    #[arg(long, default_value_t = false)]
    pub publication: bool,

    /// Create a plot suitable for publication, but with less fluff
    #[arg(long, default_value_t = false)]
    pub minimal: bool,

    /// Title of the plot
    #[arg(long, default_value = "Repeat composition")]
    pub title: String,

    /// TSV file with sample information (name, group columns)
    #[arg(long)]
    pub sampleinfo: Option<String>,

    /// Parse somatic VCFs, expects a SEQS format field
    #[arg(long, default_value_t = false)]
    pub somatic: bool,

    /// Size of markers to plot
    #[arg(long, default_value_t = 3)]
    pub size: u32,

    /// Only plot the longest allele per individual
    #[arg(long, default_value_t = false)]
    pub longest_only: bool,

    /// Corner of the legend
    #[arg(long, value_enum, default_value = "bottomright")]
    pub legend_corner: LegendCorner,

    /// Height of the plot in pixels. Default: scales with the number of rows
    /// (so labels stay legible for large cohorts).
    #[arg(long)]
    pub height: Option<u32>,

    /// Width of the plot
    #[arg(long)]
    pub width: Option<u32>,

    /// Use a tsv (with name, sequence and group columns) as input
    #[arg(short = 't', long)]
    pub table: Option<String>,
}

impl Args {
    /// Validate mutually-exclusive and dependent options, mirroring the checks
    /// at the end of `get_args()`.
    pub fn validate(&self) -> Result<(), String> {
        if let Some(names) = &self.names {
            let n = names.split(',').count();
            if n != self.vcf.len() {
                return Err(format!(
                    "Number of names does not match number of VCFs\nNames: {:?}\nVCFs: {:?}",
                    names.split(',').collect::<Vec<_>>(),
                    self.vcf
                ));
            }
        }
        if self.table.is_some() && !self.vcf.is_empty() {
            return Err("Please provide either VCFs or a table, not both".to_string());
        }
        if self.table.is_none() && self.vcf.is_empty() {
            return Err("Please provide either VCFs or a table".to_string());
        }
        if self.publication && self.minimal {
            return Err("Please specify either --publication or --minimal, not both".to_string());
        }
        if self.longest_only && self.somatic {
            return Err("--longest_only is not supported with --somatic".to_string());
        }
        if let Some(t) = self.collapse {
            if !(0.0..1.0).contains(&t) {
                return Err("--collapse threshold must be between 0 and 1".to_string());
            }
        }
        Ok(())
    }
}
