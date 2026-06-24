//! Core data model shared across the pipeline.

/// A single sequence observation: one allele (or read, in somatic mode) of a
/// repeat for a given sample. Mirrors a row of the pandas DataFrame in the
/// original `aSTRonaut.py`.
#[derive(Debug, Clone)]
pub struct RepeatRecord {
    /// Repeat coordinate `chr:pos` (empty string for `--table` input).
    pub coords: String,
    /// Sample name.
    pub sample: String,
    /// Allele label, e.g. `Allele1`, `Allele2`, `Allele1_3`, `Outlier_0`.
    pub allele: String,
    /// The repeat sequence. Only records with a sequence are kept.
    pub sequence: String,
    /// Reference allele sequence at this locus (VCF REF column). Same for every
    /// record of a locus; empty for `--table` input. Used by `-k auto`.
    pub reference: String,
    /// Whether this sample is a `case` (drives the `>>>` annotation).
    pub case: bool,
}

impl RepeatRecord {
    /// Identifier used on the y-axis: `sample_allele`, or just `sample` when
    /// allele labels are hidden.
    pub fn identifier(&self, hide_allele_label: bool) -> String {
        if hide_allele_label {
            self.sample.clone()
        } else {
            format!("{}_{}", self.sample, self.allele)
        }
    }
}
