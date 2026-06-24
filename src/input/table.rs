//! TSV table input. Ports the `--table` branch of `parse_input`.

use std::fs;

use anyhow::{bail, Result};

use crate::cli::Args;
use crate::model::RepeatRecord;

/// Split a tab-separated line into trimmed fields.
fn fields(line: &str) -> Vec<&str> {
    line.split('\t').map(|s| s.trim_end_matches('\r')).collect()
}

/// Parse a TSV with `name`, `sequence` and optional `group` columns.
pub fn parse_table(path: &str, args: &Args) -> Result<Vec<RepeatRecord>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();

    let header = match lines.next() {
        Some(h) => fields(h),
        None => bail!("Empty table"),
    };
    let name_idx = header
        .iter()
        .position(|&c| c == "name")
        .ok_or_else(|| anyhow::anyhow!("No 'name' column found in table"))?;
    let seq_idx = header
        .iter()
        .position(|&c| c == "sequence")
        .ok_or_else(|| anyhow::anyhow!("No 'sequence' column found in table"))?;
    let group_idx = header.iter().position(|&c| c == "group");

    // Collect rows first so we can validate the group column like the Python code.
    let mut rows: Vec<(String, String, Option<String>)> = Vec::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let cols = fields(line);
        let name = cols.get(name_idx).copied().unwrap_or("").to_string();
        let sequence = cols.get(seq_idx).copied().unwrap_or("").to_string();
        let group = group_idx.and_then(|i| cols.get(i).map(|s| s.to_string()));
        rows.push((name, sequence, group));
    }

    if group_idx.is_some() {
        let has_case = rows.iter().any(|(_, _, g)| g.as_deref() == Some("case"));
        if !has_case {
            bail!("'case' is not a value in the 'group' column of the table!");
        }
    }

    let mut records = Vec::new();
    for (name, sequence, group) in rows {
        if sequence.len() > args.minlen {
            let case = group.as_deref() == Some("case");
            records.push(RepeatRecord {
                coords: String::new(),
                sample: name,
                allele: "Allele1".to_string(),
                sequence,
                reference: String::new(),
                case,
            });
        }
    }
    Ok(records)
}
