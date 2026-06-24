//! Input dispatch: VCF or TSV table, plus optional sample-info annotation.
//! Ports `parse_input` and the `--sampleinfo` handling from `main`.

mod table;
mod vcf;

use std::collections::HashMap;
use std::fs;

use anyhow::{bail, Result};

use crate::cli::Args;
use crate::model::RepeatRecord;

/// Parse all input into repeat records and apply sample-info case annotation.
pub fn parse_input(args: &Args) -> Result<Vec<RepeatRecord>> {
    let mut records = if let Some(table) = &args.table {
        table::parse_table(table, args)?
    } else {
        let names: Option<Vec<String>> = args
            .names
            .as_ref()
            .map(|n| n.split(',').map(|s| s.to_string()).collect());
        let mut all = Vec::new();
        for (i, path) in args.vcf.iter().enumerate() {
            let name = names.as_ref().map(|v| v[i].as_str());
            all.extend(vcf::parse_vcf(path, args, name)?);
        }
        all
    };

    if let Some(sampleinfo) = &args.sampleinfo {
        let cases = load_sampleinfo(sampleinfo)?;
        for rec in &mut records {
            rec.case = cases.get(&rec.sample).copied().unwrap_or(false);
        }
    }

    Ok(records)
}

/// Load `name -> is_case` from a sample-info TSV (columns `name`, `group`).
/// Errors if no `case` value is present, mirroring the Python check.
fn load_sampleinfo(path: &str) -> Result<HashMap<String, bool>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let header: Vec<&str> = match lines.next() {
        Some(h) => h.split('\t').map(|s| s.trim_end_matches('\r')).collect(),
        None => bail!("Empty sample info file"),
    };
    let name_idx = header
        .iter()
        .position(|&c| c == "name")
        .ok_or_else(|| anyhow::anyhow!("No 'name' column in sample info file"))?;
    let group_idx = header
        .iter()
        .position(|&c| c == "group")
        .ok_or_else(|| anyhow::anyhow!("No 'group' column in sample info file"))?;

    let mut map = HashMap::new();
    let mut saw_case = false;
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').map(|s| s.trim_end_matches('\r')).collect();
        let name = cols.get(name_idx).copied().unwrap_or("").to_string();
        let group = cols.get(group_idx).copied().unwrap_or("");
        let is_case = group == "case";
        saw_case |= is_case;
        map.insert(name, is_case);
    }
    if !saw_case {
        bail!("'case' is not a value in the 'group' column of the sample info file!");
    }
    Ok(map)
}
