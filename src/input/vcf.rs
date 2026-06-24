//! VCF parsing with noodles. Ports `parse_vcf` and `parse_alts`.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use noodles_vcf::variant::record::samples::series::Value as SampleValue;
use noodles_vcf::variant::record::AlternateBases as _;
use noodles_vcf::{self as vcf};

use crate::cli::Args;
use crate::model::RepeatRecord;

/// Open a possibly gzip/BGZF-compressed VCF as a line reader. BGZF is a
/// concatenation of gzip members, which `MultiGzDecoder` handles transparently.
fn open_reader(path: &str) -> Result<Box<dyn BufRead>> {
    let mut magic = [0u8; 2];
    {
        let mut probe = File::open(path).with_context(|| format!("opening {path}"))?;
        let _ = probe.read(&mut magic)?;
    }
    let file = File::open(path).with_context(|| format!("opening {path}"))?;
    if magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Map a phased genotype to the two allele sequences. Ports `parse_alts`.
/// Position `None` (missing) or `Some(0)` (reference) yields `None`; otherwise
/// the corresponding alternate allele.
fn allele_seq(position: Option<usize>, alts: &[String]) -> Option<String> {
    match position {
        None | Some(0) => None,
        Some(n) => alts.get(n - 1).cloned(),
    }
}

/// Default sample name from a VCF path, stripping `.vcf.gz` / `.vcf`.
fn default_name(path: &str) -> String {
    let base = Path::new(path)
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| path.to_string());
    base.strip_suffix(".vcf.gz")
        .or_else(|| base.strip_suffix(".vcf"))
        .unwrap_or(&base)
        .to_string()
}

/// Parse a single VCF file into repeat records. Ports `parse_vcf`.
pub fn parse_vcf(path: &str, args: &Args, name: Option<&str>) -> Result<Vec<RepeatRecord>> {
    let name = name
        .map(|s| s.to_string())
        .unwrap_or_else(|| default_name(path));

    let mut reader = vcf::io::Reader::new(open_reader(path)?);
    let header = reader
        .read_header()
        .with_context(|| format!("reading header of {path}"))?;

    let mut calls: Vec<RepeatRecord> = Vec::new();

    for result in reader.records() {
        let record = result.with_context(|| format!("reading a record from {path}"))?;

        let chrom = record.reference_sequence_name();
        let pos = match record.variant_start() {
            Some(p) => usize::from(p?),
            None => continue,
        };
        let coords = format!("{chrom}:{pos}");

        if let Some(target) = &args.repeat {
            if &coords != target {
                continue;
            }
        }

        // Reference allele sequence (shared by every record at this locus).
        let reference = record.reference_bases().to_string();

        // Alternate alleles.
        let mut alts: Vec<String> = Vec::new();
        for alt in record.alternate_bases().iter() {
            alts.push(alt?.to_string());
        }

        // Genotype of the (single) sample: first two allele positions.
        let mut positions: [Option<usize>; 2] = [None, None];
        if let Some(sample) = record.samples().get_index(0) {
            for field in sample.iter(&header) {
                let (key, value) = field?;
                if key == "GT" {
                    if let Some(SampleValue::Genotype(genotype)) = value {
                        for (i, allele) in genotype.iter().take(2).enumerate() {
                            positions[i] = allele?.0;
                        }
                    }
                    break;
                }
            }
        }

        let sequences = [
            allele_seq(positions[0], &alts),
            allele_seq(positions[1], &alts),
        ];

        // Somatic read-level sequences from the SEQS INFO field.
        let somatic_sequences: Vec<String> = if args.somatic {
            match info_string(&record, &header, "SEQS")? {
                Some(s) => s.split(',').map(|x| x.to_string()).collect(),
                None => {
                    eprintln!("WARNING: No SEQS field found in VCF, skipping");
                    return Ok(Vec::new());
                }
            }
        } else {
            Vec::new()
        };

        for (allele_idx, seq) in sequences.iter().enumerate() {
            let Some(seq) = seq else { continue };
            if seq.len() <= args.minlen {
                continue;
            }
            let label = if allele_idx == 0 {
                "Allele1"
            } else {
                "Allele2"
            };
            if args.somatic {
                if let Some(group) = somatic_sequences.get(allele_idx) {
                    for (i, s) in group.split(':').enumerate() {
                        if s.len() > args.minlen {
                            calls.push(RepeatRecord {
                                coords: coords.clone(),
                                sample: name.clone(),
                                allele: format!("{label}_{i}"),
                                sequence: s.to_string(),
                                reference: reference.clone(),
                                case: false,
                            });
                        }
                    }
                }
                // Outliers are added only once (with the first allele), as in Python.
                if allele_idx == 0 {
                    if let Some(outliers) = info_string(&record, &header, "OUTLIERS")? {
                        for (i, s) in outliers.split(',').enumerate() {
                            if s.len() > args.minlen {
                                calls.push(RepeatRecord {
                                    coords: coords.clone(),
                                    sample: name.clone(),
                                    allele: format!("Outlier_{i}"),
                                    sequence: s.to_string(),
                                    reference: reference.clone(),
                                    case: false,
                                });
                            }
                        }
                    }
                }
            } else {
                calls.push(RepeatRecord {
                    coords: coords.clone(),
                    sample: name.clone(),
                    allele: label.to_string(),
                    sequence: seq.clone(),
                    reference: reference.clone(),
                    case: false,
                });
            }
        }
    }

    if args.longest_only && !calls.is_empty() {
        let longest = calls
            .into_iter()
            .max_by_key(|c| c.sequence.len())
            .expect("non-empty");
        return Ok(vec![longest]);
    }

    Ok(calls)
}

/// Fetch a string-typed INFO field value, if present.
fn info_string(record: &vcf::Record, header: &vcf::Header, key: &str) -> Result<Option<String>> {
    use noodles_vcf::variant::record::info::field::Value as InfoValue;
    match record.info().get(header, key) {
        Some(result) => match result? {
            Some(InfoValue::String(s)) => Ok(Some(s.into_owned())),
            _ => Ok(None),
        },
        None => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use std::io::Write;

    fn write_vcf(tag: &str, body: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("astronaut_test_{}_{}.vcf", std::process::id(), tag));
        let header = "##fileformat=VCFv4.2\n\
##INFO=<ID=SEQS,Number=1,Type=String,Description=\"\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(header.as_bytes()).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        path
    }

    fn args() -> Args {
        Args::parse_from(["aSTRonaut", "dummy.vcf"])
    }

    #[test]
    fn parses_two_alleles() {
        // ALT1 length 30, ALT2 length 10; minlen default 20 -> only Allele1 kept.
        let alt1 = "A".repeat(30);
        let alt2 = "C".repeat(10);
        let path = write_vcf(
            "two_alleles",
            &format!("chr1\t100\t.\tT\t{alt1},{alt2}\t.\t.\t.\tGT\t1|2\n"),
        );
        let records = parse_vcf(path.to_str().unwrap(), &args(), Some("S1")).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].coords, "chr1:100");
        assert_eq!(records[0].sample, "S1");
        assert_eq!(records[0].allele, "Allele1");
        assert_eq!(records[0].sequence.len(), 30);
    }

    #[test]
    fn reference_allele_is_skipped() {
        // Genotype 0|1 -> first phase is reference (skipped), second is ALT1.
        let alt1 = "G".repeat(25);
        let path = write_vcf(
            "ref_skip",
            &format!("chr2\t5\t.\tT\t{alt1}\t.\t.\t.\tGT\t0|1\n"),
        );
        let records = parse_vcf(path.to_str().unwrap(), &args(), Some("S2")).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].allele, "Allele2");
        assert_eq!(records[0].coords, "chr2:5");
    }
}
