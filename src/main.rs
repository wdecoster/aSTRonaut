//! aSTRonaut: stand-alone tandem-repeat sequence composition plots.
//! Rust port of `aSTRonaut.py`.

mod cli;
mod collapse;
mod input;
mod kmer;
mod model;
mod palette;
mod plot;

use std::collections::HashSet;
use std::fs;
use std::path::Path;

use anyhow::{anyhow, bail, Result};
use clap::Parser;
use kuva::prelude::*;

use cli::Args;
use model::RepeatRecord;

fn main() {
    let args = Args::parse();
    if let Err(e) = run(&args) {
        // `{:#}` prints the full anyhow context chain, e.g.
        // "opening sample.vcf.gz: No such file or directory".
        eprintln!("ERROR: {e:#}");
        std::process::exit(1);
    }
}

fn run(args: &Args) -> Result<()> {
    args.validate().map_err(|e| anyhow!(e))?;

    let records = input::parse_input(args)?;
    if records.is_empty() {
        bail!("No sequences to plot");
    }

    // Unique repeat coordinates in first-appearance order.
    let mut coords_order: Vec<String> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    for r in &records {
        if seen.insert(r.coords.clone()) {
            coords_order.push(r.coords.clone());
        }
    }

    let ext = Path::new(&args.out)
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();

    match ext.as_str() {
        "html" => {
            // Interactive SVGs concatenated into a single HTML document.
            let mut body = String::new();
            for coords in &coords_order {
                let scene = build_scene_for(&records, coords, args, true);
                let svg = SvgBackend.render_scene(&scene);
                body.push_str(&finalize_interactive_svg(svg));
                body.push('\n');
            }
            let help = if args.collapse.is_some() {
                "Hover a segment for <b>count · motif · position</b>. \
Use the <b>search</b> box (bottom-left) to highlight a motif (e.g. \
<code>AAGGG</code>). The colour legend is a static key."
            } else {
                "Hover a segment for <b>sample · motif · position</b>. \
Use the <b>search</b> box (bottom-left) to filter rows: type a sample name \
(e.g. <code>HG00096</code>) to find that individual, or a motif (e.g. \
<code>AAGGG</code>) to highlight where it occurs. The colour legend is a static key."
            };
            let html = format!(
                "<!DOCTYPE html>\n<html>\n<head><meta charset=\"utf-8\"><title>aSTRonaut</title>\n\
<style>body{{font-family:DejaVu Sans,Arial,sans-serif}}.astro-help{{color:#444;font-size:13px;margin:8px 12px;max-width:60em}}</style>\n\
</head>\n<body>\n\
<p class=\"astro-help\">{help}</p>\n\
{body}</body>\n</html>\n"
            );
            fs::write(&args.out, html)?;
        }
        "svg" | "png" | "pdf" => {
            let multi = coords_order.len() > 1;
            for coords in &coords_order {
                let out_path = if multi {
                    let safe = coords.replace(':', "_");
                    let base = args
                        .out
                        .strip_suffix(&format!(".{ext}"))
                        .unwrap_or(&args.out);
                    format!("{base}_{safe}.{ext}")
                } else {
                    args.out.clone()
                };
                write_static(&records, coords, args, &ext, &out_path)?;
            }
        }
        other => bail!("Unsupported output format '.{other}'. Use html, svg, png, or pdf"),
    }

    Ok(())
}

/// Vertical space (px) reserved at the bottom of interactive SVGs for the UI
/// strip, matching the amount kuva reserves only when the height is auto.
const UI_STRIP_PAD: f64 = 32.0;

/// Post-process an interactive SVG to fix two kuva quirks we hit because we set
/// an explicit (row-scaled) height:
///
/// 1. **Inert legend** — kuva wires click-to-toggle onto `g.legend-entry`, but
///    our run-length bar segments carry no per-motif group, so clicking only
///    greyed the legend without changing the plot. Renaming the class drops the
///    click handler and the pointer cursor, making the legend a static key.
/// 2. **UI strip overlap** — kuva only reserves bottom space for the
///    search/save/info strip when the height is auto-computed. We grow each
///    canvas by [`UI_STRIP_PAD`] and slide the strip into that new band so it
///    sits clear of the x-axis label and ticks.
fn finalize_interactive_svg(svg: String) -> String {
    let svg = svg.replace("class=\"legend-entry\"", "class=\"legend-key\"");
    // Wrap the UI strip (search foreignObject … info icon) in a downward shift.
    let svg = svg.replace(
        "<foreignObject",
        &format!("<g class=\"astro-ui\" transform=\"translate(0,{UI_STRIP_PAD})\"><foreignObject"),
    );
    let svg = svg.replace("<text id=\"kuva-readout\"", "</g><text id=\"kuva-readout\"");
    grow_svg_heights(svg, UI_STRIP_PAD)
}

/// Add `pad` to the `height` attribute of every root `<svg>` opening tag,
/// leaving inner `height=` attributes (foreignObject, rects) untouched.
fn grow_svg_heights(svg: String, pad: f64) -> String {
    let mut out = String::with_capacity(svg.len() + 64);
    let mut rest = svg.as_str();
    while let Some(pos) = rest.find("<svg") {
        out.push_str(&rest[..pos]);
        let after = &rest[pos..];
        let Some(tag_end) = after.find('>') else {
            out.push_str(after);
            return out;
        };
        let tag = &after[..tag_end];
        if let Some(hpos) = tag.find("height=\"") {
            let vstart = hpos + "height=\"".len();
            if let Some(vlen) = tag[vstart..].find('"') {
                if let Ok(num) = tag[vstart..vstart + vlen].parse::<f64>() {
                    out.push_str(&tag[..hpos]);
                    out.push_str(&format!("height=\"{}\"", num + pad));
                    out.push_str(&tag[vstart + vlen + 1..]);
                    out.push('>');
                    rest = &after[tag_end + 1..];
                    continue;
                }
            }
        }
        out.push_str(&after[..tag_end + 1]);
        rest = &after[tag_end + 1..];
    }
    out.push_str(rest);
    out
}

/// Build the renderable scene for one repeat coordinate.
fn build_scene_for(
    records: &[RepeatRecord],
    coords: &str,
    args: &Args,
    interactive: bool,
) -> kuva::render::render::Scene {
    let repeat_records: Vec<RepeatRecord> = records
        .iter()
        .filter(|r| r.coords == coords)
        .cloned()
        .collect();
    let kmers = kmer::get_kmers(&repeat_records, args, coords);
    plot::build_scene(&repeat_records, &kmers, coords, args, interactive)
}

/// Render one repeat to a static file in the given format.
fn write_static(
    records: &[RepeatRecord],
    coords: &str,
    args: &Args,
    ext: &str,
    path: &str,
) -> Result<()> {
    let scene = build_scene_for(records, coords, args, false);
    let bytes = match ext {
        "svg" => SvgBackend.render_scene(&scene).into_bytes(),
        "png" => PngBackend::new()
            .with_scale(2.0)
            .render_scene(&scene)
            .map_err(|e| anyhow!(e))?,
        "pdf" => PdfBackend::new()
            .render_scene(&scene)
            .map_err(|e| anyhow!(e))?,
        _ => unreachable!(),
    };
    fs::write(path, bytes)?;
    Ok(())
}
