//! Plot construction with kuva. Ports `plot_sequence`.

use std::collections::HashMap;

use kuva::prelude::*;
use kuva::render::render::Scene;

use crate::cli::{Args, LegendCorner, SortMode};
use crate::collapse::{collapse, Cluster};
use crate::model::RepeatRecord;
use crate::palette::{colors_for, GRAY};

/// Case / control fill colours for the collapsed count panel.
const CASE_COLOR: &str = "rgb(204, 102, 119)";
const CONTROL_COLOR: &str = "rgb(68, 119, 170)";

/// Map the legend-corner option to a kuva legend position.
fn legend_position(corner: LegendCorner) -> LegendPosition {
    match corner {
        LegendCorner::TopRight => LegendPosition::InsideTopRight,
        LegendCorner::BottomRight => LegendPosition::InsideBottomRight,
        LegendCorner::Outside => LegendPosition::OutsideRightTop,
    }
}

/// Apply the shared publication / minimal styling to a layout.
fn apply_style(mut layout: Layout, args: &Args) -> Layout {
    if args.publication || args.minimal {
        let mut theme = Theme::light();
        theme.background = "none".to_string();
        layout = layout.with_theme(theme).with_scale(1.5);
        if args.publication {
            layout = layout.with_box_axes();
        }
    }
    layout
}

/// Per-row pixel height used for the auto height calculation.
fn per_row_px(args: &Args) -> f64 {
    if args.hide_labels {
        (args.size as f64 * 2.0 + 4.0).max(8.0)
    } else {
        22.0
    }
}

/// Per-allele fraction of each selected motif (in `kmers` order), used to
/// cluster rows by composition in `SortMode::Motif`. Uses the first motif's
/// length as k (exact when all selected motifs share a length, as they do for
/// fixed/auto k; approximate for mixed-length `--motifs`).
fn motif_key(seq: &str, kmers: &[&str]) -> Vec<f64> {
    if kmers.is_empty() {
        return Vec::new();
    }
    let k = kmers[0].len();
    let fracs = crate::kmer::count_kmers(seq, k);
    kmers
        .iter()
        .map(|m| {
            let canon = crate::kmer::get_rotations(m).0;
            fracs.get(&canon).copied().unwrap_or(0.0)
        })
        .collect()
}

/// Order two sequences by motif composition: descending fraction of the first
/// motif, then the second, and so on (the same ordering used by
/// `SortMode::Motif`).
fn motif_cmp(a: &str, b: &str, kmers: &[&str]) -> std::cmp::Ordering {
    let ka = motif_key(a, kmers);
    let kb = motif_key(b, kmers);
    for (x, y) in ka.iter().zip(kb.iter()) {
        match y.partial_cmp(x).unwrap_or(std::cmp::Ordering::Equal) {
            std::cmp::Ordering::Equal => continue,
            o => return o,
        }
    }
    std::cmp::Ordering::Equal
}

/// Color every nucleotide of `seq` by the kmer that covers it, replicating the
/// successive string-replacement semantics of the Python code. Returns one
/// color token per nucleotide.
fn color_sequence(seq: &str, kmer_colors: &[(String, String)]) -> Vec<String> {
    let mut s = seq.to_string();
    for (kmer, color) in kmer_colors {
        let replacement = format!("{};", color).repeat(kmer.len());
        s = s.replace(kmer.as_str(), &replacement);
    }
    // Replace any remaining bare nucleotides with the gray ("other") token.
    let gray_token = format!("{};", GRAY);
    for nuc in ['A', 'C', 'T', 'G'] {
        s = s.replace(nuc, &gray_token);
    }
    let trimmed = s.trim_end_matches(';');
    if trimmed.is_empty() {
        Vec::new()
    } else {
        trimmed.split(';').map(|x| x.to_string()).collect()
    }
}

/// Run-length encode a sequence's per-nucleotide colors into
/// `(length, color, motif_label)` segments.
fn run_length(
    seq: &str,
    kmer_colors: &[(String, String)],
    color_to_kmer: &HashMap<String, String>,
) -> Vec<(usize, String, String)> {
    let mut runs: Vec<(usize, String, String)> = Vec::new();
    for color in color_sequence(seq, kmer_colors) {
        match runs.last_mut() {
            Some((len, c, _)) if *c == color => *len += 1,
            _ => {
                let motif = color_to_kmer
                    .get(&color)
                    .cloned()
                    .unwrap_or_else(|| "other".to_string());
                runs.push((1, color, motif));
            }
        }
    }
    runs
}

/// Build the renderable scene for a single repeat: the normal per-allele
/// barcode, or — under `--collapse` — a two-panel figure (barcode of unique
/// sequences + a case/control allele-count histogram).
pub fn build_scene(
    records: &[RepeatRecord],
    kmers: &[String],
    repeat: &str,
    args: &Args,
    interactive: bool,
) -> Scene {
    if let Some(threshold) = args.collapse {
        build_collapsed_scene(records, kmers, repeat, args, interactive, threshold)
    } else {
        let (plots, layout) = build_plot(records, kmers, repeat, args, interactive);
        render_multiple(plots, layout)
    }
}

/// Build the kuva plots and layout for a single repeat. Ports `plot_sequence`.
fn build_plot(
    records: &[RepeatRecord],
    kmers: &[String],
    repeat: &str,
    args: &Args,
    interactive: bool,
) -> (Vec<Plot>, Layout) {
    let colors = colors_for(kmers.len());
    // kmer -> color, in the order kmers were selected.
    let kmer_colors: Vec<(String, String)> =
        kmers.iter().cloned().zip(colors.iter().cloned()).collect();
    let color_to_kmer: HashMap<String, String> = kmer_colors
        .iter()
        .map(|(k, c)| (c.clone(), k.clone()))
        .collect();

    // Order rows according to the chosen sort mode.
    let mut ordered: Vec<&RepeatRecord> = records.iter().collect();
    match args.sort {
        SortMode::Alphabetic => ordered.sort_by(|a, b| b.sample.cmp(&a.sample)),
        SortMode::Length => ordered.sort_by_key(|r| r.sequence.len()),
        SortMode::Motif => {
            // Cluster by motif composition: sort by sequence length first, then
            // stably by the per-allele fraction of each selected motif (in
            // selection order, descending). Alleles dominated by the same motif
            // end up adjacent.
            ordered.sort_by_key(|r| r.sequence.len());
            let kmer_list: Vec<&str> = kmer_colors.iter().map(|(k, _)| k.as_str()).collect();
            let keys: HashMap<*const RepeatRecord, Vec<f64>> = ordered
                .iter()
                .map(|r| {
                    (
                        *r as *const RepeatRecord,
                        motif_key(&r.sequence, &kmer_list),
                    )
                })
                .collect();
            ordered.sort_by(|a, b| {
                let ka = &keys[&(*a as *const RepeatRecord)];
                let kb = &keys[&(*b as *const RepeatRecord)];
                for (x, y) in ka.iter().zip(kb.iter()) {
                    match y.partial_cmp(x).unwrap_or(std::cmp::Ordering::Equal) {
                        std::cmp::Ordering::Equal => continue,
                        o => return o,
                    }
                }
                std::cmp::Ordering::Equal
            });
        }
    }

    // One row per record (allele or somatic read), as a run-length-encoded
    // horizontal stacked bar: consecutive nucleotides with the same motif
    // collapse into a single segment, keeping the SVG small even for long
    // expansions. Each segment carries a tooltip with the sample, motif and
    // position span — so the interactive search box can locate a sample by name.
    let mut bar = BarPlot::new().with_horizontal(true).with_stacked();
    let mut tooltip_labels: Vec<String> = Vec::new();
    for rec in &ordered {
        let ident = rec.identifier(args.hide_allele_label);
        let label = if args.hide_labels {
            String::new()
        } else {
            ident.clone()
        };
        let mut pos = 0usize;
        let segments: Vec<(f64, String)> = run_length(&rec.sequence, &kmer_colors, &color_to_kmer)
            .into_iter()
            .map(|(len, color, motif)| {
                tooltip_labels.push(format!("{ident} · {motif} · {}-{}", pos, pos + len - 1));
                pos += len;
                (len as f64, color)
            })
            .collect();
        bar = bar.with_group(label, segments);
    }
    if interactive {
        bar = bar.with_tooltips().with_tooltip_labels(tooltip_labels);
    }
    let plots: Vec<Plot> = vec![Plot::Bar(bar)];

    // Manual legend: one swatch per selected motif (in selection order) + "other".
    let mut legend_entries: Vec<LegendEntry> = kmer_colors
        .iter()
        .map(|(k, c)| LegendEntry {
            label: k.clone(),
            color: c.clone(),
            shape: LegendShape::Rect,
            dasharray: None,
        })
        .collect();
    legend_entries.push(LegendEntry {
        label: "other".to_string(),
        color: GRAY.to_string(),
        shape: LegendShape::Rect,
        dasharray: None,
    });

    let title = if args.publication || args.minimal {
        args.title.clone()
    } else {
        format!("Repeat at {repeat}")
    };
    let y_label = if args.somatic {
        "Reads"
    } else {
        "Sample/allele"
    };

    // Height: explicit override, otherwise scale with the number of rows. When
    // labels are shown we need ~22 px/row to keep them legible; when hidden the
    // rows only need to clear the bar thickness, so they can be packed tighter.
    let n_rows = ordered.len();
    let height = args.height.map(|h| h as f64).unwrap_or_else(|| {
        let rows = n_rows as f64;
        let per_row = if args.hide_labels {
            (args.size as f64 * 2.0 + 4.0).max(8.0)
        } else {
            22.0
        };
        (rows * per_row + 150.0).max(800.0)
    });

    let mut layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_x_label("Nucleotide position in repeat")
        .with_y_label(y_label)
        .with_height(height)
        .with_tick_size(args.label_size)
        .with_legend_entries(legend_entries)
        // Tightly bound the y-axis to the rows; otherwise kuva rounds the range
        // outward to "nice" numbers, leaving a large empty band above the top row.
        .with_y_axis_min(0.5)
        .with_y_axis_max(n_rows as f64 + 0.5);

    if let Some(width) = args.width {
        layout = layout.with_width(width as f64);
    }

    // Legend placement.
    let legend_pos = match args.legend_corner {
        LegendCorner::TopRight => LegendPosition::InsideTopRight,
        LegendCorner::BottomRight => LegendPosition::InsideBottomRight,
        LegendCorner::Outside => LegendPosition::OutsideRightTop,
    };
    layout = layout.with_legend_position(legend_pos);

    // Publication / minimal styling.
    if args.publication || args.minimal {
        let mut theme = Theme::light();
        theme.background = "none".to_string();
        // Scale fonts, margins and legend swatches together (uniform enlargement)
        // rather than only the body font, which would leave the legend swatches
        // tiny next to oversized text.
        layout = layout.with_theme(theme).with_scale(1.5);
        if args.publication {
            layout = layout.with_box_axes();
        } else {
            // minimal: drop the y-axis title.
            layout = layout.with_y_label("");
        }
    }

    // Case annotations: a bold ">>>" to the left of the first nucleotide.
    for (i, rec) in ordered.iter().enumerate() {
        if rec.case {
            let y = (i + 1) as f64;
            layout = layout.with_annotation(
                TextAnnotation::new(">>>", -1.0, y)
                    .with_color("black")
                    .with_font_size(10),
            );
        }
    }

    if interactive {
        layout = layout.with_interactive();
    }

    (plots, layout)
}

/// Build the two-panel collapsed figure: a run-length barcode of the unique
/// representative sequences (left) beside a stacked case/control allele-count
/// histogram (right), sharing the y-axis. Rows are ordered by count, most
/// abundant at the top.
fn build_collapsed_scene(
    records: &[RepeatRecord],
    kmers: &[String],
    repeat: &str,
    args: &Args,
    interactive: bool,
    threshold: f64,
) -> Scene {
    let colors = colors_for(kmers.len());
    let kmer_colors: Vec<(String, String)> =
        kmers.iter().cloned().zip(colors.iter().cloned()).collect();
    let color_to_kmer: HashMap<String, String> = kmer_colors
        .iter()
        .map(|(k, c)| (c.clone(), k.clone()))
        .collect();

    let refs: Vec<&RepeatRecord> = records.iter().collect();
    let mut clusters = collapse(&refs, threshold);
    // Within equal counts, order clusters by the chosen --sort key (length,
    // representative sequence, or motif composition) instead of an arbitrary
    // string order, so same-count rows are grouped sensibly.
    let kmer_list: Vec<&str> = kmer_colors.iter().map(|(k, _)| k.as_str()).collect();
    match args.sort {
        SortMode::Length => clusters.sort_by_key(|c| c.representative.len()),
        SortMode::Alphabetic => clusters.sort_by(|a, b| a.representative.cmp(&b.representative)),
        SortMode::Motif => {
            clusters.sort_by(|a, b| motif_cmp(&a.representative, &b.representative, &kmer_list))
        }
    }
    // Stable sort by count (descending) keeps the secondary key within ties.
    clusters.sort_by_key(|c| std::cmp::Reverse(c.count));
    // `collapse` returns most-abundant-first; reverse so the highest count is
    // the top row (y categories run bottom-to-top).
    let display: Vec<&Cluster> = clusters.iter().rev().collect();
    let n = display.len();
    let has_cases = display.iter().any(|c| c.case_count > 0);

    // ── Barcode panel ────────────────────────────────────────────────────────
    let mut barcode = BarPlot::new().with_horizontal(true).with_stacked();
    let mut barcode_tips: Vec<String> = Vec::new();
    for c in &display {
        let label = if args.hide_labels {
            String::new()
        } else {
            format!("{}×", c.count)
        };
        let mut pos = 0usize;
        let segments: Vec<(f64, String)> =
            run_length(&c.representative, &kmer_colors, &color_to_kmer)
                .into_iter()
                .map(|(len, color, motif)| {
                    barcode_tips.push(format!(
                        "{} alleles · {} · {}-{}",
                        c.count,
                        motif,
                        pos,
                        pos + len - 1
                    ));
                    pos += len;
                    (len as f64, color)
                })
                .collect();
        barcode = barcode.with_group(label, segments);
    }
    if interactive {
        barcode = barcode.with_tooltips().with_tooltip_labels(barcode_tips);
    }

    // ── Count panel (stacked case / control) ─────────────────────────────────
    let mut countbar = BarPlot::new().with_horizontal(true).with_stacked();
    let mut count_tips: Vec<String> = Vec::new();
    for c in &display {
        countbar = countbar.with_group(
            String::new(),
            vec![
                (c.case_count as f64, CASE_COLOR.to_string()),
                (c.control_count as f64, CONTROL_COLOR.to_string()),
            ],
        );
        count_tips.push(format!("{} case alleles", c.case_count));
        count_tips.push(format!("{} control alleles", c.control_count));
    }
    if interactive {
        countbar = countbar.with_tooltips().with_tooltip_labels(count_tips);
    }

    // ── Legends ──────────────────────────────────────────────────────────────
    let mut motif_legend: Vec<LegendEntry> = kmer_colors
        .iter()
        .map(|(k, c)| LegendEntry {
            label: k.clone(),
            color: c.clone(),
            shape: LegendShape::Rect,
            dasharray: None,
        })
        .collect();
    motif_legend.push(LegendEntry {
        label: "other".to_string(),
        color: GRAY.to_string(),
        shape: LegendShape::Rect,
        dasharray: None,
    });
    let count_legend: Vec<LegendEntry> = if has_cases {
        vec![
            LegendEntry {
                label: "case".to_string(),
                color: CASE_COLOR.to_string(),
                shape: LegendShape::Rect,
                dasharray: None,
            },
            LegendEntry {
                label: "control".to_string(),
                color: CONTROL_COLOR.to_string(),
                shape: LegendShape::Rect,
                dasharray: None,
            },
        ]
    } else {
        vec![LegendEntry {
            label: "alleles".to_string(),
            color: CONTROL_COLOR.to_string(),
            shape: LegendShape::Rect,
            dasharray: None,
        }]
    };

    // ── Sizing ───────────────────────────────────────────────────────────────
    let height = args
        .height
        .map(|h| h as f64)
        .unwrap_or_else(|| (n as f64 * per_row_px(args) + 150.0).max(800.0));
    let barcode_w = args.width.map(|w| w as f64).unwrap_or(600.0);
    let count_w = 200.0;

    let title = if args.publication || args.minimal {
        args.title.clone()
    } else {
        format!("Repeat at {repeat}")
    };

    let barcode_plots = vec![Plot::Bar(barcode)];
    let count_plots = vec![Plot::Bar(countbar)];

    // The title goes on the Figure (above both panels), not on a single panel —
    // a panel-level title would consume top margin only in that cell and shift
    // its rows down relative to the other panel.
    let mut barcode_layout = Layout::auto_from_plots(&barcode_plots)
        .with_x_label("Nucleotide position in repeat")
        .with_y_label("Unique sequences (count ×)")
        .with_tick_size(args.label_size)
        .with_legend_entries(motif_legend)
        .with_legend_position(legend_position(args.legend_corner))
        .with_y_axis_min(0.5)
        .with_y_axis_max(n as f64 + 0.5);
    barcode_layout = apply_style(barcode_layout, args);

    let mut count_layout = Layout::auto_from_plots(&count_plots)
        .with_x_label("Allele count")
        .with_y_label("")
        .with_tick_size(args.label_size)
        .with_legend_entries(count_legend)
        .with_legend_position(LegendPosition::InsideTopRight)
        .with_y_axis_min(0.5)
        .with_y_axis_max(n as f64 + 0.5);
    count_layout = apply_style(count_layout, args);

    if interactive {
        barcode_layout = barcode_layout.with_interactive();
        count_layout = count_layout.with_interactive();
    }

    let mut scene = Figure::new(1, 2)
        .with_title(title)
        .with_plots(vec![barcode_plots, count_plots])
        .with_layouts(vec![barcode_layout, count_layout])
        .with_shared_y_all()
        .with_col_width(0, barcode_w)
        .with_col_width(1, count_w)
        .with_row_height(0, height)
        .render();
    // Figure::render does not propagate scene-level interactivity (only per-panel
    // tooltips), so the search/save UI strip is dropped. Re-enable it here so the
    // collapsed view keeps the search box (filter clusters by motif).
    if interactive {
        scene.interactive = true;
    }
    scene
}
