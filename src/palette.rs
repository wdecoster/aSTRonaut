//! Color palettes, ported from Plotly's qualitative `Safe` and `Light24`.

/// Gray used for nucleotides not belonging to any selected kmer ("other").
pub const GRAY: &str = "rgb(128, 128, 128)";

/// Plotly `qualitative.Safe` (Paul Tol's colour-blind-safe palette). Used when
/// there are at most 10 kmers.
pub const SAFE: [&str; 11] = [
    "rgb(136, 204, 238)",
    "rgb(204, 102, 119)",
    "rgb(221, 204, 119)",
    "rgb(17, 119, 51)",
    "rgb(51, 34, 136)",
    "rgb(170, 68, 153)",
    "rgb(68, 170, 153)",
    "rgb(153, 153, 51)",
    "rgb(136, 34, 85)",
    "rgb(102, 17, 0)",
    "rgb(136, 136, 136)",
];

/// Plotly `qualitative.Light24`. Used when there are more than 10 kmers.
pub const LIGHT24: [&str; 24] = [
    "#FD3216", "#00FE35", "#6A76FC", "#FED4C4", "#FE00CE", "#0DF9FF", "#F6F926", "#FF9616",
    "#479B55", "#EEA6FB", "#DC587D", "#D626FF", "#6E899C", "#00B5F7", "#B68E00", "#C9FBE5",
    "#FF0092", "#22FFA7", "#E3EE9E", "#86CE00", "#BC7196", "#7E7DCD", "#FC6955", "#E48F72",
];

/// Convert a `#rrggbb` hex string to an `rgb(r, g, b)` string. Ports
/// `hex_to_rgb`. This conversion matters: the sequence-coloring step replaces
/// leftover `A/C/T/G` characters, and hex tokens can contain `A`/`C`, so colors
/// must be in `rgb(...)` form (which contains no uppercase nucleotide letters).
pub fn hex_to_rgb(hex: &str) -> String {
    let h = hex.trim_start_matches('#');
    let r = u8::from_str_radix(&h[0..2], 16).unwrap_or(0);
    let g = u8::from_str_radix(&h[2..4], 16).unwrap_or(0);
    let b = u8::from_str_radix(&h[4..6], 16).unwrap_or(0);
    format!("rgb({}, {}, {})", r, g, b)
}

/// Return the colour list to use given the number of kmers. Mirrors the
/// `if len(kmers) > 10` branch in `plot_sequence`. All colours are returned as
/// `rgb(...)` strings.
pub fn colors_for(num_kmers: usize) -> Vec<String> {
    if num_kmers > 10 {
        if num_kmers > LIGHT24.len() {
            eprintln!("WARNING: Not enough colors defined!\n");
        }
        LIGHT24.iter().map(|s| hex_to_rgb(s)).collect()
    } else {
        SAFE.iter().map(|s| s.to_string()).collect()
    }
}
