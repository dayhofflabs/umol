//! Style-related lint helpers used by rules.

use crate::diagnostics::{Category, Code, Diagnostic, DiagnosticsReport, Span};
use indexmap::IndexMap;

pub fn lint_style_percent_single_digit(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut i = 0usize;
    while i + 2 < bytes.len() {
        if bytes[i] == b'%' && bytes[i + 1] == b'0' && (b'1'..=b'9').contains(&bytes[i + 2]) {
            report.push(Diagnostic::warning(
                Code("STYLE_UNNECESSARY_PERCENT_RING_INDEX"),
                Category::Style,
                Span::new(i, i + 3),
                "Prefer single-digit ring number for 1..9",
            ));
            i += 3;
            continue;
        }
        i += 1;
    }
}

pub fn lint_style_bonds(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        match bytes[i] {
            b'[' => {
                let mut j = i + 1;
                while j < bytes.len() && bytes[j] != b']' { j += 1; }
                i = j.saturating_add(1);
                continue;
            }
            b':' => {
                let prev = input[..i]
                    .as_bytes()
                    .iter()
                    .rfind(|&&b| b != b' ' && b != b'\t' && b != b'\n' && b != b'\r')
                    .copied();
                let next = input[i + 1..]
                    .as_bytes()
                    .iter()
                    .find(|&&b| b != b' ' && b != b'\t' && b != b'\n' && b != b'\r')
                    .copied();
                let is_arom = |b: u8| matches!(b, b'b' | b'c' | b'n' | b'o' | b'p' | b's');
                if prev.map(is_arom).unwrap_or(false) && next.map(is_arom).unwrap_or(false) {
                    report.push(Diagnostic::warning(
                        Code("STYLE_EXPLICIT_AROMATIC_BOND"),
                        Category::Style,
                        Span::new(i, i + 1),
                        "Avoid explicit ':' between aromatic atoms",
                    ));
                }
            }
            b'-' => {
                let prev = input[..i]
                    .as_bytes()
                    .iter()
                    .rfind(|&&b| !b.is_ascii_whitespace())
                    .copied();
                let next = input[i + 1..]
                    .as_bytes()
                    .iter()
                    .find(|&&b| !b.is_ascii_whitespace())
                    .copied();
                let is_arom = |b: u8| matches!(b, b'b' | b'c' | b'n' | b'o' | b'p' | b's');
                if !(prev.map(is_arom).unwrap_or(false) && next.map(is_arom).unwrap_or(false)) {
                    report.push(Diagnostic::warning(
                        Code("STYLE_EXPLICIT_SINGLE_BOND"),
                        Category::Style,
                        Span::new(i, i + 1),
                        "Avoid explicit '-' when default applies",
                    ));
                }
            }
            _ => {}
        }
        i += 1;
    }
}

pub fn lint_ring_style(input: &str, report: &mut DiagnosticsReport) {
    let mut i = 0usize;
    let bytes = input.as_bytes();
    let mut used: Vec<u32> = Vec::new();
    let mut counts: IndexMap<u32, u32> = IndexMap::new();
    while i < bytes.len() {
        match bytes[i] {
            b'[' => {
                let mut j = i + 1;
                while j < bytes.len() && bytes[j] != b']' { j += 1; }
                i = j + 1;
                continue;
            }
            b'%' => {
                if i + 2 < bytes.len() && bytes[i + 1].is_ascii_digit() && bytes[i + 2].is_ascii_digit() {
                    let val = (bytes[i + 1] - b'0') as u32 * 10 + (bytes[i + 2] - b'0') as u32;
                    used.push(val);
                    *counts.entry(val).or_insert(0) += 1;
                    i += 3;
                    continue;
                }
            }
            b'0'..=b'9' => {
                let val = (bytes[i] - b'0') as u32;
                used.push(val);
                *counts.entry(val).or_insert(0) += 1;
                i += 1;
                continue;
            }
            _ => {}
        }
        i += 1;
    }
    if !used.is_empty() {
        let mut set = used.clone();
        set.sort_unstable();
        set.dedup();
        let first = set[0];
        if first != 1 {
            report.push(Diagnostic::warning(
                Code("STYLE_FIRST_RING_NOT_ONE"),
                Category::Style,
                Span::new(0, 0),
                "Prefer starting ring numbering at 1",
            ));
        }
        if set.len() >= 2 {
            let mut prev = set[0];
            for &v in &set[1..] {
                if v > prev + 1 {
                    report.push(Diagnostic::warning(
                        Code("STYLE_NONCONSECUTIVE_RING_NUMBERING"),
                        Category::Style,
                        Span::new(0, 0),
                        "Prefer consecutive ring numbering",
                    ));
                    break;
                }
                prev = v;
            }
        }
        for (_k, c) in counts.iter() {
            if *c > 2 {
                report.push(Diagnostic::warning(
                    Code("STYLE_REUSED_RING_INDICES"),
                    Category::Style,
                    Span::new(0, 0),
                    "Avoid reusing the same ring number",
                ));
                break;
            }
        }
    }
}

pub fn lint_dot_positions(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() && bytes[i].is_ascii_whitespace() { i += 1; }
    if i < bytes.len() && bytes[i] == b'.' {
        report.push(Diagnostic::error(
            Code("SYN_LEADING_DOT"),
            Category::Syn,
            Span::new(i, i + 1),
            "Leading dot",
        ));
    }
    let mut j = bytes.len();
    while j > 0 && bytes[j - 1].is_ascii_whitespace() { j -= 1; }
    if j > 0 && bytes[j - 1] == b'.' {
        report.push(Diagnostic::error(
            Code("SYN_TRAILING_DOT"),
            Category::Syn,
            Span::new(j - 1, j),
            "Trailing dot",
        ));
    }
    let mut k = 0usize;
    while k + 1 < bytes.len() {
        if bytes[k] == b'.' && bytes[k + 1] == b'.' {
            report.push(Diagnostic::error(
                Code("SYN_MULTIPLE_DOTS"),
                Category::Syn,
                Span::new(k, k + 2),
                "Multiple dots",
            ));
            break;
        }
        k += 1;
    }
}


