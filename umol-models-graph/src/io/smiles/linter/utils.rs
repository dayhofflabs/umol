//! Shared SMILES linter utilities: bracket helpers and style helpers.

use bstr::ByteSlice;
use indexmap::IndexMap;
use memchr::{memchr, memchr2_iter, memchr_iter};
use regex::Regex;

use crate::diagnostics::{Category, Code, Diagnostic, DiagnosticsReport, Span};
use crate::io::smiles::parser::utils::{
    parse_bracket as parse_bracket_parser, BracketField, BracketFields,
};

/// Returns the index of the first ']' at or after `idx`.
pub fn find_closing_bracket(bytes: &[u8], idx: usize) -> Option<usize> {
    bytes
        .get(idx..)
        .and_then(|r| memchr(b']', r))
        .map(|o| idx + o)
}

pub fn is_bare_organic(s: &str) -> bool {
    matches!(
        s,
        "B" | "C" | "N" | "O" | "S" | "P" | "F" | "Cl" | "Br" | "I"
    )
}

pub fn inner_contains_h1(s: &str) -> bool {
    s.as_bytes().contains_str(b"H1")
}

pub fn find_subslice(hay: &str, needle: &str) -> Option<(usize, usize)> {
    hay.find(needle).map(|s| (s, s + needle.len()))
}

pub fn find_charge_plus_minus_one(s: &str) -> Option<(usize, usize)> {
    let bytes = s.as_bytes();
    for idx in memchr2_iter(b'+', b'-', bytes) {
        if idx + 1 < bytes.len() && bytes[idx + 1] == b'1' {
            let end = idx + 2;
            if end >= bytes.len() || !bytes[end].is_ascii_digit() {
                return Some((idx, end));
            }
        }
    }
    None
}

/// Returns true if bracket tail fields are out of canonical order.
/// Canonical order: Chirality -> H count -> Charge -> Class.
/// Uses the shared bracket parser once and checks the order of parsed fields.
pub fn bracket_order_misordered(s: &str) -> bool {
    let (_elem, _iso, tails) = parse_bracket_parser(s);
    let mut prev: Option<usize> = None;
    for f in tails.iter() {
        let id = match f {
            BracketField::Chiral(_) => 0,
            BracketField::HydrogenCount(_) => 1,
            BracketField::Charge(_) => 2,
            BracketField::Class(_) => 3,
        };
        if let Some(p) = prev {
            if id < p {
                return true;
            }
        }
        prev = Some(id);
    }
    false
}

pub fn find_h_two_digits(s: &str) -> Option<(usize, usize)> {
    let re = Regex::new(r"H\d{2,}").unwrap();
    re.find(s).map(|m| (m.start(), m.end()))
}

/// Detects invalid class index after ':' in a bracket tail.
/// Returns (start, end, is_negative) where:
/// - is_negative = true if immediately followed by '-'
/// - otherwise true if non-digit or end-of-input after ':'
pub fn invalid_class_index(s: &str) -> Option<(usize, usize, bool)> {
    let bytes = s.as_bytes();
    for idx in memchr_iter(b':', bytes) {
        let next = idx + 1;
        if next >= bytes.len() {
            return Some((idx, idx + 1, false));
        }
        match bytes[next] {
            b'-' => return Some((idx, next + 1, true)),
            d if !d.is_ascii_digit() => return Some((idx, next + 1, false)),
            _ => {}
        }
    }
    None
}

pub fn lint_trailing_bond(s: &str, report: &mut DiagnosticsReport) {
    let trimmed = s.as_bytes().trim_end_with(|c| c.is_ascii_whitespace());
    if trimmed.is_empty() {
        return;
    }
    if matches!(
        trimmed[trimmed.len() - 1],
        b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\'
    ) {
        report.push(Diagnostic::error(
            Code("SYN_TRAILING_BOND"),
            Category::Syn,
            Span::new(trimmed.len() - 1, trimmed.len()),
            "Trailing bond symbol",
        ));
    }
}

pub fn lint_dot_before_ring(s: &str, report: &mut DiagnosticsReport) {
    // Find ".<digit>" or ".%<two digits>", but not inside brackets
    let re = Regex::new(r"\.(?:[0-9]|%[0-9]{2})").unwrap();
    let mut i = 0usize;
    let bytes = s.as_bytes();
    while i < bytes.len() {
        if bytes[i] == b'[' {
            let mut j = i + 1;
            while j < bytes.len() && bytes[j] != b']' {
                j += 1;
            }
            i = j.saturating_add(1);
            continue;
        }
        if let Some(m) = re.find(&s[i..]) {
            let ss = i + m.start();
            let ee = i + m.end();
            report.push(Diagnostic::error(
                Code("SYN_DOT_BEFORE_RING"),
                Category::Syn,
                Span::new(ss, ee),
                "Dot before ring index is invalid",
            ));
            i = ee;
            continue;
        }
        i += 1;
    }
}

pub fn lint_intertoken_whitespace(s: &str, report: &mut DiagnosticsReport) {
    let bytes = s.as_bytes();
    let mut last_non_ws: isize = -1;
    for (i, &b) in bytes.iter().enumerate() {
        if !matches!(b, b' ' | b'\t' | b'\n' | b'\r') {
            last_non_ws = i as isize;
        }
    }
    if last_non_ws < 0 {
        return;
    }
    let mut i = 0usize;
    while i < (last_non_ws as usize) {
        if matches!(bytes[i], b' ' | b'\t' | b'\n' | b'\r') {
            let start = i;
            while i < bytes.len() && matches!(bytes[i], b' ' | b'\t' | b'\n' | b'\r') {
                i += 1;
            }
            report.push(Diagnostic::error(
                Code("LEX_INTERTOKEN_WHITESPACE"),
                Category::Lex,
                Span::new(start, i),
                "Inter-token whitespace is not allowed",
            ));
            continue;
        }
        i += 1;
    }
}

pub fn lint_dot_positions(s: &str, report: &mut DiagnosticsReport) {
    // Leading dot
    let re_lead = Regex::new(r"^\s*\.").unwrap();
    if let Some(m) = re_lead.find(s) {
        let start = m.start();
        report.push(Diagnostic::error(
            Code("SYN_LEADING_DOT"),
            Category::Syn,
            Span::new(start + (m.as_str().len() - 1), start + m.as_str().len()),
            "Leading dot",
        ));
    }
    // Trailing dot
    let re_trail = Regex::new(r"\.\s*$").unwrap();
    if let Some(m) = re_trail.find(s) {
        report.push(Diagnostic::error(
            Code("SYN_TRAILING_DOT"),
            Category::Syn,
            Span::new(m.start(), m.start() + 1),
            "Trailing dot",
        ));
    }
    // Multiple consecutive dots
    let re_multi = Regex::new(r"\.\.").unwrap();
    if let Some(m) = re_multi.find(s) {
        report.push(Diagnostic::error(
            Code("SYN_MULTIPLE_DOTS"),
            Category::Syn,
            Span::new(m.start(), m.end()),
            "Multiple dots",
        ));
    }
}

pub fn lint_style_percent_single_digit(input: &str, report: &mut DiagnosticsReport) {
    let re = Regex::new(r"%(0[1-9])").unwrap();
    for m in re.find_iter(input) {
        let start = m.start();
        report.push(Diagnostic::warning(
            Code("STYLE_UNNECESSARY_PERCENT_RING_INDEX"),
            Category::Style,
            Span::new(start, start + 3),
            "Prefer single-digit ring number for 1..9",
        ));
    }
}

pub fn lint_style_bonds(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        match bytes[i] {
            b'[' => {
                if let Some(off) = memchr(b']', &bytes[i + 1..]) {
                    i = i + 1 + off + 1;
                } else {
                    break;
                }
                continue;
            }
            b':' => {
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
                while j < bytes.len() && bytes[j] != b']' {
                    j += 1;
                }
                i = j + 1;
                continue;
            }
            b'%' => {
                if i + 2 < bytes.len()
                    && bytes[i + 1].is_ascii_digit()
                    && bytes[i + 2].is_ascii_digit()
                {
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

pub fn parse_bracket(inner: &str) -> BracketFields {
    let (elem, iso, tails) = parse_bracket_parser(inner);
    let mut parsed = BracketFields::default();
    parsed.element = elem;
    parsed.isotope = iso;
    for f in tails {
        match f {
            BracketField::Chiral(_) => {}
            BracketField::HydrogenCount(v) => parsed.hcount = Some(v),
            BracketField::Charge(q) => parsed.charge = Some(q),
            BracketField::Class(c) => parsed.class = Some(c),
        }
    }
    parsed
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(b"[CH3]".as_slice(), 1, Some(4))]
    #[case(b"[CH3".as_slice(), 1, None)]
    fn test_find_closing_bracket(
        #[case] s: &[u8],
        #[case] start: usize,
        #[case] expect: Option<usize>,
    ) {
        assert_eq!(find_closing_bracket(s, start), expect);
    }

    #[rstest]
    #[case("+1", Some((0, 2)))]
    #[case("-1]", Some((0, 2)))]
    #[case("+10", None)]
    #[case("+2", None)]
    fn test_find_charge_plus_minus_one(#[case] s: &str, #[case] expect: Option<(usize, usize)>) {
        assert_eq!(find_charge_plus_minus_one(s), expect);
    }

    #[rstest]
    #[case(":", Some((0, 1, false)))]
    #[case(":x", Some((0, 2, false)))]
    #[case(":-1", Some((0, 2, true)))]
    #[case(":12", None)]
    fn test_invalid_class_index(#[case] s: &str, #[case] expect: Option<(usize, usize, bool)>) {
        assert_eq!(invalid_class_index(s), expect);
    }

    #[rstest]
    // Misordered examples
    #[case("*H@", true)]
    #[case("CH@", true)]
    #[case("*+1H", true)]
    #[case("C:1+1", true)]
    // Ordered examples
    #[case("*@H", false)]
    #[case("C@H", false)]
    #[case("C@H+1:1", false)]
    #[case("CH+1:2", false)]
    fn test_bracket_order_misordered(#[case] inner: &str, #[case] misordered: bool) {
        assert_eq!(bracket_order_misordered(inner), misordered);
    }

    #[rstest]
    #[case("C-  ", true)]
    #[case("C-", true)]
    #[case("C :", true)]
    #[case("C ", false)]
    #[case("C", false)]
    fn test_lint_trailing_bond(#[case] s: &str, #[case] is_err: bool) {
        let mut rep = DiagnosticsReport::new();
        lint_trailing_bond(s, &mut rep);
        assert_eq!(
            rep.diagnostics
                .iter()
                .any(|d| d.code.0 == "SYN_TRAILING_BOND"),
            is_err
        );
    }
}
