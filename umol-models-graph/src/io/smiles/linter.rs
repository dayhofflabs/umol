//! Linting utilities for SMILES: collect diagnostics from lexing/parsing.

use lalrpop_util::ParseError;
use logos::Logos;

use super::lexer::{Lexer, Token};
use super::parser::grammar::MoleculeParser;
use crate::diagnostics::{Category, Code, Diagnostic, DiagnosticsReport, Span};
use crate::io::smiles::state::ParseState;
use umol_data::{isotope::Isotope as KnownIsotope, Element};

// Initial linter for lexical/syntactic errors.
pub fn lint_smiles(input: &str) -> DiagnosticsReport {
    let mut report = DiagnosticsReport::new();

    // Pass 1: lexical sweep
    let mut it = Lexer::new(input).peekable();
    while let Some(item) = it.next() {
        match item {
            Ok((start, tok, end)) => {
                // Percent token post-check is already enforced by lexer; nothing to do here for now
                let _ = (tok, start, end);
            }
            Err(_e) => {
                // Logos emits Token::Error via Lexer iterator; map to LEX_INVALID_TOKEN for safety
                // Fallback span: if we have a next ok token, use its start; otherwise 0..input.len()
                let span = Span::new(0, input.len());
                report.push(Diagnostic::error(
                    Code("LEX_INVALID_TOKEN"),
                    Category::Lex,
                    span,
                    "Invalid token",
                ));
            }
        }
    }

    // Secondary pass directly over raw logos to catch explicit Error tokens with spans
    for (res, span) in Token::lexer(input).spanned() {
        if res.is_err() {
            let slice = &input[span.start..span.end];
            if slice == "%" {
                report.push(Diagnostic::error(
                    Code("LEX_BAD_PERCENT_FORM"),
                    Category::Lex,
                    Span::new(span.start, span.end),
                    "'%' not followed by two digits",
                ));
            } else {
                report.push(Diagnostic::error(
                    Code("LEX_INVALID_TOKEN"),
                    Category::Lex,
                    Span::new(span.start, span.end),
                    "Invalid token",
                ));
            }
        }
    }

    // Quick lexical/syntactic heuristics
    lint_intertoken_whitespace(input, &mut report);
    lint_trailing_bond(input, &mut report);
    lint_dot_before_ring(input, &mut report);
    lint_dot_positions(input, &mut report);
    lint_style_bonds(input, &mut report);
    lint_ring_style(input, &mut report);

    // Style warnings from Standard Form (string-based heuristics)
    lint_style_percent_single_digit(input, &mut report);
    lint_brackets(input, &mut report);

    report
}

// Run the parser and translate lalrpop errors to diagnostics; also runs the lexical/style lint.
pub fn lint_smiles_parse(input: &str) -> DiagnosticsReport {
    let mut report = lint_smiles(input);

    // If lexical errors present, parsing is likely to cascade; still attempt a parse for location.
    let mut state = ParseState::default();
    let parser = MoleculeParser::new();
    let lexer = Lexer::new(input);
    let result = parser.parse(&mut state, lexer);
    if let Err(err) = result {
        match err {
            ParseError::InvalidToken { location } => {
                report.push(Diagnostic::error(
                    Code("SYN_UNEXPECTED_TOKEN"),
                    Category::Syn,
                    Span::new(location, location.saturating_add(1)),
                    "Unexpected token",
                ));
            }
            ParseError::UnrecognizedToken {
                token: (l, _tok, r),
                expected,
            } => {
                let mut d = Diagnostic::error(
                    Code("SYN_UNEXPECTED_TOKEN"),
                    Category::Syn,
                    Span::new(l, r),
                    "Unexpected token",
                );
                if !expected.is_empty() {
                    d = d.with_details(format!("expected one of: {}", expected.join(", ")));
                }
                report.push(d);
            }
            ParseError::UnrecognizedEof { location, expected } => {
                let mut d = Diagnostic::error(
                    Code("SYN_UNEXPECTED_TOKEN"),
                    Category::Syn,
                    Span::new(location, location),
                    "Unexpected end of input",
                );
                if !expected.is_empty() {
                    d = d.with_details(format!("expected one of: {}", expected.join(", ")));
                }
                report.push(d);
            }
            ParseError::ExtraToken {
                token: (l, _tok, r),
            } => {
                report.push(Diagnostic::error(
                    Code("SYN_UNEXPECTED_TOKEN"),
                    Category::Syn,
                    Span::new(l, r),
                    "Extra token",
                ));
            }
            ParseError::User { .. } => {
                report.push(Diagnostic::error(
                    Code("SYN_UNEXPECTED_TOKEN"),
                    Category::Syn,
                    Span::new(0, input.len()),
                    "Parse rejected",
                ));
            }
        }
    }
    // Merge any diagnostics emitted by the parser/state (e.g., ring rules)
    for d in state.diagnostics.into_iter() {
        report.push(d);
    }
    report
}

fn lint_style_percent_single_digit(input: &str, report: &mut DiagnosticsReport) {
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

fn lint_brackets(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == b'[' {
            if let Some(close) = find_closing_bracket(bytes, i + 1) {
                // Slice inside brackets
                let inner = &input[i + 1..close];
                // Parse inner fields (best-effort microparser)
                let parsed = parse_bracket_inner(inner);
                // Hydrogen element must not carry H-count
                if matches!(parsed.element, Some(Element::H)) && parsed.hcount.is_some() {
                    report.push(Diagnostic::error(
                        Code("BRKT_H_ON_H"),
                        Category::Brkt,
                        Span::new(i, close + 1),
                        "Hydrogen element must not have an H-count",
                    ));
                }
                // H-count two digits already handled; extra: H-count exceeds element's max implicit H
                if let (Some(elem), Some(h)) = (parsed.element, parsed.hcount) {
                    if h as u8 > elem.max_implicit_hydrogens() {
                        report.push(Diagnostic::warning(
                            Code("NUM_HCOUNT_EXCEEDS_MAX_IMPLICIT"),
                            Category::Num,
                            Span::new(i, close + 1),
                            "H-count exceeds element's max implicit hydrogens",
                        ));
                    }
                }
                // Class upper bound (max 4 digits)
                if let Some(class) = parsed.class {
                    if class > 9999 {
                        report.push(Diagnostic::error(
                            Code("NUM_CLASS_TOO_LARGE"),
                            Category::Num,
                            Span::new(i, close + 1),
                            "Atom class must be <= 9999",
                        ));
                    }
                }
                // Charge absolute limit
                if let Some(q) = parsed.charge {
                    if q.unsigned_abs() > 15 {
                        report.push(Diagnostic::error(
                            Code("NUM_CHARGE_OUT_OF_RANGE"),
                            Category::Num,
                            Span::new(i, close + 1),
                            "Absolute charge must be <= 15",
                        ));
                    }
                    if let Some(elem) = parsed.element {
                        let (min_q, max_q) = elem.charge_bounds();
                        if q < min_q as i32 || q > max_q as i32 {
                            report.push(Diagnostic::warning(
                                Code("NUM_CHARGE_OUTSIDE_ELEMENT_RANGE"),
                                Category::Num,
                                Span::new(i, close + 1),
                                "Charge outside element-supported bounds",
                            ));
                        }
                        if q > 0 && (q as u8) > elem.valence_electrons() {
                            report.push(Diagnostic::warning(
                                Code("NUM_CHARGE_EXCEEDS_VALENCE_ELECTRONS"),
                                Category::Num,
                                Span::new(i, close + 1),
                                "Positive charge exceeds valence electrons",
                            ));
                        }
                    }
                }
                // Isotope numeric limits and catalog
                if let Some(isotope) = parsed.isotope {
                    if isotope > 999 {
                        report.push(Diagnostic::error(
                            Code("NUM_ISOTOPE_TOO_LARGE"),
                            Category::Num,
                            Span::new(i, close + 1),
                            "Isotope mass number must be <= 999",
                        ));
                    } else if isotope > 0 {
                        if let Some(elem) = parsed.element {
                            if !KnownIsotope::is_catalogued(elem, isotope) {
                                report.push(Diagnostic::warning(
                                    Code("NUM_ISOTOPE_UNCATALOGUED"),
                                    Category::Num,
                                    Span::new(i, close + 1),
                                    "Isotope is not catalogued",
                                ));
                            }
                        }
                    }
                }
                // STYLE_BRACKET_ORGANIC: [C],[N],... with nothing else
                if is_bare_organic(inner) {
                    report.push(Diagnostic::warning(
                        Code("STYLE_BRACKET_ORGANIC"),
                        Category::Style,
                        Span::new(i, close + 1),
                        "Prefer bare organic atom over bracketed form",
                    ));
                }
                // STYLE_HCOUNT_ONE_SIMPLE: contains H1
                if inner_contains_h1(inner) {
                    // Span of 'H1' inside the brackets
                    if let Some((h_start, h_end)) = find_subslice(inner, "H1") {
                        report.push(Diagnostic::warning(
                            Code("STYLE_HCOUNT_ONE_SIMPLE"),
                            Category::Style,
                            Span::new(i + 1 + h_start, i + 1 + h_end),
                            "Prefer 'H' over 'H1'",
                        ));
                    } else {
                        report.push(Diagnostic::warning(
                            Code("STYLE_HCOUNT_ONE_SIMPLE"),
                            Category::Style,
                            Span::new(i, close + 1),
                            "Prefer 'H' over 'H1'",
                        ));
                    }
                }
                // STYLE_CHARGE_SIGN_SIMPLE: +1 or -1
                if let Some((c_start, c_end)) = find_charge_plus_minus_one(inner) {
                    report.push(Diagnostic::warning(
                        Code("STYLE_CHARGE_SIGN_SIMPLE"),
                        Category::Style,
                        Span::new(i + 1 + c_start, i + 1 + c_end),
                        "Prefer [+]/[-] over [+1]/[-1]",
                    ));
                }
                // BRKT_HCOUNT_TWO_DIGITS: H followed by two digits
                if let Some((h2s, h2e)) = find_h_two_digits(inner) {
                    report.push(Diagnostic::error(
                        Code("BRKT_HCOUNT_TWO_DIGITS"),
                        Category::Brkt,
                        Span::new(i + 1 + h2s, i + 1 + h2e),
                        "Hydrogen count must be a single digit",
                    ));
                }
                // BRKT_EMPTY_CLASS and NUM_CLASS_NEGATIVE
                if let Some((cs, ce, neg)) = find_class_issues(inner) {
                    if neg {
                        report.push(Diagnostic::error(
                            Code("NUM_CLASS_NEGATIVE"),
                            Category::Num,
                            Span::new(i + 1 + cs, i + 1 + ce),
                            "Atom class must be non-negative",
                        ));
                    } else {
                        report.push(Diagnostic::error(
                            Code("BRKT_EMPTY_CLASS"),
                            Category::Brkt,
                            Span::new(i + 1 + cs, i + 1 + ce),
                            "Class field ':' must be followed by digits",
                        ));
                    }
                }
                // STYLE_BRKT_ORDER: preferred order [chirality][H][charge][class]
                if bracket_order_misordered(inner) {
                    report.push(Diagnostic::warning(
                        Code("STYLE_BRKT_ORDER"),
                        Category::Style,
                        Span::new(i, close + 1),
                        "Prefer [chirality][H][charge][class] ordering",
                    ));
                }
                i = close + 1;
                continue;
            }
        }
        i += 1;
    }
}

fn find_closing_bracket(bytes: &[u8], mut idx: usize) -> Option<usize> {
    while idx < bytes.len() {
        if bytes[idx] == b']' {
            return Some(idx);
        }
        idx += 1;
    }
    None
}

fn is_bare_organic(inner: &str) -> bool {
    matches!(
        inner,
        "B" | "C" | "N" | "O" | "S" | "P" | "F" | "Cl" | "Br" | "I"
    )
}

fn inner_contains_h1(inner: &str) -> bool {
    // ensure 'H1' appears as a field token; simple substring check is sufficient in ASCII
    inner.contains("H1")
}

fn find_subslice(hay: &str, needle: &str) -> Option<(usize, usize)> {
    hay.find(needle).map(|s| (s, s + needle.len()))
}

fn find_charge_plus_minus_one(inner: &str) -> Option<(usize, usize)> {
    // Look for +1 or -1 not followed by another digit
    let bytes = inner.as_bytes();
    let mut i = 0usize;
    while i + 1 < bytes.len() {
        if (bytes[i] == b'+' || bytes[i] == b'-') && bytes[i + 1] == b'1' {
            let end = i + 2;
            if end >= bytes.len() || !bytes[end].is_ascii_digit() {
                return Some((i, end));
            }
        }
        i += 1;
    }
    None
}

fn bracket_order_misordered(inner: &str) -> bool {
    // Identify first occurrence indices for groups; -1 if absent
    let chiral_idx = find_first_chiral(inner);
    let h_idx = inner.find('H').map(|x| x as isize).unwrap_or(-1);
    let charge_idx = find_first_charge(inner);
    let class_idx = inner.find(':').map(|x| x as isize).unwrap_or(-1);
    // If fewer than 2 fields present (excluding element/isotope), do not warn
    let present = [chiral_idx, h_idx, charge_idx, class_idx]
        .into_iter()
        .filter(|&i| i >= 0)
        .count();
    if present < 2 {
        return false;
    }
    // Preferred non-decreasing order: chirality <= H <= charge <= class (when present)
    let mut last = -1isize;
    for idx in [chiral_idx, h_idx, charge_idx, class_idx] {
        if idx >= 0 {
            if idx < last {
                return true;
            }
            last = idx;
        }
    }
    false
}

fn find_first_chiral(inner: &str) -> isize {
    let patterns = ["@@", "@TH", "@AL", "@SP", "@TB", "@OH", "@"]; // order ensures '@@' before '@'
    let mut best: Option<usize> = None;
    for p in patterns.iter() {
        if let Some(i) = inner.find(p) {
            best = Some(best.map_or(i, |b| b.min(i)));
        }
    }
    best.map(|x| x as isize).unwrap_or(-1)
}

fn find_first_charge(inner: &str) -> isize {
    // '+' or '-' or '++'/'--' or '+d'/'-d'
    let mut best: Option<usize> = None;
    for p in ["++", "--", "+", "-"] {
        if let Some(i) = inner.find(p) {
            best = Some(best.map_or(i, |b| b.min(i)));
        }
    }
    best.map(|x| x as isize).unwrap_or(-1)
}

fn lint_trailing_bond(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut end = bytes.len();
    while end > 0 && matches!(bytes[end - 1], b' ' | b'\t' | b'\n' | b'\r') {
        end -= 1;
    }
    if end == 0 {
        return;
    }
    match bytes[end - 1] {
        b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\' => {
            report.push(Diagnostic::error(
                Code("SYN_TRAILING_BOND"),
                Category::Syn,
                Span::new(end - 1, end),
                "Trailing bond symbol",
            ));
        }
        _ => {}
    }
}

fn lint_dot_before_ring(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == b'.' {
            if i + 1 < bytes.len() {
                if bytes[i + 1].is_ascii_digit() {
                    report.push(Diagnostic::error(
                        Code("SYN_DOT_BEFORE_RING"),
                        Category::Syn,
                        Span::new(i, i + 2),
                        "Dot before ring index is invalid",
                    ));
                } else if bytes[i + 1] == b'%' {
                    // allow short span when incomplete; extend if two digits follow
                    let mut end = i + 2;
                    if i + 3 < bytes.len()
                        && bytes[i + 2].is_ascii_digit()
                        && bytes[i + 3].is_ascii_digit()
                    {
                        end = i + 4;
                    }
                    report.push(Diagnostic::error(
                        Code("SYN_DOT_BEFORE_RING"),
                        Category::Syn,
                        Span::new(i, end),
                        "Dot before ring index is invalid",
                    ));
                }
            }
        }
        i += 1;
    }
}

fn lint_intertoken_whitespace(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
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

fn find_h_two_digits(inner: &str) -> Option<(usize, usize)> {
    let bytes = inner.as_bytes();
    let mut i = 0usize;
    while i + 2 < bytes.len() {
        if bytes[i] == b'H' && bytes[i + 1].is_ascii_digit() && bytes[i + 2].is_ascii_digit() {
            return Some((i, i + 3));
        }
        i += 1;
    }
    None
}

fn find_class_issues(inner: &str) -> Option<(usize, usize, bool)> {
    let bytes = inner.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == b':' {
            if i + 1 >= bytes.len() {
                return Some((i, i + 1, false));
            }
            if bytes[i + 1] == b'-' {
                return Some((i, i + 2, true));
            }
            if !bytes[i + 1].is_ascii_digit() {
                return Some((i, i + 2, false));
            }
        }
        i += 1;
    }
    None
}

#[derive(Default, Debug, Clone, Copy)]
struct BracketParsed {
    element: Option<Element>,
    isotope: Option<u32>,
    hcount: Option<u32>,
    charge: Option<i32>,
    class: Option<u32>,
}

fn parse_bracket_inner(inner: &str) -> BracketParsed {
    let bytes = inner.as_bytes();
    let mut idx = 0usize;
    let mut parsed = BracketParsed::default();
    // Isotope (digits)
    while idx < bytes.len() && bytes[idx].is_ascii_digit() {
        idx += 1;
    }
    if idx > 0 {
        if let Some(v) = atoi::atoi::<u32>(&bytes[..idx]) { parsed.isotope = Some(v); }
    }
    // Symbol
    if idx < bytes.len() {
        if bytes[idx] == b'*' {
            idx += 1; // unknown element
        } else if bytes[idx].is_ascii_alphabetic() {
            let mut len = 1usize;
            if idx + 1 < bytes.len() && bytes[idx + 1].is_ascii_lowercase() {
                len = 2;
            }
            let sym = &inner[idx..idx + len];
            if let Some(el) = Element::from_symbol(sym) { parsed.element = Some(el); }
            idx += len;
        }
    }
    // Scan rest for fields
    let rest = &inner[idx..];
    // H-count: find standalone 'H' optionally followed by one digit
    if let Some(pos) = rest.find('H') {
        let hb = rest.as_bytes();
        let hpos = pos;
        let mut hval: u32 = 1;
        if hpos + 1 < hb.len() && hb[hpos + 1].is_ascii_digit() {
            hval = (hb[hpos + 1] - b'0') as u32;
            if hpos + 2 < hb.len() && hb[hpos + 2].is_ascii_digit() {
                // two digits handled elsewhere as error
            }
        }
        parsed.hcount = Some(hval);
    }
    // Class: ':' digits
    if let Some(cpos) = rest.find(':') {
        let s = &rest[cpos + 1..];
        let mut len = 0usize;
        for &b in s.as_bytes() { if b.is_ascii_digit() { len += 1; } else { break; } }
        if len > 0 { if let Some(v) = atoi::atoi::<u32>(s[..len].as_bytes()) { parsed.class = Some(v); } }
    }
    // Charge: handle '++','--','+','-','+d','-d','+dd','-dd'
    if let Some(pos) = rest.find('+') {
        let s = &rest[pos..];
        if s.starts_with("++") { parsed.charge = Some(2); }
        else {
            let mut val = 1i32;
            let digits = &s[1..];
            if let Some(d) = atoi::atoi::<u32>(digits.as_bytes()) { val = d as i32; }
            parsed.charge = Some(val);
        }
    } else if let Some(pos) = rest.find('-') {
        let s = &rest[pos..];
        if s.starts_with("--") { parsed.charge = Some(-2); }
        else {
            let mut val = -1i32;
            let digits = &s[1..];
            if let Some(d) = atoi::atoi::<u32>(digits.as_bytes()) { val = -(d as i32); }
            parsed.charge = Some(val);
        }
    }
    parsed
}

fn lint_style_bonds(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    for i in 0..bytes.len() {
        if bytes[i] == b':' {
            // Check neighbors roughly for aromatic atoms
            let prev = input[..i].as_bytes().iter().rfind(|&&b| b != b' ' && b != b'\t' && b != b'\n' && b != b'\r').copied();
            let next = input[i + 1..].as_bytes().iter().find(|&&b| b != b' ' && b != b'\t' && b != b'\n' && b != b'\r').copied();
            let is_arom = |b: u8| matches!(b, b'b'|b'c'|b'n'|b'o'|b'p'|b's');
            if prev.map(is_arom).unwrap_or(false) && next.map(is_arom).unwrap_or(false) {
                report.push(Diagnostic::warning(
                    Code("STYLE_EXPLICIT_AROMATIC_BOND"),
                    Category::Style,
                    Span::new(i, i + 1),
                    "Avoid explicit ':' between aromatic atoms",
                ));
            }
        } else if bytes[i] == b'-' {
            // Warn on explicit single bond unless between aromatic atoms
            let prev = input[..i].as_bytes().iter().rfind(|&&b| !b.is_ascii_whitespace()).copied();
            let next = input[i + 1..].as_bytes().iter().find(|&&b| !b.is_ascii_whitespace()).copied();
            let is_arom = |b: u8| matches!(b, b'b'|b'c'|b'n'|b'o'|b'p'|b's');
            if !(prev.map(is_arom).unwrap_or(false) && next.map(is_arom).unwrap_or(false)) {
                report.push(Diagnostic::warning(
                    Code("STYLE_EXPLICIT_SINGLE_BOND"),
                    Category::Style,
                    Span::new(i, i + 1),
                    "Avoid explicit '-' when default applies",
                ));
            }
        }
    }
}

fn lint_ring_style(input: &str, report: &mut DiagnosticsReport) {
    // Ignore digits inside brackets to avoid isotope/class confusion
    let mut i = 0usize;
    let bytes = input.as_bytes();
    let mut used: Vec<u32> = Vec::new();
    let mut counts: indexmap::IndexMap<u32, u32> = indexmap::IndexMap::new();
    while i < bytes.len() {
        match bytes[i] {
            b'[' => {
                // skip to closing ']'
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
        if first != 1 { report.push(Diagnostic::warning(Code("STYLE_FIRST_RING_NOT_ONE"), Category::Style, Span::new(0, 0), "Prefer starting ring numbering at 1")); }
        if set.len() >= 2 {
            let mut prev = set[0];
            for &v in &set[1..] {
                if v > prev + 1 { report.push(Diagnostic::warning(Code("STYLE_NONCONSECUTIVE_RING_NUMBERING"), Category::Style, Span::new(0, 0), "Prefer consecutive ring numbering")); break; }
                prev = v;
            }
        }
        for (_k, c) in counts.iter() {
            if *c > 2 { report.push(Diagnostic::warning(Code("STYLE_REUSED_RING_INDICES"), Category::Style, Span::new(0, 0), "Avoid reusing the same ring number")); break; }
        }
    }
}

fn lint_dot_positions(input: &str, report: &mut DiagnosticsReport) {
    let bytes = input.as_bytes();
    // leading dot
    let mut i = 0usize;
    while i < bytes.len() && bytes[i].is_ascii_whitespace() { i += 1; }
    if i < bytes.len() && bytes[i] == b'.' {
        report.push(Diagnostic::error(Code("SYN_LEADING_DOT"), Category::Syn, Span::new(i, i + 1), "Leading dot"));
    }
    // trailing dot
    let mut j = bytes.len();
    while j > 0 && bytes[j - 1].is_ascii_whitespace() { j -= 1; }
    if j > 0 && bytes[j - 1] == b'.' {
        report.push(Diagnostic::error(Code("SYN_TRAILING_DOT"), Category::Syn, Span::new(j - 1, j), "Trailing dot"));
    }
    // multiple dots '..'
    let mut k = 0usize;
    while k + 1 < bytes.len() {
        if bytes[k] == b'.' && bytes[k + 1] == b'.' {
            report.push(Diagnostic::error(Code("SYN_MULTIPLE_DOTS"), Category::Syn, Span::new(k, k + 2), "Multiple dots"));
            break;
        }
        k += 1;
    }
}
