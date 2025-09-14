//! Linting utilities for SMILES: collect diagnostics from lexing/parsing.

use lalrpop_util::ParseError;
use logos::Logos;

use super::lexer::{Lexer, Token};
use super::parser::grammar::MoleculeParser;
use crate::diagnostics::{Category, Code, Diagnostic, DiagnosticsReport, Span};
use crate::io::smiles::state::ParseState;

// Initial lexer-only lint to demonstrate diagnostics flow.
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
                Code("STYLE_SINGLE_DIGIT_RING"),
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
                // STYLE_BARE_ORGANIC: [C],[N],... with nothing else
                if is_bare_organic(inner) {
                    report.push(Diagnostic::warning(
                        Code("STYLE_BARE_ORGANIC"),
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
