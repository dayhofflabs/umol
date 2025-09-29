//! Linting for SMILES: collect diagnostics from lexing/parsing.

mod context;
mod emitter;
mod registry;
mod rules;
pub use context::LintContext;
pub use emitter::{DiagnosticCandidate, Emitter, Scope};
pub use registry::{LintEngine, RuleRegistry};
pub use rules::Rule;

use super::parser::parse_smiles;
use crate::diagnostics::{Category, Code, DiagnosticsReport, Severity, Span};
use crate::io::smiles::ParseError;
use crate::io::smiles::parser::utils::{is_valid_bracket_inner, parse_bracket, BracketField};


// SMILES linter, runs post-parse
pub fn lint_smiles(input: &str) -> DiagnosticsReport {
    let ctx = LintContext::new(input);
    let mut report = DiagnosticsReport::new();
    let mut emitter = Emitter::new(&mut report);

    // Map parser errors into diagnostics
    let parse_res = parse_smiles(input.as_bytes());
    if let Err(err) = &parse_res {
        match *err {
            // Ring diagnostics
            ParseError::RingUnclosed { open_pos } => {
                emitter.candidate(DiagnosticCandidate {
                    code: Code("RING_UNCLOSED"),
                    category: Category::Ring,
                    severity: Severity::Error,
                    span: Span::new(open_pos, open_pos),
                    message: "Ring index opened but not closed",
                    scope: Scope::Global,
                });
            }
            ParseError::RingSelfLoop { pos } => {
                emitter.candidate(DiagnosticCandidate {
                    code: Code("RING_SELF_LOOP"),
                    category: Category::Ring,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Ring closure creates a self-loop",
                    scope: Scope::Global,
                });
            }
            ParseError::RingTwoMember { pos } => {
                emitter.candidate(DiagnosticCandidate {
                    code: Code("RING_TWO_MEMBER"),
                    category: Category::Ring,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Ring closure creates a two-member ring",
                    scope: Scope::Global,
                });
            }
            ParseError::RingBondDirConflict { pos, .. } => {
                emitter.candidate(DiagnosticCandidate {
                    code: Code("RING_BOND_DIR_CONFLICT"),
                    category: Category::Ring,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Conflicting ring bond directions",
                    scope: Scope::Global,
                });
            }
            ParseError::RingBondOrderConflict { pos, .. } => {
                emitter.candidate(DiagnosticCandidate {
                    code: Code("RING_BOND_ORDER_CONFLICT"),
                    category: Category::Ring,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Conflicting ring bond orders",
                    scope: Scope::Global,
                });
            }

            // Bracket balance
            ParseError::UnbalancedOpenBracket { pos } => {
                emitter.candidate(DiagnosticCandidate {
                    code: Code("PARSER_UNBALANCED_OPEN_BRACKET"),
                    category: Category::Bracket,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Unbalanced '[' bracket",
                    scope: Scope::Bracket { start: pos, end: pos.saturating_add(1) },
                });
            }
            ParseError::UnbalancedCloseBracket { pos } => {
                emitter.candidate(DiagnosticCandidate {
                    code: Code("PARSER_UNBALANCED_CLOSE_BRACKET"),
                    category: Category::Bracket,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Unbalanced ']' bracket",
                    scope: Scope::Bracket { start: pos, end: pos.saturating_add(1) },
                });
            }

            // Branch/group
            ParseError::UnbalancedBranchOpen { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("BRCH_UNCLOSED"), category: Category::Branch, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Open '(' not closed", scope: Scope::Global });
            }
            ParseError::UnbalancedBranchClose { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("BRCH_UNEXPECTED_CLOSE"), category: Category::Branch, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Unexpected ')'", scope: Scope::Global });
            }
            ParseError::EmptyBranch { pos } | ParseError::EmptyGroup { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("BRCH_EMPTY_BRANCH"), category: Category::Branch, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Empty branch/group", scope: Scope::Global });
            }

            // Bond/dot
            ParseError::LeadingBond { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("LEX_LEADING_BOND"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Leading bond token", scope: Scope::Global });
            }
            ParseError::ConsecutiveBond { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("LEX_CONSECUTIVE_BONDS"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Consecutive bond tokens", scope: Scope::Global });
            }
            ParseError::TrailingBond { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("LEX_TRAILING_BOND"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Trailing bond token", scope: Scope::Global });
                // Derive dangling bond before component end or ')'
                if next_non_ws_is(input.as_bytes(), pos, b')') || next_non_ws_is(input.as_bytes(), pos, b'.') || pos + 1 >= input.len() {
                    emitter.candidate(DiagnosticCandidate { code: Code("BRCH_DANGLING_BOND"), category: Category::Branch, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Dangling bond before group end or component boundary", scope: Scope::Global });
                }
            }
            ParseError::LeadingDot { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("LEX_LEADING_DOT"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Leading dot", scope: Scope::Global });
            }
            ParseError::ConsecutiveDot { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("LEX_MULTIPLE_DOTS"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Consecutive dots", scope: Scope::Global });
            }
            ParseError::TrailingDot { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("LEX_TRAILING_DOT"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Trailing dot", scope: Scope::Global });
            }

            // Ring leading and index invalid
            ParseError::LeadingRing { pos } => {
                if pos > 0 && input.as_bytes()[pos - 1] == b'.' {
                    emitter.candidate(DiagnosticCandidate { code: Code("LEX_DOT_BEFORE_RING"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos - 1, pos + 1), message: "Dot before ring index", scope: Scope::Global });
                }
                emitter.candidate(DiagnosticCandidate { code: Code("LEX_LEADING_RING"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Leading ring index", scope: Scope::Global });
            }
            ParseError::RingIndexInvalid { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("LEX_RING_INDEX_INVALID"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Invalid percent ring index", scope: Scope::Global });
            }

            // Whitespace/comments (strict mode)
            ParseError::InvalidWhitespace { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("LEX_INTERTOKEN_WHITESPACE"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Inter-token whitespace is not allowed", scope: Scope::Global });
            }
            ParseError::InvalidComment { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("LEX_COMMENT"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(2)), message: "Comments are not allowed", scope: Scope::Global });
            }
            ParseError::UnterminatedBlockComment { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("LEX_UNTERMINATED_BLOCK_COMMENT"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, input.len()), message: "Unterminated block comment", scope: Scope::Global });
            }

            _ => {}
        }
    }

    // Style and numeric advisories
    run_style_and_numeric_checks(&ctx, &mut emitter, parse_res.is_ok());

    emitter.flush();
    report
}

fn next_non_ws_is(bytes: &[u8], pos: usize, ch: u8) -> bool {
    let mut i = pos + 1;
    while i < bytes.len() && matches!(bytes[i], b' ' | b'\t' | b'\n' | b'\r') { i += 1; }
    i < bytes.len() && bytes[i] == ch
}

fn run_style_and_numeric_checks(ctx: &LintContext, emit: &mut Emitter, only_when_parse_ok: bool) {
    let input = ctx.input;
    let bytes = input.as_bytes();

    // Percent ring index style: %01..%09
    let mut i = 0usize;
    while i + 2 < bytes.len() {
        if bytes[i] == b'%' && bytes[i + 1] == b'0' && (bytes[i + 2] >= b'1' && bytes[i + 2] <= b'9') {
            emit.candidate(DiagnosticCandidate { code: Code("STYLE_UNNECESSARY_PERCENT_RING_INDEX"), category: Category::Style, severity: Severity::Warning, span: Span::new(i, i + 3), message: "Prefer single-digit ring index for 1..9", scope: Scope::Global });
        }
        i += 1;
    }

    // Bracket-based style/num checks
    let mut j = 0usize;
    while j < bytes.len() {
        if bytes[j] == b'[' {
            let start = j;
            let mut k = j + 1;
            while k < bytes.len() && bytes[k] != b']' { k += 1; }
            if k >= bytes.len() { break; }
            if let Ok(inner) = std::str::from_utf8(&bytes[start + 1..k]) {
                // Only check style against valid bracket substrings
                if is_valid_bracket_inner(inner) {
                    let (elem_opt, iso_opt, fields) = parse_bracket(inner);

                    // STYLE_BRKT_ORGANIC
                    if fields.is_empty() && iso_opt.is_none() {
                        if let Some(elem) = elem_opt {
                            if is_organic_subset(elem) {
                                emit.candidate(DiagnosticCandidate { code: Code("STYLE_BRKT_ORGANIC"), category: Category::Style, severity: Severity::Warning, span: Span::new(start, k + 1), message: "Prefer bare organic atom over bracketed form", scope: Scope::Bracket { start, end: k + 1 } });
                            }
                        }
                    }

                    // STYLE_BRKT_ORDER
                    let mut last_rank = 0u8;
                    let mut ordered = true;
                    for f in &fields {
                        let r = field_rank(f);
                        if r < last_rank { ordered = false; break; }
                        last_rank = r;
                    }
                    if !ordered {
                        emit.candidate(DiagnosticCandidate { code: Code("STYLE_BRKT_ORDER"), category: Category::Style, severity: Severity::Warning, span: Span::new(start, k + 1), message: "Prefer [chirality][H][charge][class] ordering", scope: Scope::Bracket { start, end: k + 1 } });
                    }

                    // NUM_ISOTOPE_TOO_LARGE
                    if let Some(iso) = iso_opt { if iso >= 1000 { emit.candidate(DiagnosticCandidate { code: Code("NUM_ISOTOPE_TOO_LARGE"), category: Category::Num, severity: Severity::Error, span: Span::new(start + 1, k), message: "Isotope mass number too large", scope: Scope::Bracket { start, end: k + 1 } }); } }
                }

                // STYLE_CHARGE_SIGN_SIMPLE: detect +1 / -1 literally
                if contains_sign_one(inner) {
                    emit.candidate(DiagnosticCandidate { code: Code("STYLE_CHARGE_SIGN_SIMPLE"), category: Category::Style, severity: Severity::Warning, span: Span::new(start, k + 1), message: "Prefer [+]/[-] over [+1]/[-1]", scope: Scope::Bracket { start, end: k + 1 } });
                }
                // STYLE_HCOUNT_ONE_SIMPLE: detect H1 literally
                if contains_h1(inner) {
                    emit.candidate(DiagnosticCandidate { code: Code("STYLE_HCOUNT_ONE_SIMPLE"), category: Category::Style, severity: Severity::Warning, span: Span::new(start, k + 1), message: "Prefer H over H1", scope: Scope::Bracket { start, end: k + 1 } });
                }

                // NUM_CHIRAL_OUT_OF_RANGE: detect @TBn/@OHn over limits by raw pattern
                if chiral_param_out_of_range(inner) {
                    emit.candidate(DiagnosticCandidate { code: Code("NUM_CHIRAL_OUT_OF_RANGE"), category: Category::Num, severity: Severity::Error, span: Span::new(start + 1, k), message: "Chirality parameter out of range", scope: Scope::Bracket { start, end: k + 1 } });
                }
            }
            j = k + 1;
            continue;
        }
        j += 1;
    }

    // Ring numbering style checks only when parsing succeeded to avoid noise
    if only_when_parse_ok {
        let seq = ring_indices_sequence(bytes);
        if let Some((first_num, s, e)) = seq.first().copied() {
            if first_num != 1 {
                emit.candidate(DiagnosticCandidate { code: Code("STYLE_FIRST_RING_NOT_ONE"), category: Category::Style, severity: Severity::Warning, span: Span::new(s, e), message: "Prefer starting ring numbering at 1", scope: Scope::Global });
            }
        }
        // Non-consecutive jumps
        let mut last: Option<u32> = None;
        for (num, s, e) in seq.into_iter() {
            if let Some(p) = last {
                if num > p + 1 { emit.candidate(DiagnosticCandidate { code: Code("STYLE_NONCONSECUTIVE_RING_NUMBERING"), category: Category::Style, severity: Severity::Warning, span: Span::new(s, e), message: "Non-consecutive ring numbering", scope: Scope::Global }); break; }
            }
            last = Some(num);
        }
    }
}

fn ring_indices_sequence(bytes: &[u8]) -> Vec<(u32, usize, usize)> {
    let mut res = Vec::new();
    let mut i = 0usize;
    let n = bytes.len();
    let mut in_brkt = false;
    while i < n {
        let b = bytes[i];
        if b == b'[' { in_brkt = true; i += 1; continue; }
        if b == b']' { in_brkt = false; i += 1; continue; }
        if in_brkt { i += 1; continue; }
        if b == b'%' {
            if i + 2 < n && bytes[i + 1].is_ascii_digit() && bytes[i + 2].is_ascii_digit() {
                let num = ((bytes[i + 1] - b'0') as u32) * 10 + (bytes[i + 2] - b'0') as u32;
                res.push((num, i, i + 3));
                i += 3; continue;
            }
        } else if b.is_ascii_digit() {
            let num = (b - b'0') as u32;
            res.push((num, i, i + 1));
            i += 1; continue;
        }
        i += 1;
    }
    res
}

fn field_rank(f: &BracketField) -> u8 {
    match f { BracketField::Chiral(_) => 0, BracketField::HydrogenCount(_) => 1, BracketField::Charge(_) => 2, BracketField::Class(_) => 3 }
}

fn is_organic_subset(elem: umol_data::Element) -> bool {
    use umol_data::Element::*;
    matches!(elem, B | C | N | O | P | S | F | Cl | Br | I)
}

fn contains_sign_one(inner: &str) -> bool {
    let b = inner.as_bytes();
    let mut i = 0usize;
    while i < b.len() {
        let c = b[i];
        if c == b'+' || c == b'-' {
            let mut j = i + 1;
            let mut val: i32 = 0;
            let mut has_digit = false;
            while j < b.len() && b[j].is_ascii_digit() {
                has_digit = true;
                val = val.saturating_mul(10) + (b[j] - b'0') as i32;
                j += 1;
            }
            if has_digit && val == 1 { return true; }
            i = j; continue;
        }
        i += 1;
    }
    false
}

fn contains_h1(inner: &str) -> bool {
    let b = inner.as_bytes();
    let mut i = 0usize;
    while i + 1 < b.len() {
        if b[i] == b'H' && b[i + 1] == b'1' {
            if i + 2 >= b.len() || !b[i + 2].is_ascii_digit() { return true; }
        }
        i += 1;
    }
    false
}

fn chiral_param_out_of_range(inner: &str) -> bool {
    let b = inner.as_bytes();
    let mut i = 0usize;
    while i + 3 < b.len() {
        if b[i] == b'@' && b[i + 1] == b'T' && b[i + 2] == b'B' {
            let mut j = i + 3; let mut val: u32 = 0; let mut has_digit=false;
            while j < b.len() && b[j].is_ascii_digit() { has_digit=true; val = val.saturating_mul(10) + (b[j] - b'0') as u32; j += 1; }
            if has_digit && val > 20 { return true; }
            i = j; continue;
        }
        if b[i] == b'@' && b[i + 1] == b'O' && b[i + 2] == b'H' {
            let mut j = i + 3; let mut val: u32 = 0; let mut has_digit=false;
            while j < b.len() && b[j].is_ascii_digit() { has_digit=true; val = val.saturating_mul(10) + (b[j] - b'0') as u32; j += 1; }
            if has_digit && val > 30 { return true; }
            i = j; continue;
        }
        i += 1;
    }
    false
}
