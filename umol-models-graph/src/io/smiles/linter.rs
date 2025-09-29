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
                    code: Code("BRKT_UNBALANCED_OPEN"),
                    category: Category::Bracket,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Unbalanced '[' bracket",
                    scope: Scope::Bracket { start: pos, end: pos.saturating_add(1) },
                });
            }
            ParseError::UnbalancedCloseBracket { pos } => {
                emitter.candidate(DiagnosticCandidate {
                    code: Code("BRKT_UNBALANCED_CLOSE"),
                    category: Category::Bracket,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Unbalanced ']' bracket",
                    scope: Scope::Bracket { start: pos, end: pos.saturating_add(1) },
                });
            }
            ParseError::InvalidBracket { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("BRKT_INVALID"), category: Category::Bracket, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Invalid bracket atom", scope: Scope::Global });
            }
            ParseError::BracketHCountTwoDigits { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("BRKT_HCOUNT_TWO_DIGITS"), category: Category::Bracket, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "H-count must be one digit", scope: Scope::Global });
            }
            ParseError::BracketEmptyClass { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("BRKT_EMPTY_CLASS"), category: Category::Bracket, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Empty atom class field", scope: Scope::Global });
            }
            ParseError::BracketDuplicateField { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("BRKT_DUP_FIELD"), category: Category::Bracket, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Duplicate bracket field", scope: Scope::Global });
            }
            ParseError::BracketHOnH { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("BRKT_H_ON_H"), category: Category::Bracket, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Hydrogen field on hydrogen element", scope: Scope::Global });
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
                    // Treat as leading dot instead of dot-before-ring
                    emitter.candidate(DiagnosticCandidate { code: Code("LEX_LEADING_DOT"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos - 1, pos), message: "Leading dot", scope: Scope::Global });
                } else {
                    emitter.candidate(DiagnosticCandidate { code: Code("LEX_LEADING_RING"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Leading ring index", scope: Scope::Global });
                }
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
            ParseError::TopLevelGroupTrailing { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("BRCH_GROUP_TRAILING"), category: Category::Branch, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Group trailing content", scope: Scope::Global });
            }
            ParseError::UnsupportedToken { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("LEX_INVALID_TOKEN"), category: Category::Lex, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Invalid token", scope: Scope::Global });
            }
            ParseError::FieldOutsideBracket { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("BRKT_FIELD_OUTSIDE"), category: Category::Bracket, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Bracket-only field outside brackets", scope: Scope::Global });
            }
            ParseError::GroupLeadingConnector { pos } => {
                emitter.candidate(DiagnosticCandidate { code: Code("GRP_LEADING_CONNECTOR"), category: Category::Branch, severity: Severity::Error, span: Span::new(pos, pos.saturating_add(1)), message: "Group begins with a connector", scope: Scope::Global });
            }
        }
    }

    // Style and numeric advisories
    run_style_and_numeric_checks(&ctx, &mut emitter, parse_res.is_ok());

    emitter.flush();
    report
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

    // Bracket-based style/num checks are deferred to a later pass to avoid re-parsing.

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

// Bracket-style helper functions removed along with bracket re-parsing.
