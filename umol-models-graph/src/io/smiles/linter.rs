//! Linting for SMILES: collect diagnostics from lexing/parsing.

mod context;
mod emitter;
pub use context::LintContext;
pub use emitter::{Emitter, LintConfig};

use super::checker::{
    check_aromaticity, check_stereo_double, check_topology, check_valence, AromaticityConfig,
    AromaticityModel, ValenceConfig, ValenceModel,
};
use super::diagnostics::{Category, Code, Diagnostic, DiagnosticsReport, Severity, Span};
use super::parser::parse_smiles;
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
                emitter.emit(Diagnostic {
                    code: Code("RING_UNCLOSED"),
                    category: Category::Ring,
                    severity: Severity::Error,
                    span: Span::new(open_pos, open_pos),
                    message: "Ring index opened but not closed",
                    details: None,
                });
            }
            ParseError::RingSelfLoop { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("RING_SELF_LOOP"),
                    category: Category::Ring,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Ring closure creates a self-loop",
                    details: None,
                });
            }
            ParseError::RingTwoMember { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("RING_TWO_MEMBER"),
                    category: Category::Ring,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Ring closure creates a two-member ring",
                    details: None,
                });
            }
            ParseError::RingMultipleRings { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("RING_MULTIPLE_RINGS"),
                    category: Category::Ring,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Multiple ring closures between the same atom pair",
                    details: None,
                });
            }
            ParseError::RingBondDirConflict { pos, .. } => {
                emitter.emit(Diagnostic {
                    code: Code("RING_BOND_DIR_CONFLICT"),
                    category: Category::Ring,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Conflicting ring bond directions",
                    details: None,
                });
            }
            ParseError::RingBondOrderConflict { pos, .. } => {
                emitter.emit(Diagnostic {
                    code: Code("RING_BOND_ORDER_CONFLICT"),
                    category: Category::Ring,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Conflicting ring bond orders",
                    details: None,
                });
            }

            // Bracket balance
            ParseError::UnbalancedOpenBracket { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("BRKT_UNBALANCED_OPEN"),
                    category: Category::Bracket,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Unbalanced '[' bracket",
                    details: None,
                });
            }
            ParseError::UnbalancedCloseBracket { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("BRKT_UNBALANCED_CLOSE"),
                    category: Category::Bracket,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Unbalanced ']' bracket",
                    details: None,
                });
            }
            ParseError::InvalidBracket { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("BRKT_INVALID"),
                    category: Category::Bracket,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Invalid bracket atom",
                    details: None,
                });
            }
            ParseError::BracketEmptyClass { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("BRKT_EMPTY_CLASS"),
                    category: Category::Bracket,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Empty atom class field",
                    details: None,
                });
            }
            ParseError::BracketChiralityOutOfRange { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("BRKT_CHIRAL_OUT_OF_RANGE"),
                    category: Category::Bracket,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Chirality descriptor parameter out of range",
                    details: None,
                });
            }
            ParseError::BracketDuplicateField { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("BRKT_DUP_FIELD"),
                    category: Category::Bracket,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Duplicate bracket field",
                    details: None,
                });
            }
            ParseError::BracketHOnH { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("BRKT_H_ON_H"),
                    category: Category::Bracket,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Hydrogen field on hydrogen element",
                    details: None,
                });
            }

            // Branch/group
            ParseError::UnbalancedBranchOpen { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("BRCH_UNCLOSED"),
                    category: Category::Branch,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Open '(' not closed",
                    details: None,
                });
            }
            ParseError::UnbalancedBranchClose { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("BRCH_UNEXPECTED_CLOSE"),
                    category: Category::Branch,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Unexpected ')'",
                    details: None,
                });
            }
            ParseError::EmptyBranch { pos } | ParseError::EmptyGroup { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("BRCH_EMPTY_BRANCH"),
                    category: Category::Branch,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Empty branch/group",
                    details: None,
                });
            }

            // Bond/dot
            ParseError::LeadingBond { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("LEX_LEADING_BOND"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Leading bond token",
                    details: None,
                });
            }
            ParseError::ConsecutiveBond { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("LEX_CONSECUTIVE_BONDS"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Consecutive bond tokens",
                    details: None,
                });
            }
            ParseError::TrailingBond { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("LEX_TRAILING_BOND"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Trailing bond token",
                    details: None,
                });
            }
            ParseError::LeadingDot { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("LEX_LEADING_DOT"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Leading dot",
                    details: None,
                });
            }
            ParseError::ConsecutiveDot { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("LEX_MULTIPLE_DOTS"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Consecutive dots",
                    details: None,
                });
            }
            ParseError::TrailingDot { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("LEX_TRAILING_DOT"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Trailing dot",
                    details: None,
                });
            }

            // Ring leading and index invalid
            ParseError::LeadingRing { pos } => {
                if pos > 0 && input.as_bytes()[pos - 1] == b'.' {
                    // Treat as leading dot instead of dot-before-ring
                    emitter.emit(Diagnostic {
                        code: Code("LEX_LEADING_DOT"),
                        category: Category::Lex,
                        severity: Severity::Error,
                        span: Span::new(pos - 1, pos),
                        message: "Leading dot",
                        details: None,
                    });
                } else {
                    emitter.emit(Diagnostic {
                        code: Code("LEX_LEADING_RING"),
                        category: Category::Lex,
                        severity: Severity::Error,
                        span: Span::new(pos, pos.saturating_add(1)),
                        message: "Leading ring index",
                        details: None,
                    });
                }
            }
            ParseError::RingIndexInvalid { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("LEX_RING_INDEX_INVALID"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Invalid percent ring index",
                    details: None,
                });
            }

            // Whitespace/comments (strict mode)
            ParseError::InvalidWhitespace { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("LEX_INTERTOKEN_WHITESPACE"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Inter-token whitespace is not allowed",
                    details: None,
                });
            }
            ParseError::InvalidComment { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("LEX_COMMENT"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(2)),
                    message: "Comments are not allowed",
                    details: None,
                });
            }
            ParseError::UnterminatedBlockComment { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("LEX_UNTERMINATED_BLOCK_COMMENT"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(pos, input.len()),
                    message: "Unterminated block comment",
                    details: None,
                });
            }
            ParseError::TopLevelGroupTrailing { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("BRCH_GROUP_TRAILING"),
                    category: Category::Branch,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Group trailing content",
                    details: None,
                });
            }
            ParseError::UnsupportedToken { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("LEX_INVALID_TOKEN"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Invalid token",
                    details: None,
                });
            }
            ParseError::InvalidElement { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("LEX_INVALID_ELEMENT"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Invalid element symbol",
                    details: None,
                });
            }
            ParseError::FieldOutsideBracket { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("BRKT_FIELD_OUTSIDE"),
                    category: Category::Bracket,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Bracket-only field outside brackets",
                    details: None,
                });
            }
            ParseError::GroupLeadingDot { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("GRP_LEADING_DOT"),
                    category: Category::Branch,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Group begins with a dot",
                    details: None,
                });
            }
            ParseError::GroupLeadingBond { pos } => {
                emitter.emit(Diagnostic {
                    code: Code("GRP_LEADING_BOND"),
                    category: Category::Branch,
                    severity: Severity::Error,
                    span: Span::new(pos, pos.saturating_add(1)),
                    message: "Group begins with a bond",
                    details: None,
                });
            }
        }
    }

    // Style and numeric advisories
    run_style_and_numeric_checks(&ctx, &mut emitter, parse_res.is_ok());

    // Post-parse: topology + stereo checks when parsing succeeds (after emitter releases &mut report)
    if let Ok(ref mol) = parse_res {
        let input_len = input.len();
        check_topology(mol, None, &mut report, input_len);
        let v_cfg = ValenceConfig::default();
        let v_model = ValenceModel::simple_organic();
        check_valence(mol, None, &mut report, input_len, &v_model, &v_cfg);
        check_stereo_double(mol, None, &mut report, input_len);
        // Aromaticity verification scaffold (HMO/Clar config only)
        let a_cfg = AromaticityConfig::default();
        let a_model = AromaticityModel::default();
        let _ = check_aromaticity(mol, None, &mut report, input_len, &a_model, &a_cfg);
    }

    report
}

fn run_style_and_numeric_checks(ctx: &LintContext, emit: &mut Emitter, only_when_parse_ok: bool) {
    let input = ctx.input;
    let bytes = input.as_bytes();

    // Percent ring index style: %01..%09
    let mut i = 0usize;
    while i + 2 < bytes.len() {
        if bytes[i] == b'%'
            && bytes[i + 1] == b'0'
            && (bytes[i + 2] >= b'1' && bytes[i + 2] <= b'9')
        {
            emit.emit(Diagnostic {
                code: Code("STYLE_UNNECESSARY_PERCENT_RING_INDEX"),
                category: Category::Style,
                severity: Severity::Warning,
                span: Span::new(i, i + 3),
                message: "Prefer single-digit ring index for 1..9",
                details: None,
            });
        }
        i += 1;
    }

    // Ring numbering style checks only when parsing succeeded to avoid noise
    if only_when_parse_ok {
        let seq = ring_indices_sequence(bytes);
        if let Some((first_num, s, e)) = seq.first().copied() {
            if first_num != 1 {
                emit.emit(Diagnostic {
                    code: Code("STYLE_FIRST_RING_NOT_ONE"),
                    category: Category::Style,
                    severity: Severity::Warning,
                    span: Span::new(s, e),
                    message: "Prefer starting ring numbering at 1",
                    details: None,
                });
            }
        }
        // Non-consecutive jumps
        let mut last: Option<u32> = None;
        for (num, s, e) in seq.into_iter() {
            if let Some(p) = last {
                if num > p + 1 {
                    emit.emit(Diagnostic {
                        code: Code("STYLE_NONCONSECUTIVE_RING_NUMBERING"),
                        category: Category::Style,
                        severity: Severity::Warning,
                        span: Span::new(s, e),
                        message: "Non-consecutive ring numbering",
                        details: None,
                    });
                    break;
                }
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
        if b == b'[' {
            in_brkt = true;
            i += 1;
            continue;
        }
        if b == b']' {
            in_brkt = false;
            i += 1;
            continue;
        }
        if in_brkt {
            i += 1;
            continue;
        }
        if b == b'%' {
            if i + 2 < n && bytes[i + 1].is_ascii_digit() && bytes[i + 2].is_ascii_digit() {
                let num = ((bytes[i + 1] - b'0') as u32) * 10 + (bytes[i + 2] - b'0') as u32;
                res.push((num, i, i + 3));
                i += 3;
                continue;
            }
        } else if b.is_ascii_digit() {
            let num = (b - b'0') as u32;
            res.push((num, i, i + 1));
            i += 1;
            continue;
        }
        i += 1;
    }
    res
}
