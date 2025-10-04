use crate::io::smiles::diagnostics::{Category, Diagnostic, DiagnosticCode, Severity, Span};
use crate::io::smiles::ParseError;

pub fn map_parse_error(err: &ParseError, input: &str) -> Diagnostic {
    let input_len = input.len();
    // match *err {
    //     // Ring diagnostics
    //     ParseError::UnbalancedRingIndex { open_pos } => Diagnostic {
    //         code: DiagnosticCode::RingUnclosed,
    //         category: Category::Ring,
    //         severity: Severity::Error,
    //         span: Span::new(open_pos, open_pos),
    //         message: "Ring index opened but not closed",
    //         details: None,
    //     },
    //     ParseError::SelfLoopRing { pos } => Diagnostic {
    //         code: DiagnosticCode::TopoSelfLoop,
    //         category: Category::Ring,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Ring closure creates a self-loop",
    //         details: None,
    //     },
    //     ParseError::TwoMemberRing { pos } => Diagnostic {
    //         code: DiagnosticCode::TopoParallelEdges,
    //         category: Category::Ring,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Ring closure creates a two-member ring",
    //         details: None,
    //     },
    //     ParseError::MultipleRings { pos } => Diagnostic {
    //         code: DiagnosticCode::TopoParallelEdges,
    //         category: Category::Ring,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Multiple ring closures between the same atom pair",
    //         details: None,
    //     },
    //     ParseError::RingBondDirConflict { pos, .. } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Ring,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Conflicting ring bond directions",
    //         details: None,
    //     },
    //     ParseError::RingBondOrderConflict { pos, .. } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Ring,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Conflicting ring bond orders",
    //         details: None,
    //     },

    //     // Bracket balance and fields
    //     ParseError::UnbalancedOpenBracket { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Bracket,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Unbalanced '[' bracket",
    //         details: None,
    //     },
    //     ParseError::UnbalancedCloseBracket { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Bracket,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Unbalanced ']' bracket",
    //         details: None,
    //     },
    //     ParseError::InvalidBracket { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Bracket,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Invalid bracket atom",
    //         details: None,
    //     },
    //     ParseError::MissingClassIndex { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Bracket,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Empty atom class field",
    //         details: None,
    //     },
    //     ParseError::ChiralityOutOfRange { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Bracket,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Chirality descriptor parameter out of range",
    //         details: None,
    //     },
    //     ParseError::DuplicateBracketField { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Bracket,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Duplicate bracket field",
    //         details: None,
    //     },
    //     // FIX THIS
    //     ParseError::BracketHwithH { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Bracket,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Hydrogen field on hydrogen element",
    //         details: None,
    //     },
    //     // FIX THIS
    //     ParseError::StrayBracketField { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Bracket,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Bracket-only field outside brackets",
    //         details: None,
    //     },

    //     // Branch/group
    //     // FIX THIS
    //     ParseError::UnbalancedOpenParen { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Branch,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Open '(' not closed",
    //         details: None,
    //     },
    //     // FIX THIS
    //     ParseError::UnbalancedCloseParen { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Branch,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Unexpected ')'",
    //         details: None,
    //     },
    //     // FIX THIS
    //     ParseError::EmptyBranch { pos } | ParseError::EmptyGroup { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Branch,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Empty branch/group",
    //         details: None,
    //     },
    //     // FIX THIS
    //     ParseError::NonfinalGroup { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Branch,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Group trailing content",
    //         details: None,
    //     },
    //     // FIX THIS
    //     ParseError::GroupLeadingDot { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Branch,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Group begins with a dot",
    //         details: None,
    //     },
    //     // FIX THIS
    //     ParseError::GroupLeadingBond { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Branch,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Group begins with a bond",
    //         details: None,
    //     },

    //     // Bond/dot
    //     // FIX THIS
    //     ParseError::LeadingBondSymbol { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Lex,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Leading bond token",
    //         details: None,
    //     },
    //     // FIX THIS
    //     ParseError::ConsecutiveBonds { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Lex,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Consecutive bond tokens",
    //         details: None,
    //     },
    //     // FIX THIS
    //     ParseError::TrailingBondSymbol { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Lex,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Trailing bond token",
    //         details: None,
    //     },
    //     // FIX THIS
    //     ParseError::LeadingDot { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Lex,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Leading dot",
    //         details: None,
    //     },
    //     // FIX THIS
    //     ParseError::ConsecutiveDots { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Lex,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Consecutive dots",
    //         details: None,
    //     },
    //     // FIX THIS
    //     ParseError::TrailingDot { pos } => Diagnostic {
    //         code: DiagnosticCode::ParseError,
    //         category: Category::Lex,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Trailing dot",
    //         details: None,
    //     },

    //     // Ring leading and index invalid
    //     // FIX THIS
    //     ParseError::LeadingRing { pos } => {
    //         if pos > 0 && input.as_bytes()[pos - 1] == b'.' {
    //             Diagnostic {
    //                 code: DiagnosticCode::ParseError,
    //                 category: Category::Lex,
    //                 severity: Severity::Error,
    //                 span: Span::new(pos - 1, pos),
    //                 message: "Leading dot",
    //                 details: None,
    //             }
    //         } else {
    //             Diagnostic {
    //                 code: DiagnosticCode::ParseError,
    //                 category: Category::Lex,
    //                 severity: Severity::Error,
    //                 span: Span::new(pos, pos.saturating_add(1)),
    //                 message: "Leading ring index",
    //                 details: None,
    //             }
    //         }
    //     }
    //     ParseError::InvalidRingIndex { pos } => Diagnostic {
    //         code: DiagnosticCode::RingIndexInvalid,
    //         category: Category::Lex,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Invalid percent ring index",
    //         details: None,
    //     },

    //     // Whitespace/comments (strict mode)
    //     ParseError::InvalidWhitespace { pos } => Diagnostic {
    //         code: DiagnosticCode::LexIntertokenWhitespace,
    //         category: Category::Lex,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Inter-token whitespace is not allowed",
    //         details: None,
    //     },
    //     ParseError::InvalidComment { pos } => Diagnostic {
    //         code: DiagnosticCode::InvalidComment,
    //         category: Category::Lex,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(2)),
    //         message: "Comments are not allowed",
    //         details: None,
    //     },
    //     ParseError::UnterminatedBlockComment { pos } => Diagnostic {
    //         code: DiagnosticCode::LexUnterminatedBlockComment,
    //         category: Category::Lex,
    //         severity: Severity::Error,
    //         span: Span::new(pos, input_len),
    //         message: "Unterminated block comment",
    //         details: None,
    //     },

    //     // Generic tokens
    //     ParseError::UnsupportedToken { pos } => Diagnostic {
    //         code: DiagnosticCode::UnsupportedToken,
    //         category: Category::Lex,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Invalid token",
    //         details: None,
    //     },
    //     ParseError::InvalidElement { pos } => Diagnostic {
    //         code: DiagnosticCode::LexInvalidElement,
    //         category: Category::Lex,
    //         severity: Severity::Error,
    //         span: Span::new(pos, pos.saturating_add(1)),
    //         message: "Invalid element symbol",
    //         details: None,
    //     },
    // }
    Diagnostic {
        code: DiagnosticCode::InternalError,
        category: Category::Lex,
        severity: Severity::Error,
        span: Span::new(0, input_len),
        message: "Invalid token",
        details: None,
    }
}


