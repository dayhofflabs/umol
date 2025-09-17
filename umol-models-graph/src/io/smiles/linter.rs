//! Linting for SMILES: collect diagnostics from lexing/parsing.

use lalrpop_util::ParseError;

mod utils;
mod context;
mod emitter;
mod registry;
mod rules;
pub use context::LintContext;
pub use emitter::{DiagnosticCandidate, Emitter, Scope};
pub use registry::{LintEngine, RuleRegistry};
pub use rules::{Phase, Rule};

use super::lexer::Lexer;
use super::parser::grammar::MoleculeParser;
use crate::diagnostics::{Category, Code, Diagnostic, DiagnosticsReport, Span};
use crate::io::smiles::state::ParseState;

// Initial linter for lexical/syntactic errors.
pub fn lint_smiles(input: &str) -> DiagnosticsReport {
    // Run through the rule engine (legacy rule preserves current behavior)
    let ctx = LintContext::new(input);
    let mut report = DiagnosticsReport::new();
    let mut registry = RuleRegistry::new();
    registry.register(&rules::LEX_ERRORS_RULE);
    registry.register(&rules::WS_RULE);
    registry.register(&rules::TRAILING_BOND_RULE);
    registry.register(&rules::DOT_RULES);
    registry.register(&rules::STYLE_PCT_RULE);
    registry.register(&rules::BOND_STYLE_RULE);
    registry.register(&rules::RING_STYLE_RULE);
    registry.register(&rules::BRACKET_RULE);
    let engine = LintEngine::new(registry);
    engine.run(&ctx, &mut report);
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
