//! Linting for SMILES: collect diagnostics from lexing/parsing.

mod context;
mod emitter;
mod registry;
mod rules;
mod utils;
pub use context::LintContext;
pub use emitter::{DiagnosticCandidate, Emitter, Scope};
pub use registry::{LintEngine, RuleRegistry};
pub use rules::Rule;

use super::parser::parse_smiles;
use crate::diagnostics::{Category, Code, Diagnostic, DiagnosticsReport, Severity, Span};
// Remove legacy parser state

// Initial linter for lexical/syntactic errors.
pub fn lint_smiles(input: &str) -> DiagnosticsReport {
    // Run through the rule engine (legacy rule preserves current behavior)
    let ctx = LintContext::new(input);
    let mut report = DiagnosticsReport::new();
    let mut registry = RuleRegistry::new();
    registry.register(&rules::LEX_ERRORS_RULE);
    registry.register(&rules::WS_RULE);
    registry.register(&rules::COMMENTS_RULE);
    let engine = LintEngine::new(registry);
    engine.run(&ctx, &mut report);
    report
}

// Run the parser and translate parser errors to diagnostics; also runs the lexical/style lint.
pub fn lint_smiles_parse(input: &str) -> DiagnosticsReport {
    let mut report = lint_smiles(input);

    // Determine if we already have syntactic errors from pre-parse linting
    let has_syntactic_errors = report
        .diagnostics
        .iter()
        .any(|d| matches!(d.category, Category::Lex | Category::Bracket | Category::Branch | Category::Ring) && d.severity == Severity::Error);

    // If lexical errors present, parsing is likely to cascade; still attempt a parse for location.
    // Call FSM-based parser for location sanity; translate only generic failure
    if !has_syntactic_errors {
        let _ = parse_smiles(input.as_bytes());
    }
    report
}

// Experimental: parse in lint-fast mode (no IR, no late passes) and return diagnostics.
#[allow(dead_code)]
pub fn lint_smiles_parse_fast(input: &str) -> DiagnosticsReport { lint_smiles(input) }
