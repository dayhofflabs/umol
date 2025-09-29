//! Linting for SMILES: collect diagnostics from lexing/parsing.

// Removed lalrpop dependency

mod context;
mod emitter;
mod registry;
mod rules;
mod utils;
pub use context::LintContext;
pub use emitter::{DiagnosticCandidate, Emitter, Scope};
pub use registry::{LintEngine, RuleRegistry};
pub use rules::Rule;

// Use FSM-based parser via io::smiles::parser
use super::parser::parse_smiles;
use crate::diagnostics::{Category, Code, Diagnostic, DiagnosticsReport, Severity, Span};
use crate::io::smiles::state::{ParseState, ParserMode};

// Initial linter for lexical/syntactic errors.
pub fn lint_smiles(input: &str) -> DiagnosticsReport {
    // Run through the rule engine (legacy rule preserves current behavior)
    let ctx = LintContext::new(input);
    let mut report = DiagnosticsReport::new();
    let mut registry = RuleRegistry::new();
    registry.register(&rules::LEX_ERRORS_RULE);
    registry.register(&rules::WS_RULE);
    registry.register(&rules::COMMENTS_RULE);
    registry.register(&rules::TRAILING_BOND_RULE);
    registry.register(&rules::DOT_RULES);
    registry.register(&rules::LEADING_BOND_RULE);
    registry.register(&rules::LEADING_RING_RULE);
    registry.register(&rules::CONSECUTIVE_BONDS_RULE);
    registry.register(&rules::TOP_GRP_TRAILING_RULE);
    registry.register(&rules::STYLE_PCT_RULE);
    registry.register(&rules::BRANCH_RULE);
    registry.register(&rules::BOND_STYLE_RULE);
    registry.register(&rules::GROUP_PARENS_RULE);
    registry.register(&rules::MIXED_AFTER_RULE);
    registry.register(&rules::RING_STYLE_RULE);
    registry.register(&rules::RING_NUMBERING_RULE);
    registry.register(&rules::RING_ERRORS_RULE);
    registry.register(&rules::BRACKET_RULE);
    registry.register(&rules::AROM_INCONSISTENT_RULE);
    let engine = LintEngine::new(registry);
    engine.run(&ctx, &mut report);
    report
}

// Run the parser and translate lalrpop errors to diagnostics; also runs the lexical/style lint.
pub fn lint_smiles_parse(input: &str) -> DiagnosticsReport {
    let mut report = lint_smiles(input);

    // Determine if we already have syntactic errors from pre-parse linting
    let has_syntactic_errors = report
        .diagnostics
        .iter()
        .any(|d| matches!(d.category, Category::Lex | Category::Bracket | Category::Branch | Category::Ring) && d.severity == Severity::Error);

    // If lexical errors present, parsing is likely to cascade; still attempt a parse for location.
    let mut state = ParseState::default();
    // Call FSM-based parser for location sanity; translate only generic failure
    if !has_syntactic_errors {
        let _ = parse_smiles(input.as_bytes());
    }
    // Merge any diagnostics emitted by the parser/state (e.g., ring rules)
    for d in state.diagnostics.into_iter() {
        // Allow only parser diagnostics that are not covered by pre-parse syntactic lints
        let c = d.code.0;
        let allow = match d.category {
            Category::Stereo | Category::Internal => true,
            Category::Ring => matches!(c, "RING_SELF_LOOP" | "RING_TWO_MEMBER" | "RING_BOND_DIR_CONFLICT" | "RING_BOND_ORDER_CONFLICT"),
            // Drop Lex/Bracket/Branch/Style/Num from parser in favor of lint-only
            Category::Lex | Category::Bracket | Category::Branch | Category::Style | Category::Num | Category::Syn => false,
        };
        if allow { report.push(d); }
    }
    report
}

// Experimental: parse in lint-fast mode (no IR, no late passes) and return diagnostics.
#[allow(dead_code)]
pub fn lint_smiles_parse_fast(input: &str) -> DiagnosticsReport { lint_smiles(input) }
