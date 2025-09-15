//! Lexical rules for SMILES linting.

use logos::Logos;

use super::{Phase, Rule, RuleMeta};
use crate::diagnostics::{Category, Severity, Span};
use crate::io::smiles::linter::bracket::{
    lint_dot_before_ring, lint_dot_positions, lint_intertoken_whitespace, lint_trailing_bond,
};
use crate::io::smiles::linter::emitter::{DiagnosticCandidate, Emitter, Scope};
use crate::io::smiles::linter::LintContext;

pub struct LexErrorsRule;
static META_LEX: RuleMeta = RuleMeta {
    id: "LEX_ERRORS",
    category: Category::Lex,
    default_severity: Severity::Error,
};
impl Rule for LexErrorsRule {
    fn meta(&self) -> &'static RuleMeta {
        &META_LEX
    }
    fn phase(&self) -> Phase {
        Phase::Lex
    }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let logos_lexer = crate::io::smiles::lexer::Token::lexer(ctx.input);
        for (res, span) in logos_lexer.spanned() {
            if res.is_err() {
                let slice = &ctx.input[span.start..span.end];
                if slice == "%" {
                    emit.candidate(DiagnosticCandidate {
                        code: crate::diagnostics::Code("LEX_BAD_PERCENT_FORM"),
                        category: Category::Lex,
                        severity: Severity::Error,
                        span: Span::new(span.start, span.end),
                        message: "'%' not followed by two digits",
                        scope: Scope::Global,
                    });
                } else {
                    emit.candidate(DiagnosticCandidate {
                        code: crate::diagnostics::Code("LEX_INVALID_TOKEN"),
                        category: Category::Lex,
                        severity: Severity::Error,
                        span: Span::new(span.start, span.end),
                        message: "Invalid token",
                        scope: Scope::Global,
                    });
                }
            }
        }
    }
}
pub static LEX_ERRORS_RULE: LexErrorsRule = LexErrorsRule;

pub struct WhitespaceRule;
static META_WS: RuleMeta = RuleMeta {
    id: "LEX_WS",
    category: Category::Lex,
    default_severity: Severity::Error,
};
impl Rule for WhitespaceRule {
    fn meta(&self) -> &'static RuleMeta {
        &META_WS
    }
    fn phase(&self) -> Phase {
        Phase::Lex
    }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let mut tmp = crate::diagnostics::DiagnosticsReport::new();
        lint_intertoken_whitespace(ctx.input, &mut tmp);
        for d in tmp.diagnostics {
            emit.candidate(DiagnosticCandidate {
                code: d.code,
                category: d.category,
                severity: d.severity,
                span: d.span,
                message: Box::<str>::leak(d.message.into()),
                scope: Scope::Global,
            });
        }
    }
}
pub static WS_RULE: WhitespaceRule = WhitespaceRule;

pub struct TrailingBondRule;
static META_TB: RuleMeta = RuleMeta {
    id: "SYN_TB",
    category: Category::Syn,
    default_severity: Severity::Error,
};
impl Rule for TrailingBondRule {
    fn meta(&self) -> &'static RuleMeta {
        &META_TB
    }
    fn phase(&self) -> Phase {
        Phase::Lex
    }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let mut tmp = crate::diagnostics::DiagnosticsReport::new();
        lint_trailing_bond(ctx.input, &mut tmp);
        for d in tmp.diagnostics {
            emit.candidate(DiagnosticCandidate {
                code: d.code,
                category: d.category,
                severity: d.severity,
                span: d.span,
                message: Box::<str>::leak(d.message.into()),
                scope: Scope::Global,
            });
        }
    }
}
pub static TRAILING_BOND_RULE: TrailingBondRule = TrailingBondRule;

pub struct DotRules;
static META_DOT: RuleMeta = RuleMeta {
    id: "SYN_DOTS",
    category: Category::Syn,
    default_severity: Severity::Error,
};
impl Rule for DotRules {
    fn meta(&self) -> &'static RuleMeta {
        &META_DOT
    }
    fn phase(&self) -> Phase {
        Phase::Lex
    }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let mut tmp = crate::diagnostics::DiagnosticsReport::new();
        lint_dot_before_ring(ctx.input, &mut tmp);
        lint_dot_positions(ctx.input, &mut tmp);
        for d in tmp.diagnostics {
            emit.candidate(DiagnosticCandidate {
                code: d.code,
                category: d.category,
                severity: d.severity,
                span: d.span,
                message: Box::<str>::leak(d.message.into()),
                scope: Scope::Global,
            });
        }
    }
}
pub static DOT_RULES: DotRules = DotRules;
