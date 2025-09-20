//! Lexical rules for SMILES linting.

use logos::Logos;

use super::{Phase, Rule, RuleMeta};
use crate::diagnostics::{Category, Code, Severity, Span};
use crate::io::smiles::lexer::Token;
use crate::io::smiles::linter::emitter::{DiagnosticCandidate, Emitter, Scope};
use crate::io::smiles::linter::LintContext;
use crate::io::smiles::segment::Segment;

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
        let logos_lexer = Token::lexer(ctx.input);
        for (res, span) in logos_lexer.spanned() {
            if res.is_err() {
                let slice = &ctx.input[span.start..span.end];
                if slice == "%" {
                    emit.candidate(DiagnosticCandidate {
                        code: Code("LEX_BAD_PERCENT_FORM"),
                        category: Category::Lex,
                        severity: Severity::Error,
                        span: Span::new(span.start, span.end),
                        message: "'%' not followed by two digits",
                        scope: Scope::Global,
                    });
                } else {
                    emit.candidate(DiagnosticCandidate {
                        code: Code("LEX_INVALID_TOKEN"),
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
        let segs = ctx.segments();
        if let Some(last_non_ws) = segs
            .iter()
            .rposition(|seg| !matches!(seg, Segment::WhitespaceBlock { .. }))
        {
            for (i, seg) in segs.iter().enumerate() {
                if i > last_non_ws {
                    break;
                }
                if let Segment::WhitespaceBlock { span } = seg {
                    emit.candidate(DiagnosticCandidate {
                        code: Code("LEX_INTERTOKEN_WHITESPACE"),
                        category: Category::Lex,
                        severity: Severity::Error,
                        span: *span,
                        message: "Inter-token whitespace is not allowed",
                        scope: Scope::Global,
                    });
                }
            }
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
        let mut last_bond_span = None;
        for seg in ctx.segments().iter() {
            match seg {
                Segment::WhitespaceBlock { .. } => {}
                Segment::Bond { span, .. } => last_bond_span = Some(*span),
                _ => last_bond_span = None,
            }
        }
        if let Some(span) = last_bond_span {
            emit.candidate(DiagnosticCandidate {
                code: Code("SYN_TRAILING_BOND"),
                category: Category::Syn,
                severity: Severity::Error,
                span,
                message: "Trailing bond symbol",
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
        // Dot before ring
        let mut last_dot: Option<Span> = None;
        for seg in ctx.segments().iter() {
            match seg {
                Segment::WhitespaceBlock { .. } => {}
                Segment::ComponentSeparator { span } => last_dot = Some(*span),
                Segment::RingClosure { span, .. } => {
                    if let Some(dot) = last_dot.take() {
                        if dot.end == span.start {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("SYN_DOT_BEFORE_RING"),
                                category: Category::Syn,
                                severity: Severity::Error,
                                span: Span::new(dot.start, span.end),
                                message: "Dot before ring index is invalid",
                                scope: Scope::Global,
                            });
                        }
                    }
                }
                _ => last_dot = None,
            }
        }
        // Dot positions: leading, trailing, multiple
        let segs = ctx.segments();
        if let Some(i) = segs
            .iter()
            .position(|seg| !matches!(seg, Segment::WhitespaceBlock { .. }))
        {
            if let Segment::ComponentSeparator { span } = segs[i] {
                emit.candidate(DiagnosticCandidate {
                    code: Code("SYN_LEADING_DOT"),
                    category: Category::Syn,
                    severity: Severity::Error,
                    span,
                    message: "Leading dot",
                    scope: Scope::Global,
                });
            }
        }
        if let Some(i) = segs
            .iter()
            .rposition(|seg| !matches!(seg, Segment::WhitespaceBlock { .. }))
        {
            if let Segment::ComponentSeparator { span } = segs[i] {
                emit.candidate(DiagnosticCandidate {
                    code: Code("SYN_TRAILING_DOT"),
                    category: Category::Syn,
                    severity: Severity::Error,
                    span,
                    message: "Trailing dot",
                    scope: Scope::Global,
                });
            }
        }
        for w in segs.windows(2) {
            if let [Segment::ComponentSeparator { span: s1 }, Segment::ComponentSeparator { span: s2 }] = w
            {
                emit.candidate(DiagnosticCandidate {
                    code: Code("SYN_MULTIPLE_DOTS"),
                    category: Category::Syn,
                    severity: Severity::Error,
                    span: Span::new(s1.start, s2.end),
                    message: "Multiple dots",
                    scope: Scope::Global,
                });
                break;
            }
        }
    }
}
pub static DOT_RULES: DotRules = DotRules;
