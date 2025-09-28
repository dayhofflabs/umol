//! Lexical rules for SMILES linting.

use logos::Logos;

use super::{Rule, RuleMeta};
use crate::diagnostics::{Category, Code, Severity, Span};
use crate::io::smiles::iterators::{BondKind, Segment};
use crate::io::smiles::lexer::Token;
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
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let logos_lexer = Token::lexer(ctx.input.as_bytes());
        for (res, span) in logos_lexer.spanned() {
            if res.is_err() {
                let slice = &ctx.input[span.start..span.end];
                if slice == "%" {
                    emit.candidate(DiagnosticCandidate {
                        code: Code("LEX_RING_INDEX_INVALID"),
                        category: Category::Lex,
                        severity: Severity::Error,
                        span: Span::new(span.start, span.end),
                        message: "Percent ring index must be exactly two digits",
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
                code: Code("LEX_TRAILING_BOND"),
                category: Category::Lex,
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
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        // Dot before ring
        let mut last_dot: Option<Span> = None;
        for seg in ctx.segments().iter() {
            match seg {
                Segment::WhitespaceBlock { .. } => {}
                Segment::NewComponent { span } => last_dot = Some(*span),
                Segment::RingClosure { span, .. } => {
                    if let Some(dot) = last_dot.take() {
                        if dot.end == span.start {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("LEX_DOT_BEFORE_RING"),
                                category: Category::Lex,
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
            if let Segment::NewComponent { span } = segs[i] {
                emit.candidate(DiagnosticCandidate {
                    code: Code("LEX_LEADING_DOT"),
                    category: Category::Lex,
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
            if let Segment::NewComponent { span } = segs[i] {
                emit.candidate(DiagnosticCandidate {
                    code: Code("LEX_TRAILING_DOT"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span,
                    message: "Trailing dot",
                    scope: Scope::Global,
                });
            }
        }
        for w in segs.windows(2) {
            if let [Segment::NewComponent { span: s1 }, Segment::NewComponent { span: s2 }] =
                w
            {
                emit.candidate(DiagnosticCandidate {
                    code: Code("LEX_MULTIPLE_DOTS"),
                    category: Category::Lex,
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

// Leading bond
pub struct LeadingBondRule;
static META_LEX_LEAD_BOND: RuleMeta = RuleMeta { id: "LEX_LEADING_BOND_RULE", category: Category::Lex, default_severity: Severity::Error };
impl Rule for LeadingBondRule {
    fn meta(&self) -> &'static RuleMeta { &META_LEX_LEAD_BOND }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        if let Some(seg) = ctx.segments().iter().find(|s| !matches!(s, Segment::WhitespaceBlock { .. })) {
            if let Segment::Bond { span, .. } = seg {
                emit.candidate(DiagnosticCandidate { code: Code("LEX_LEADING_BOND"), category: Category::Lex, severity: Severity::Error, span: *span, message: "Leading bond symbol", scope: Scope::Global });
            }
        }
    }
}
pub static LEADING_BOND_RULE: LeadingBondRule = LeadingBondRule;

// Leading ring
pub struct LeadingRingRule;
static META_LEX_LEAD_RING: RuleMeta = RuleMeta { id: "LEX_LEADING_RING_RULE", category: Category::Lex, default_severity: Severity::Error };
impl Rule for LeadingRingRule {
    fn meta(&self) -> &'static RuleMeta { &META_LEX_LEAD_RING }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        if let Some(seg) = ctx.segments().iter().find(|s| !matches!(s, Segment::WhitespaceBlock { .. })) {
            if let Segment::RingClosure { span, .. } = seg {
                emit.candidate(DiagnosticCandidate { code: Code("LEX_LEADING_RING"), category: Category::Lex, severity: Severity::Error, span: *span, message: "Leading ring index", scope: Scope::Global });
            }
        }
    }
}
pub static LEADING_RING_RULE: LeadingRingRule = LeadingRingRule;

// Consecutive bonds
pub struct ConsecutiveBondsRule;
static META_LEX_CONS_BONDS: RuleMeta = RuleMeta { id: "LEX_CONSEC_BONDS_RULE", category: Category::Lex, default_severity: Severity::Error };
impl Rule for ConsecutiveBondsRule {
    fn meta(&self) -> &'static RuleMeta { &META_LEX_CONS_BONDS }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let segs = ctx.segments();
        let mut last_bond_span: Option<Span> = None;
        for seg in segs.iter() {
            match seg {
                Segment::WhitespaceBlock { .. } => {}
                Segment::Bond { span, .. } => {
                    if let Some(prev) = last_bond_span.take() {
                        emit.candidate(DiagnosticCandidate { code: Code("LEX_CONSECUTIVE_BONDS"), category: Category::Lex, severity: Severity::Error, span: Span::new(prev.start, span.end), message: "Consecutive bond symbols", scope: Scope::Global });
                    } else {
                        last_bond_span = Some(*span);
                    }
                }
                _ => last_bond_span = None,
            }
        }
    }
}
pub static CONSECUTIVE_BONDS_RULE: ConsecutiveBondsRule = ConsecutiveBondsRule;

// Top-level group trailing
pub struct TopGroupTrailingRule;
static META_TOP_GRP: RuleMeta = RuleMeta { id: "LEX_TOP_GRP_TRAILING_RULE", category: Category::Lex, default_severity: Severity::Error };
impl Rule for TopGroupTrailingRule {
    fn meta(&self) -> &'static RuleMeta { &META_TOP_GRP }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let segs = ctx.segments();
        let mut depth: i32 = 0;
        for idx in 0..segs.len() {
            match segs[idx] {
                Segment::WhitespaceBlock { .. } => {}
                Segment::BranchOpen { .. } => depth += 1,
                Segment::BranchClose { span } => {
                    depth -= 1;
                    if depth == 0 {
                        let mut j = idx + 1;
                        while j < segs.len() && matches!(segs[j], Segment::WhitespaceBlock { .. }) { j += 1; }
                        if j < segs.len() && !matches!(segs[j], Segment::NewComponent { .. }) {
                            emit.candidate(DiagnosticCandidate { code: Code("TOP_GRP_TRAILING"), category: Category::Lex, severity: Severity::Error, span, message: "Top-level group followed by non-dot", scope: Scope::Global });
                        }
                    }
                }
                _ => {}
            }
        }
    }
}
pub static TOP_GRP_TRAILING_RULE: TopGroupTrailingRule = TopGroupTrailingRule;

pub struct InconsistentAromaticityRule;
static META_AROM: RuleMeta = RuleMeta {
    id: "LEX_AROM_CONSISTENCY",
    category: Category::Lex,
    default_severity: Severity::Warning,
};

fn is_aromatic_simple(raw: &str) -> bool {
    matches!(raw, "b" | "c" | "n" | "o" | "p" | "s" | "se" | "as")
}

fn is_aromatic_bracket(inner: &str) -> bool {
    // Very lightweight check: skip leading digits, then look at the element symbol start(s)
    let bytes = inner.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() && bytes[i].is_ascii_digit() { i += 1; }
    if i >= bytes.len() { return false; }
    // Two-letter aromatic: se, as
    if i + 1 < bytes.len() {
        let two = &inner[i..i + 2].to_ascii_lowercase();
        if two == "se" || two == "as" { return true; }
    }
    // Single-letter aromatic: b,c,n,o,p,s
    let ch = bytes[i] as char;
    matches!(ch, 'b' | 'c' | 'n' | 'o' | 'p' | 's')
}

impl Rule for InconsistentAromaticityRule {
    fn meta(&self) -> &'static RuleMeta { &META_AROM }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let segs = ctx.segments();
        for i in 0..segs.len() {
            let (span, kind) = match segs[i] {
                Segment::Bond { span, kind } => (span, kind),
                _ => continue,
            };
            // Only consider explicit non-aromatic vs aromatic bonds
            let is_non_aromatic_bond = matches!(kind, BondKind::Single | BondKind::Double | BondKind::Triple | BondKind::Quadruple);
            let is_aromatic_bond = matches!(kind, BondKind::Aromatic);
            if !(is_non_aromatic_bond || is_aromatic_bond) { continue; }

            let left = (0..i).rfind(|&j| !matches!(segs[j], Segment::WhitespaceBlock { .. }));
            let right = ((i + 1)..segs.len()).find(|&j| !matches!(segs[j], Segment::WhitespaceBlock { .. }));
            let (Some(li), Some(ri)) = (left, right) else { continue };
            let left_arom = match &segs[li] {
                Segment::AtomSimple { raw, .. } => is_aromatic_simple(raw),
                Segment::AtomBracket { inner, .. } => is_aromatic_bracket(inner),
                _ => false,
            };
            let right_arom = match &segs[ri] {
                Segment::AtomSimple { raw, .. } => is_aromatic_simple(raw),
                Segment::AtomBracket { inner, .. } => is_aromatic_bracket(inner),
                _ => false,
            };

            if (is_non_aromatic_bond && left_arom && right_arom) || (is_aromatic_bond && !(left_arom && right_arom)) {
                emit.candidate(DiagnosticCandidate {
                    code: Code("LEX_INCONSISTENT_AROMATICITY"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span,
                    message: "Inconsistent aromaticity on bond",
                    scope: Scope::Global,
                });
            }
        }
    }
}
pub static AROM_INCONSISTENT_RULE: InconsistentAromaticityRule = InconsistentAromaticityRule;
