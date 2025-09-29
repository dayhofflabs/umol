//! Lexical rules for SMILES linting.

// Removed logos dependency

use super::{Rule, RuleMeta};
use crate::diagnostics::{Category, Code, Severity, Span};
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
    fn check(&self, _ctx: &LintContext, _emit: &mut Emitter) {}
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
        // Simple inter-token whitespace detector: any internal ASCII whitespace is flagged
        let bytes = ctx.input.as_bytes();
        let n = bytes.len();
        // Trim leading terminator whitespace
        let mut start = 0usize;
        while start < n && matches!(bytes[start], b' ' | b'\t' | b'\n' | b'\r') { start += 1; }
        // If all whitespace, ok
        if start == n { return; }
        // Trim trailing terminator whitespace
        let mut end = n;
        while end > start && matches!(bytes[end - 1], b' ' | b'\t' | b'\n' | b'\r') { end -= 1; }
        for i in start..end {
            let b = bytes[i];
            if matches!(b, b' ' | b'\t' | b'\n' | b'\r') {
                emit.candidate(DiagnosticCandidate {
                    code: Code("LEX_INTERTOKEN_WHITESPACE"),
                    category: Category::Lex,
                    severity: Severity::Error,
                    span: Span::new(i, i + 1),
                    message: "Inter-token whitespace is not allowed",
                    scope: Scope::Global,
                });
                break;
            }
        }
    }
}
pub static WS_RULE: WhitespaceRule = WhitespaceRule;

pub struct CommentsRule;
static META_COMMENTS: RuleMeta = RuleMeta {
    id: "LEX_COMMENTS",
    category: Category::Lex,
    default_severity: Severity::Error,
};
impl Rule for CommentsRule {
    fn meta(&self) -> &'static RuleMeta { &META_COMMENTS }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let bytes = ctx.input.as_bytes();
        let mut i = 0usize;
        while i + 1 < bytes.len() {
            let b0 = bytes[i];
            let b1 = bytes[i + 1];
            if b0 == b'/' && b1 == b'/' {
                let start = i;
                i += 2;
                while i < bytes.len() && bytes[i] != b'\n' && bytes[i] != b'\r' { i += 1; }
                emit.candidate(DiagnosticCandidate { code: Code("LEX_COMMENT"), category: Category::Lex, severity: Severity::Error, span: Span::new(start, i), message: "Comments are not allowed", scope: Scope::Global });
                continue;
            }
            if b0 == b'/' && b1 == b'*' {
                let start = i;
                i += 2;
                let mut closed = false;
                while i + 1 < bytes.len() {
                    if bytes[i] == b'*' && bytes[i + 1] == b'/' { i += 2; closed = true; break; }
                    i += 1;
                }
                if closed {
                    emit.candidate(DiagnosticCandidate { code: Code("LEX_COMMENT"), category: Category::Lex, severity: Severity::Error, span: Span::new(start, i), message: "Comments are not allowed", scope: Scope::Global });
                } else {
                    emit.candidate(DiagnosticCandidate { code: Code("LEX_UNTERMINATED_BLOCK_COMMENT"), category: Category::Lex, severity: Severity::Error, span: Span::new(start, ctx.input.len()), message: "Unterminated block comment", scope: Scope::Global });
                    break;
                }
                continue;
            }
            i += 1;
        }
    }
}
pub static COMMENTS_RULE: CommentsRule = CommentsRule;

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
    fn check(&self, _ctx: &LintContext, _emit: &mut Emitter) {}
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
    fn check(&self, _ctx: &LintContext, _emit: &mut Emitter) {}
}
pub static DOT_RULES: DotRules = DotRules;

// Leading bond
pub struct LeadingBondRule;
static META_LEX_LEAD_BOND: RuleMeta = RuleMeta { id: "LEX_LEADING_BOND_RULE", category: Category::Lex, default_severity: Severity::Error };
impl Rule for LeadingBondRule {
    fn meta(&self) -> &'static RuleMeta { &META_LEX_LEAD_BOND }
    fn check(&self, _ctx: &LintContext, _emit: &mut Emitter) {}
}
pub static LEADING_BOND_RULE: LeadingBondRule = LeadingBondRule;

// Leading ring
pub struct LeadingRingRule;
static META_LEX_LEAD_RING: RuleMeta = RuleMeta { id: "LEX_LEADING_RING_RULE", category: Category::Lex, default_severity: Severity::Error };
impl Rule for LeadingRingRule {
    fn meta(&self) -> &'static RuleMeta { &META_LEX_LEAD_RING }
    fn check(&self, _ctx: &LintContext, _emit: &mut Emitter) {}
}
pub static LEADING_RING_RULE: LeadingRingRule = LeadingRingRule;

// Consecutive bonds
pub struct ConsecutiveBondsRule;
static META_LEX_CONS_BONDS: RuleMeta = RuleMeta { id: "LEX_CONS_BONDS_RULE", category: Category::Lex, default_severity: Severity::Error };
impl Rule for ConsecutiveBondsRule {
    fn meta(&self) -> &'static RuleMeta { &META_LEX_CONS_BONDS }
    fn check(&self, _ctx: &LintContext, _emit: &mut Emitter) {}
}
pub static CONSECUTIVE_BONDS_RULE: ConsecutiveBondsRule = ConsecutiveBondsRule;

// Top-level group trailing
pub struct TopGroupTrailingRule;
static META_TOP_GRP: RuleMeta = RuleMeta { id: "LEX_TOP_GRP_TRAILING_RULE", category: Category::Lex, default_severity: Severity::Error };
impl Rule for TopGroupTrailingRule {
    fn meta(&self) -> &'static RuleMeta { &META_TOP_GRP }
    fn check(&self, _ctx: &LintContext, _emit: &mut Emitter) {}
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
    fn check(&self, _ctx: &LintContext, _emit: &mut Emitter) {}
}
pub static AROM_INCONSISTENT_RULE: InconsistentAromaticityRule = InconsistentAromaticityRule;
