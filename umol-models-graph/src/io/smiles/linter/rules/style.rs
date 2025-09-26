//! Style rules for SMILES linting.

use regex::Regex;

use super::{Rule, RuleMeta};
use crate::diagnostics::{Category, Severity};
use crate::io::smiles::iterators::{BondKind, Segment};
use crate::io::smiles::linter::emitter::{DiagnosticCandidate, Emitter, Scope};
use crate::io::smiles::linter::LintContext;

pub struct StylePercentSingleDigitRule;

static META_PCT: RuleMeta = RuleMeta {
    id: "STYLE_PCT",
    category: Category::Style,
    default_severity: Severity::Warning,
};

impl Rule for StylePercentSingleDigitRule {
    fn meta(&self) -> &'static RuleMeta {
        &META_PCT
    }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let re = Regex::new(r"%(0[1-9])").unwrap();
        for m in re.find_iter(ctx.input) {
            let start = m.start();
            emit.candidate(DiagnosticCandidate {
                code: crate::diagnostics::Code("STYLE_UNNECESSARY_PERCENT_RING_INDEX"),
                category: Category::Style,
                severity: Severity::Warning,
                span: crate::diagnostics::Span::new(start, start + 3),
                message: "Prefer single-digit ring number for 1..9",
                scope: Scope::Global,
            });
        }
    }
}
pub static STYLE_PCT_RULE: StylePercentSingleDigitRule = StylePercentSingleDigitRule;

pub struct BondStyleRule;

static META_BOND_STYLE: RuleMeta = RuleMeta {
    id: "STYLE_BONDS",
    category: Category::Style,
    default_severity: Severity::Warning,
};

impl Rule for BondStyleRule {
    fn meta(&self) -> &'static RuleMeta {
        &META_BOND_STYLE
    }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let segs = ctx.segments();
        let is_arom_str =
            |raw: &str| matches!(raw, "b" | "c" | "n" | "o" | "p" | "s" | "se" | "as");
        for i in 0..segs.len() {
            match segs[i] {
                Segment::Bond {
                    span,
                    kind: BondKind::Aromatic,
                } => {
                    let prev_arom = (0..i)
                        .rfind(|&j| !matches!(segs[j], Segment::WhitespaceBlock { .. }))
                        .and_then(|j| match &segs[j] {
                            Segment::AtomSimple { raw, .. } => Some(is_arom_str(raw)),
                            _ => Some(false),
                        })
                        .unwrap_or(false);
                    let next_arom = ((i + 1)..segs.len())
                        .find(|&j| !matches!(segs[j], Segment::WhitespaceBlock { .. }))
                        .and_then(|j| match &segs[j] {
                            Segment::AtomSimple { raw, .. } => Some(is_arom_str(raw)),
                            _ => Some(false),
                        })
                        .unwrap_or(false);
                    if prev_arom && next_arom {
                        emit.candidate(DiagnosticCandidate {
                            code: crate::diagnostics::Code("STYLE_EXPLICIT_AROMATIC_BOND"),
                            category: Category::Style,
                            severity: Severity::Warning,
                            span,
                            message: "Avoid explicit ':' between aromatic atoms",
                            scope: Scope::Global,
                        });
                    }
                }
                Segment::Bond {
                    span,
                    kind: BondKind::Single,
                } => {
                    let prev_arom = (0..i)
                        .rfind(|&j| !matches!(segs[j], Segment::WhitespaceBlock { .. }))
                        .and_then(|j| match &segs[j] {
                            Segment::AtomSimple { raw, .. } => Some(is_arom_str(raw)),
                            _ => Some(false),
                        })
                        .unwrap_or(false);
                    let next_arom = ((i + 1)..segs.len())
                        .find(|&j| !matches!(segs[j], Segment::WhitespaceBlock { .. }))
                        .and_then(|j| match &segs[j] {
                            Segment::AtomSimple { raw, .. } => Some(is_arom_str(raw)),
                            _ => Some(false),
                        })
                        .unwrap_or(false);
                    if !(prev_arom && next_arom) {
                        emit.candidate(DiagnosticCandidate {
                            code: crate::diagnostics::Code("STYLE_EXPLICIT_SINGLE_BOND"),
                            category: Category::Style,
                            severity: Severity::Warning,
                            span,
                            message: "Avoid explicit '-' when default applies",
                            scope: Scope::Global,
                        });
                    }
                }
                _ => {}
            }
        }
    }
}
pub static BOND_STYLE_RULE: BondStyleRule = BondStyleRule;

pub struct GroupParenStyleRule;

static META_GRP_PARENS: RuleMeta = RuleMeta {
    id: "STYLE_GROUP_PARENS",
    category: Category::Style,
    default_severity: Severity::Warning,
};

impl Rule for GroupParenStyleRule {
    fn meta(&self) -> &'static RuleMeta { &META_GRP_PARENS }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let segs = ctx.segments();
        // Detect top-level redundant parens: (chain) where removing outer parens does not change tokens
        // Heuristic: if entire input is exactly one group and group contents contain no top-level-only constructs
        if segs.len() >= 2 {
            let mut non_ws = segs.iter().filter(|s| !matches!(s, Segment::WhitespaceBlock { .. }));
            let first = non_ws.next();
            let last = segs.iter().rfind(|s| !matches!(s, Segment::WhitespaceBlock { .. }));
            if matches!(first, Some(Segment::BranchOpen { .. })) && matches!(last, Some(Segment::BranchClose { .. })) {
                // Treat as top-level group; warn unnecessary top-level parens
                if let (Some(Segment::BranchOpen { span: open }), Some(Segment::BranchClose { span: close })) = (first, last) {
                    emit.candidate(DiagnosticCandidate {
                        code: crate::diagnostics::Code("STYLE_UNNECESSARY_TOPLEVEL_PARENS"),
                        category: Category::Style,
                        severity: Severity::Warning,
                        span: crate::diagnostics::Span::new(open.start, close.end),
                        message: "Unnecessary top-level grouping parentheses",
                        scope: Scope::Global,
                    });
                }
            }
        }

        // Detect immediately nested redundant parens: ((...))
        let mut i = 0usize;
        while i + 3 < segs.len() {
            match (&segs[i], &segs[i + 1], &segs[i + 2], &segs[i + 3]) {
                (
                    Segment::BranchOpen { span: s1 },
                    Segment::BranchOpen { .. },
                    ..,
                    Segment::BranchClose { span: s4 },
                ) => {
                    // Skip scanning inner; this is a simple nested-open-close pattern
                    emit.candidate(DiagnosticCandidate {
                        code: crate::diagnostics::Code("STYLE_REDUNDANT_NESTED_PARENS"),
                        category: Category::Style,
                        severity: Severity::Warning,
                        span: crate::diagnostics::Span::new(s1.start, s4.end),
                        message: "Redundant nested grouping parentheses",
                        scope: Scope::Global,
                    });
                    i += 1;
                }
                _ => i += 1,
            }
        }
    }
}
pub static GROUP_PARENS_RULE: GroupParenStyleRule = GroupParenStyleRule;

pub struct MixedAfterAtomStyleRule;

static META_MIXED_AFTER: RuleMeta = RuleMeta {
    id: "STYLE_MIXED_AFTER",
    category: Category::Style,
    default_severity: Severity::Warning,
};

impl Rule for MixedAfterAtomStyleRule {
    fn meta(&self) -> &'static RuleMeta { &META_MIXED_AFTER }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let segs = ctx.segments();
        // Look for an atom followed by a sequence that mixes branches and ring closers interleaved
        // Heuristic: after AtomSimple/AtomBracket, if we see both a BranchOpen and a RingClosure before next Atom or component end, warn
        let mut i = 0usize;
        while i < segs.len() {
            let mut saw_atom = false;
            if matches!(segs[i], Segment::AtomSimple { .. } | Segment::AtomBracket { .. }) {
                saw_atom = true;
            }
            if saw_atom {
                let mut j = i + 1;
                let mut saw_branch = false;
                let mut saw_ring = false;
                while j < segs.len() {
                    match &segs[j] {
                        Segment::WhitespaceBlock { .. } => {}
                        Segment::RingClosure { .. } => { saw_ring = true; }
                        Segment::BranchOpen { .. } => { saw_branch = true; }
                        Segment::NewComponent { .. } | Segment::AtomSimple { .. } | Segment::AtomBracket { .. } => { break; }
                        Segment::BranchClose { .. } | Segment::Bond { .. } => {}
                        _ => {}
                    }
                    j += 1;
                }
                if saw_branch && saw_ring {
                    // Find span from first after-atom token to just before j
                    let start = match &segs[i + 1] { Segment::WhitespaceBlock { .. } => ctx.input.len().min(ctx.input.len()), _ => match &segs[i + 1] { Segment::BranchOpen { span } | Segment::RingClosure { span, .. } | Segment::Bond { span, .. } | Segment::NewComponent { span } => span.start, _ => i }, };
                    let end = match j.checked_sub(1).and_then(|k| match &segs[k] { Segment::BranchOpen { span } | Segment::BranchClose { span } | Segment::RingClosure { span, .. } | Segment::Bond { span, .. } | Segment::NewComponent { span } => Some(span.end), _ => None }) { Some(e) => e, None => ctx.input.len() };
                    emit.candidate(DiagnosticCandidate {
                        code: crate::diagnostics::Code("STYLE_MIXED_RINGBONDS_BRANCHES"),
                        category: Category::Style,
                        severity: Severity::Warning,
                        span: crate::diagnostics::Span::new(start, end),
                        message: "Prefer not to mix branches and ringbonds after the same atom",
                        scope: Scope::Global,
                    });
                }
                i = j;
                continue;
            }
            i += 1;
        }
    }
}
pub static MIXED_AFTER_RULE: MixedAfterAtomStyleRule = MixedAfterAtomStyleRule;

pub struct RingNumberingStyleRule;

static META_RING_NUM: RuleMeta = RuleMeta {
    id: "STYLE_RING_NUM",
    category: Category::Style,
    default_severity: Severity::Warning,
};

impl Rule for RingNumberingStyleRule {
    fn meta(&self) -> &'static RuleMeta { &META_RING_NUM }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let segs = ctx.segments();
        // STYLE_FIRST_RING_NOT_ONE / STYLE_NONCONSECUTIVE_RING_NUMBERING / STYLE_REUSED_RING_INDICES
        // Track first-seen order of ring indices per component
        let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
        let mut order: Vec<u32> = Vec::new();
        for i in 0..segs.len() {
            match segs[i] {
                Segment::RingClosure { span: _, index } => {
                    if !seen.contains(&index) {
                        seen.insert(index);
                        order.push(index);
                    }
                }
                Segment::NewComponent { .. } => {
                    // Flush order-based lints for the component
                    if !order.is_empty() {
                        // First ring not one (only warn if 1 is used anywhere but not first)
                        if order[0] != 1 && order.iter().any(|&v| v == 1) {
                            emit.candidate(DiagnosticCandidate {
                                code: crate::diagnostics::Code("STYLE_FIRST_RING_NOT_ONE"),
                                category: Category::Style,
                                severity: Severity::Warning,
                                span: crate::diagnostics::Span::new(0, 0),
                                message: "Prefer starting ring numbering at 1",
                                scope: Scope::Global,
                            });
                        }
                        // Non-consecutive numbering (consider numeric value, 0 allowed)
                        let mut prev = order[0];
                        for &v in &order[1..] {
                            if v != prev && v != prev + 1 {
                                emit.candidate(DiagnosticCandidate {
                                    code: crate::diagnostics::Code("STYLE_NONCONSECUTIVE_RING_NUMBERING"),
                                    category: Category::Style,
                                    severity: Severity::Warning,
                                    span: crate::diagnostics::Span::new(0, 0),
                                    message: "Avoid non-consecutive ring numbers in the parsing sequence",
                                    scope: Scope::Global,
                                });
                                break;
                            }
                            prev = v;
                        }
                    }
                    seen.clear();
                    order.clear();
                }
                _ => {}
            }
        }
        // Flush at end of input
        if !order.is_empty() {
            if order[0] != 1 && order.iter().any(|&v| v == 1) {
                emit.candidate(DiagnosticCandidate {
                    code: crate::diagnostics::Code("STYLE_FIRST_RING_NOT_ONE"),
                    category: Category::Style,
                    severity: Severity::Warning,
                    span: crate::diagnostics::Span::new(0, 0),
                    message: "Prefer starting ring numbering at 1",
                    scope: Scope::Global,
                });
            }
            let mut prev = order[0];
            for &v in &order[1..] {
                if v != prev && v != prev + 1 {
                    emit.candidate(DiagnosticCandidate {
                        code: crate::diagnostics::Code("STYLE_NONCONSECUTIVE_RING_NUMBERING"),
                        category: Category::Style,
                        severity: Severity::Warning,
                        span: crate::diagnostics::Span::new(0, 0),
                        message: "Avoid non-consecutive ring numbers in the parsing sequence",
                        scope: Scope::Global,
                    });
                    break;
                }
                prev = v;
            }
        }
        // STYLE_REUSED_RING_INDICES
        // Flag if the same index appears more than twice in a component (open+close is typical once, but we count appearances)
        let mut counts: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
        for i in 0..segs.len() {
            match segs[i] {
                Segment::RingClosure { index, .. } => {
                    *counts.entry(index).or_default() += 1;
                }
                Segment::NewComponent { .. } => {
                    if counts.values().any(|&v| v > 2) {
                        emit.candidate(DiagnosticCandidate {
                            code: crate::diagnostics::Code("STYLE_REUSED_RING_INDICES"),
                            category: Category::Style,
                            severity: Severity::Warning,
                            span: crate::diagnostics::Span::new(0, 0),
                            message: "Avoid reusing the same ring digit within a component",
                            scope: Scope::Global,
                        });
                    }
                    counts.clear();
                }
                _ => {}
            }
        }
        if counts.values().any(|&v| v > 2) {
            emit.candidate(DiagnosticCandidate {
                code: crate::diagnostics::Code("STYLE_REUSED_RING_INDICES"),
                category: Category::Style,
                severity: Severity::Warning,
                span: crate::diagnostics::Span::new(0, 0),
                message: "Avoid reusing the same ring digit within a component",
                scope: Scope::Global,
            });
        }
    }
}
pub static RING_NUMBERING_RULE: RingNumberingStyleRule = RingNumberingStyleRule;
