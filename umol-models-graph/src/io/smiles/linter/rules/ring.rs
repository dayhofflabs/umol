//! Ring rules for SMILES linting: style and errors in one module.

use std::collections::{HashMap, HashSet};

use indexmap::IndexMap;

use super::{Rule, RuleMeta};
use crate::diagnostics::{Category, Code, Severity, Span};
use crate::io::smiles::iterators::{BondKind, Segment};
use crate::io::smiles::linter::emitter::{DiagnosticCandidate, Emitter, Scope};
use crate::io::smiles::linter::LintContext;

pub struct RingStyleRule;

static META_RING_STYLE: RuleMeta = RuleMeta {
    id: "STYLE_RINGS",
    category: Category::Style,
    default_severity: Severity::Warning,
};

impl Rule for RingStyleRule {
    fn meta(&self) -> &'static RuleMeta {
        &META_RING_STYLE
    }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let segs = ctx.segments();
        let mut used: Vec<u32> = Vec::new();
        let mut counts: IndexMap<u32, u32> = IndexMap::new();
        for seg in segs.iter() {
            if let Segment::RingClosure { index, .. } = seg {
                used.push(*index);
                *counts.entry(*index).or_insert(0) += 1;
            }
        }
        if used.is_empty() {
            return;
        }
        let mut set = used.clone();
        set.sort_unstable();
        set.dedup();
        if set[0] != 1 {
            emit.candidate(DiagnosticCandidate {
                code: Code("STYLE_FIRST_RING_NOT_ONE"),
                category: Category::Style,
                severity: Severity::Warning,
                span: Span::new(0, 0),
                message: "Prefer starting ring numbering at 1",
                scope: Scope::Global,
            });
        }
        if set.len() >= 2 {
            let mut prev = set[0];
            for &v in &set[1..] {
                if v > prev + 1 {
                    emit.candidate(DiagnosticCandidate {
                        code: Code("STYLE_NONCONSECUTIVE_RING_NUMBERING"),
                        category: Category::Style,
                        severity: Severity::Warning,
                        span: Span::new(0, 0),
                        message: "Prefer consecutive ring numbering",
                        scope: Scope::Global,
                    });
                    break;
                }
                prev = v;
            }
        }
        for (_k, c) in counts.iter() {
            if *c > 2 {
                emit.candidate(DiagnosticCandidate {
                    code: Code("STYLE_REUSED_RING_INDICES"),
                    category: Category::Style,
                    severity: Severity::Warning,
                    span: Span::new(0, 0),
                    message: "Avoid reusing the same ring number",
                    scope: Scope::Global,
                });
                break;
            }
        }
    }
}
pub static RING_STYLE_RULE: RingStyleRule = RingStyleRule;

pub struct RingErrorsRule;
static META_RING_ERR: RuleMeta = RuleMeta {
    id: "RING_ERRORS",
    category: Category::Ring,
    default_severity: Severity::Error,
};
impl Rule for RingErrorsRule {
    fn meta(&self) -> &'static RuleMeta {
        &META_RING_ERR
    }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let segs = ctx.segments();
        let mut unmatched: HashMap<u32, Span> = HashMap::new();
        let mut dirs: HashMap<u32, HashSet<BondKind>> = HashMap::new();
        let mut prev_non_ws_is_dir: Option<BondKind> = None;

        let flush_component = |emit: &mut Emitter, unmatched: &mut HashMap<u32, Span>| {
            for (_idx, span) in unmatched.drain() {
                emit.candidate(DiagnosticCandidate {
                    code: Code("RING_UNCLOSED"),
                    category: Category::Ring,
                    severity: Severity::Error,
                    span,
                    message: "Ring index not closed within component",
                    scope: Scope::Global,
                });
            }
        };

        for seg in segs.iter() {
            match seg {
                Segment::WhitespaceBlock { .. } => {}
                Segment::ComponentSeparator { .. } => {
                    flush_component(emit, &mut unmatched);
                    prev_non_ws_is_dir = None;
                }
                Segment::Bond { kind, .. } => {
                    prev_non_ws_is_dir = match kind {
                        BondKind::Up | BondKind::Down => Some(*kind),
                        _ => None,
                    };
                }
                Segment::RingClosure { index, span } => {
                    if unmatched.remove(index).is_none() {
                        unmatched.insert(*index, *span);
                    }
                    if let Some(dir) = prev_non_ws_is_dir {
                        let set = dirs.entry(*index).or_default();
                        set.insert(dir);
                        if set.contains(&BondKind::Up) && set.contains(&BondKind::Down) {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("RING_CONFLICT_DIR"),
                                category: Category::Ring,
                                severity: Severity::Error,
                                span: *span,
                                message: "Conflicting '/' and '\\' directions for ring index",
                                scope: Scope::Global,
                            });
                        }
                    }
                    prev_non_ws_is_dir = None;
                }
                _ => {
                    prev_non_ws_is_dir = None;
                }
            }
        }
        flush_component(emit, &mut unmatched);
    }
}
pub static RING_ERRORS_RULE: RingErrorsRule = RingErrorsRule;
