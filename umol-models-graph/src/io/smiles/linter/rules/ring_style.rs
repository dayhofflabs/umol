//! Ring style rules for SMILES linting.

use super::{Phase, Rule, RuleMeta};
use crate::diagnostics::{Category, Severity, Span};
use crate::io::smiles::linter::emitter::{DiagnosticCandidate, Emitter, Scope};
use crate::io::smiles::linter::LintContext;
use crate::io::smiles::segment::Segment;
use indexmap::IndexMap;

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
    fn phase(&self) -> Phase {
        Phase::RingStyle
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
                code: crate::diagnostics::Code("STYLE_FIRST_RING_NOT_ONE"),
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
                        code: crate::diagnostics::Code("STYLE_NONCONSECUTIVE_RING_NUMBERING"),
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
                    code: crate::diagnostics::Code("STYLE_REUSED_RING_INDICES"),
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
