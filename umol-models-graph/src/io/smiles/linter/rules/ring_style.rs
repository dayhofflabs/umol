//! Ring style rules for SMILES linting.

use super::{Phase, Rule, RuleMeta};
use crate::diagnostics::{Category, DiagnosticsReport, Severity};
use crate::io::smiles::linter::emitter::{DiagnosticCandidate, Emitter, Scope};
use crate::io::smiles::linter::utils::lint_ring_style;
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
    fn phase(&self) -> Phase {
        Phase::RingStyle
    }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let mut tmp = DiagnosticsReport::new();
        lint_ring_style(ctx.input, &mut tmp);
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
pub static RING_STYLE_RULE: RingStyleRule = RingStyleRule;
