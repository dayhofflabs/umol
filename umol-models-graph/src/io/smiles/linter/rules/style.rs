//! Style rules for SMILES linting.

use super::{Phase, Rule, RuleMeta};
use crate::diagnostics::{Category, Severity};
use crate::io::smiles::linter::emitter::{DiagnosticCandidate, Emitter, Scope};
use crate::io::smiles::linter::{style as sh, LintContext};

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
    fn phase(&self) -> Phase {
        Phase::RingStyle
    }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let mut tmp = crate::diagnostics::DiagnosticsReport::new();
        sh::lint_style_percent_single_digit(ctx.input, &mut tmp);
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
    fn phase(&self) -> Phase {
        Phase::RingStyle
    }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let mut tmp = crate::diagnostics::DiagnosticsReport::new();
        sh::lint_style_bonds(ctx.input, &mut tmp);
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
pub static BOND_STYLE_RULE: BondStyleRule = BondStyleRule;
