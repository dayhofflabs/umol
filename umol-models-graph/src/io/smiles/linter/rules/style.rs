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
