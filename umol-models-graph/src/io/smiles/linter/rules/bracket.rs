//! Bracket rules for SMILES linting.

use umol_data::isotope::Isotope as KnownIsotope;
use umol_data::Element;

use super::{Phase, Rule, RuleMeta};
use crate::diagnostics::{Category, Code, Severity, Span};
use crate::io::smiles::linter::emitter::{DiagnosticCandidate, Emitter, Scope};
use crate::io::smiles::linter::{bracket as bh, LintContext};

pub struct BracketRule;
static META_BRKT: RuleMeta = RuleMeta {
    id: "BRKT_RULE",
    category: Category::Brkt,
    default_severity: Severity::Error,
};
impl Rule for BracketRule {
    fn meta(&self) -> &'static RuleMeta {
        &META_BRKT
    }
    fn phase(&self) -> Phase {
        Phase::Bracket
    }
    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let bytes = ctx.input.as_bytes();
        let mut i = 0usize;
        while i < bytes.len() {
            if bytes[i] == b'[' {
                if let Some(close) = bh::find_closing_bracket(bytes, i + 1) {
                    let inner = &ctx.input[i + 1..close];
                    let parsed = bh::parse_bracket_inner(inner);
                    let scope = Scope::Bracket {
                        start: i,
                        end: close + 1,
                    };
                    let mut had_error = false;
                    if matches!(parsed.element, Some(Element::H)) && parsed.hcount.is_some() {
                        emit.candidate(DiagnosticCandidate {
                            code: Code("BRKT_H_ON_H"),
                            category: Category::Brkt,
                            severity: Severity::Error,
                            span: Span::new(i, close + 1),
                            message: "Hydrogen element must not have an H-count",
                            scope,
                        });
                        had_error = true;
                    }
                    if let (Some(elem), Some(h)) = (parsed.element, parsed.hcount) {
                        if h as u8 > elem.max_implicit_hydrogens() {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("NUM_HCOUNT_EXCEEDS_MAX_IMPLICIT"),
                                category: Category::Num,
                                severity: Severity::Warning,
                                span: Span::new(i, close + 1),
                                message: "H-count exceeds element's max implicit hydrogens",
                                scope,
                            });
                        }
                    }
                    if let Some(class) = parsed.class {
                        if class > 9999 {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("NUM_CLASS_TOO_LARGE"),
                                category: Category::Num,
                                severity: Severity::Error,
                                span: Span::new(i, close + 1),
                                message: "Atom class must be <= 9999",
                                scope,
                            });
                            had_error = true;
                        }
                    }
                    if let Some(q) = parsed.charge {
                        if q.unsigned_abs() > 15 {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("NUM_CHARGE_OUT_OF_RANGE"),
                                category: Category::Num,
                                severity: Severity::Error,
                                span: Span::new(i, close + 1),
                                message: "Absolute charge must be <= 15",
                                scope,
                            });
                            had_error = true;
                        }
                        if let Some(elem) = parsed.element {
                            let (min_q, max_q) = elem.charge_bounds();
                            if q < min_q as i32 || q > max_q as i32 {
                                emit.candidate(DiagnosticCandidate {
                                    code: Code("NUM_CHARGE_OUTSIDE_ELEMENT_RANGE"),
                                    category: Category::Num,
                                    severity: Severity::Warning,
                                    span: Span::new(i, close + 1),
                                    message: "Charge outside element-supported bounds",
                                    scope,
                                });
                            }
                            if q > 0 && (q as u8) > elem.valence_electrons() {
                                emit.candidate(DiagnosticCandidate {
                                    code: Code("NUM_CHARGE_EXCEEDS_VALENCE_ELECTRONS"),
                                    category: Category::Num,
                                    severity: Severity::Warning,
                                    span: Span::new(i, close + 1),
                                    message: "Positive charge exceeds valence electrons",
                                    scope,
                                });
                            }
                        }
                    }
                    if let Some(isotope) = parsed.isotope {
                        if isotope > 999 {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("NUM_ISOTOPE_TOO_LARGE"),
                                category: Category::Num,
                                severity: Severity::Error,
                                span: Span::new(i, close + 1),
                                message: "Isotope mass number must be <= 999",
                                scope,
                            });
                            had_error = true;
                        } else if isotope > 0 {
                            if let Some(elem) = parsed.element {
                                if !KnownIsotope::is_catalogued(elem, isotope) {
                                    emit.candidate(DiagnosticCandidate {
                                        code: Code("NUM_ISOTOPE_UNCATALOGUED"),
                                        category: Category::Num,
                                        severity: Severity::Warning,
                                        span: Span::new(i, close + 1),
                                        message: "Isotope is not catalogued",
                                        scope,
                                    });
                                }
                            }
                        }
                    }
                    if !had_error && bh::is_bare_organic(inner) {
                        emit.candidate(DiagnosticCandidate {
                            code: Code("STYLE_BRACKET_ORGANIC"),
                            category: Category::Style,
                            severity: Severity::Warning,
                            span: Span::new(i, close + 1),
                            message: "Prefer bare organic atom over bracketed form",
                            scope,
                        });
                    }
                    if !had_error && bh::inner_contains_h1(inner) {
                        if let Some((h_start, h_end)) = bh::find_subslice(inner, "H1") {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("STYLE_HCOUNT_ONE_SIMPLE"),
                                category: Category::Style,
                                severity: Severity::Warning,
                                span: Span::new(i + 1 + h_start, i + 1 + h_end),
                                message: "Prefer 'H' over 'H1'",
                                scope,
                            });
                        } else {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("STYLE_HCOUNT_ONE_SIMPLE"),
                                category: Category::Style,
                                severity: Severity::Warning,
                                span: Span::new(i, close + 1),
                                message: "Prefer 'H' over 'H1'",
                                scope,
                            });
                        }
                    }
                    // STYLE_CHARGE_SIGN_SIMPLE: prefer [+]/[-] over [+1]/[-1]
                    if !had_error {
                        if let Some((c_start, c_end)) = bh::find_charge_plus_minus_one(inner) {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("STYLE_CHARGE_SIGN_SIMPLE"),
                                category: Category::Style,
                                severity: Severity::Warning,
                                span: Span::new(i + 1 + c_start, i + 1 + c_end),
                                message: "Prefer [+]/[-] over [+1]/[-1]",
                                scope,
                            });
                        }
                    }
                    if let Some((h2s, h2e)) = bh::find_h_two_digits(inner) {
                        emit.candidate(DiagnosticCandidate {
                            code: Code("BRKT_HCOUNT_TWO_DIGITS"),
                            category: Category::Brkt,
                            severity: Severity::Error,
                            span: Span::new(i + 1 + h2s, i + 1 + h2e),
                            message: "Hydrogen count must be a single digit",
                            scope,
                        });
                        had_error = true;
                    }
                    if let Some((cs, ce, neg)) = bh::find_class_issues(inner) {
                        if neg {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("NUM_CLASS_NEGATIVE"),
                                category: Category::Num,
                                severity: Severity::Error,
                                span: Span::new(i + 1 + cs, i + 1 + ce),
                                message: "Atom class must be non-negative",
                                scope,
                            });
                        } else {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("BRKT_EMPTY_CLASS"),
                                category: Category::Brkt,
                                severity: Severity::Error,
                                span: Span::new(i + 1 + cs, i + 1 + ce),
                                message: "Class field ':' must be followed by digits",
                                scope,
                            });
                        }
                        had_error = true;
                    }
                    if !had_error && bh::bracket_order_misordered(inner) {
                        emit.candidate(DiagnosticCandidate {
                            code: Code("STYLE_BRKT_ORDER"),
                            category: Category::Style,
                            severity: Severity::Warning,
                            span: Span::new(i, close + 1),
                            message: "Prefer [chirality][H][charge][class] ordering",
                            scope,
                        });
                    }
                    i = close + 1;
                    continue;
                }
            }
            i += 1;
        }
    }
}
pub static BRACKET_RULE: BracketRule = BracketRule;
