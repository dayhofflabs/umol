//! Emitter for SMILES linting.

use std::mem;

use crate::diagnostics::{Category, Code, Diagnostic, DiagnosticsReport, Severity, Span};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Scope {
    Global,
    Bracket { start: usize, end: usize },
}

pub struct DiagnosticCandidate {
    pub code: Code,
    pub category: Category,
    pub severity: Severity,
    pub span: Span,
    pub message: &'static str,
    pub scope: Scope,
}

#[derive(Default)]
struct GateState {
    bracket_error_ranges: Vec<(usize, usize)>,
}

pub struct Emitter<'a> {
    pub report: &'a mut DiagnosticsReport,
    gates: GateState,
    staged: Vec<DiagnosticCandidate>,
}

impl<'a> Emitter<'a> {
    pub fn new(report: &'a mut DiagnosticsReport) -> Self {
        Self {
            report,
            gates: GateState::default(),
            staged: Vec::new(),
        }
    }

    pub fn candidate(&mut self, c: DiagnosticCandidate) {
        if matches!(c.scope, Scope::Bracket { .. }) && matches!(c.severity, Severity::Error) {
            self.gates
                .bracket_error_ranges
                .push((c.span.start, c.span.end));
        }
        self.staged.push(c);
    }

    pub fn flush(&mut self) {
        let staged = mem::take(&mut self.staged);
        for c in staged.into_iter() {
            if self.suppressed(&c) {
                continue;
            }
            self.report.push(Diagnostic {
                code: c.code,
                category: c.category,
                severity: c.severity,
                span: c.span,
                message: c.message.into(),
                details: None,
            });
        }
    }

    fn suppressed(&self, c: &DiagnosticCandidate) -> bool {
        // Suppress STYLE inside a bracket scope if any Error exists in same bracket range
        if c.category == Category::Style {
            if let Scope::Bracket { start, end } = c.scope {
                for (es, ee) in &self.gates.bracket_error_ranges {
                    if *es >= start && *ee <= end {
                        return true;
                    }
                }
            }
        }
        false
    }
}
