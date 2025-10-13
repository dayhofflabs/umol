//! Emitter for SMILES linting (simple config-driven gating).

use std::collections::{HashMap, HashSet};

use super::super::diagnostics::{Diagnostic, DiagnosticCode, DiagnosticList, Severity};

#[derive(Default)]
pub struct LintConfig {
    pub enabled_codes: Option<HashSet<&'static str>>, // if Some, only these codes are emitted
    pub disabled_codes: HashSet<&'static str>,        // always suppressed codes
    pub severity_overrides: HashMap<&'static str, Severity>, // per-code severity override
}

pub struct Emitter<'a> {
    pub report: &'a mut DiagnosticList,
    pub config: LintConfig,
}

impl<'a> Emitter<'a> {
    pub fn new(report: &'a mut DiagnosticList) -> Self {
        Self {
            report,
            config: LintConfig::default(),
        }
    }

    pub fn with_config(report: &'a mut DiagnosticList, config: LintConfig) -> Self {
        Self { report, config }
    }

    pub fn emit(&mut self, mut d: Diagnostic) {
        let code_str = d.code.as_str();
        if self.config.disabled_codes.contains(code_str) {
            return;
        }
        if let Some(enabled) = &self.config.enabled_codes {
            if !enabled.contains(code_str) {
                return;
            }
        }
        if let Some(sev) = self.config.severity_overrides.get(code_str) {
            d.severity = *sev;
        }
        self.report.push(d);
    }
}
