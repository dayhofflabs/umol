//! Registry for SMILES linting rules.

use std::collections::HashMap;

use super::context::LintContext;
use super::emitter::Emitter;
use super::rules::{Phase, Rule};
use crate::diagnostics::DiagnosticsReport;

pub struct RuleRegistry {
    rules_by_phase: HashMap<Phase, Vec<&'static dyn Rule>>,
}

impl Default for RuleRegistry {
    fn default() -> Self {
        Self::new()
    }
}

impl RuleRegistry {
    pub fn new() -> Self {
        Self {
            rules_by_phase: HashMap::new(),
        }
    }
    pub fn register(&mut self, rule: &'static dyn Rule) {
        self.rules_by_phase
            .entry(rule.phase())
            .or_default()
            .push(rule);
    }
}

pub struct LintEngine {
    registry: RuleRegistry,
}

impl LintEngine {
    pub fn new(registry: RuleRegistry) -> Self {
        Self { registry }
    }
    pub fn run(&self, ctx: &LintContext, report: &mut DiagnosticsReport) {
        // Phases in order
        let phases = [
            Phase::Lex,
            Phase::Bracket,
            Phase::RingStyle,
            Phase::Parse,
            Phase::Semantic,
        ];
        for ph in phases.iter() {
            if let Some(rules) = self.registry.rules_by_phase.get(ph) {
                let mut emitter = Emitter::new(report);
                for r in rules {
                    r.check(ctx, &mut emitter);
                }
                emitter.flush();
            }
        }
    }
}
