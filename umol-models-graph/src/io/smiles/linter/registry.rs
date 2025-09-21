//! Registry for SMILES linting rules.

use super::context::LintContext;
use super::emitter::Emitter;
use super::rules::Rule;
use crate::diagnostics::DiagnosticsReport;

pub struct RuleRegistry {
    rules_in_order: Vec<&'static dyn Rule>,
}

impl Default for RuleRegistry {
    fn default() -> Self {
        Self::new()
    }
}

impl RuleRegistry {
    pub fn new() -> Self {
        Self {
            rules_in_order: Vec::new(),
        }
    }
    pub fn register(&mut self, rule: &'static dyn Rule) {
        self.rules_in_order.push(rule);
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
        let mut emitter = Emitter::new(report);
        for r in &self.registry.rules_in_order {
            r.check(ctx, &mut emitter);
        }
        emitter.flush();
    }
}
