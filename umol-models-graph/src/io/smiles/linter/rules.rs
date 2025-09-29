//! Rules for SMILES linting

use super::context::LintContext;
use super::emitter::Emitter;
use crate::diagnostics::{Category, Severity};

pub struct RuleMeta {
    pub id: &'static str,
    pub category: Category,
    pub default_severity: Severity,
}

pub trait Rule: Sync + Send {
    fn meta(&self) -> &'static RuleMeta;
    fn check(&self, ctx: &LintContext, emit: &mut Emitter);
}
