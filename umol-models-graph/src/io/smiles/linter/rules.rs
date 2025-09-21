//! Rules for SMILES linting

use super::context::LintContext;
use super::emitter::Emitter;
use crate::diagnostics::{Category, Severity};

pub mod bracket;
pub mod branch;
pub mod lexical;
pub mod ring;
pub mod style;

pub use bracket::BRACKET_RULE;
pub use branch::BRANCH_RULE;
pub use lexical::{DOT_RULES, LEX_ERRORS_RULE, TRAILING_BOND_RULE, WS_RULE};
pub use ring::{RING_ERRORS_RULE, RING_STYLE_RULE};
pub use style::{BOND_STYLE_RULE, STYLE_PCT_RULE};

pub struct RuleMeta {
    pub id: &'static str,
    pub category: Category,
    pub default_severity: Severity,
}

pub trait Rule: Sync + Send {
    fn meta(&self) -> &'static RuleMeta;
    fn check(&self, ctx: &LintContext, emit: &mut Emitter);
}
