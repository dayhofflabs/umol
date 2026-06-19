//! Ring-membership scope and single-entry fact, shared by atom/bond/dative constraints.

use super::super::value::ValueAst;

/// `All` = total ring count; `Size(s)` = count of size-`s` rings. `All` sorts first.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RingScope {
    All,
    Size(u8),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RingMembershipAst {
    pub scope: RingScope,
    pub count: ValueAst,
}

impl RingMembershipAst {
    pub fn new(scope: RingScope, count: impl Into<ValueAst>) -> Self {
        Self {
            scope,
            count: count.into(),
        }
    }
}
