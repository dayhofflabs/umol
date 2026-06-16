//! Shared operator enums used across the AST.

/// Relational operators for `ValuePredicate::Rel`. Negation-closed: every
/// operator's logical negation is another operator (`Lt`↔`Ge`, `Le`↔`Gt`,
/// `Eq`↔`Ne`). Canonicalization orients `Gt`/`Ge` into `Lt`/`Le` (operand swap),
/// so a canonical relation uses only `Lt`/`Le`/`Eq`/`Ne`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RelOp {
    Le,
    Ge,
    Eq,
    Lt,
    Gt,
    Ne,
}

/// Membership operator: `In` / `NotIn`. Negation-closed (`In`↔`NotIn`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MemOp {
    In,
    NotIn,
}
