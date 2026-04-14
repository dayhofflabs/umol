//! Utility functions for converting between AST and IR types.

use umol_shared::spin::SpinMultiplicity;
use umol_shared::value_ast::ValueAst;

use super::atom_pattern::Pattern;
use crate::ast::config::{MultiplicityMode, NumericMode, UnpairedElectronsMode};
use crate::ast::error::LoweringError;

pub(crate) fn raise_u8_ground(value: u8, mode: &NumericMode) -> Option<ValueAst> {
    match (value, mode) {
        (0, NumericMode::Zero) => None,
        (n, _) => Some(ValueAst::Lit(n as i64)),
    }
}

pub(crate) fn raise_i8_ground(value: i8, mode: &NumericMode) -> Option<ValueAst> {
    match (value, mode) {
        (0, NumericMode::Zero) => None,
        (n, _) => Some(ValueAst::Lit(n as i64)),
    }
}

pub(crate) fn raise_u8_pattern(pattern: Pattern<u8>, mode: &NumericMode) -> Option<ValueAst> {
    match (pattern, mode) {
        (Pattern::Is(0), NumericMode::Zero) => None,
        (Pattern::Any, NumericMode::Required) => None,
        (Pattern::Any, NumericMode::Zero) => Some(ValueAst::Undetermined),
        (Pattern::Is(n), _) => Some(ValueAst::Lit(n as i64)),
    }
}

pub(crate) fn raise_i8_pattern(pattern: Pattern<i8>, mode: &NumericMode) -> Option<ValueAst> {
    match (pattern, mode) {
        (Pattern::Is(0), NumericMode::Zero) => None,
        (Pattern::Any, NumericMode::Required) => None,
        (Pattern::Any, NumericMode::Zero) => Some(ValueAst::Undetermined),
        (Pattern::Is(n), _) => Some(ValueAst::Lit(n as i64)),
    }
}

/// Coupled spin-state lowering for unpaired electrons (u) and multiplicity (m).
///
/// Resolution is sequential: u resolves first, m second. Each undetermined field
/// falls back to its mode:
///
///   u_ast Undetermined  → Zero: Is(0), Required: Any, Derived: derive from m_ast
///   u_ast Lit(n)        → Is(n)
///
///   m_ast Undetermined  → Required: Any, Derived: derive from resolved u
///   m_ast Lit(n)        → Is(SpinMultiplicity(n))
///
/// The asymmetry matters: u in Derived mode reads raw m_ast, but m in Derived
/// mode reads the already-resolved u pattern. This lets `C#s3` (m=3, u absent)
/// derive u=2, then confirm m=3 matches u+1.
pub(crate) fn lower_spin(
    u_ast: ValueAst,
    m_ast: ValueAst,
    u_mode: &UnpairedElectronsMode,
    m_mode: &MultiplicityMode,
) -> Result<(Pattern<u8>, Pattern<SpinMultiplicity>), LoweringError> {
    // Step 1: resolve unpaired electrons
    let unpaired_electrons = match &u_ast {
        ValueAst::Undetermined => match u_mode {
            UnpairedElectronsMode::Zero => Pattern::Is(0),
            UnpairedElectronsMode::Required => Pattern::Any,
            UnpairedElectronsMode::Derived => match &m_ast {
                ValueAst::Lit(m) => {
                    let m_val = u8::try_from(*m).map_err(|_| LoweringError::NonGround {
                        field: "multiplicity",
                    })?;
                    if m_val < 1 {
                        return Err(LoweringError::InvalidMultiplicity(0));
                    }
                    Pattern::Is(m_val - 1)
                }
                _ => Pattern::Any,
            },
        },
        ValueAst::Lit(n) => {
            Pattern::Is(u8::try_from(*n).map_err(|_| LoweringError::NonGround {
                field: "unpaired_electrons",
            })?)
        }
        _ => {
            return Err(LoweringError::NonGround {
                field: "unpaired_electrons",
            })
        }
    };

    // Step 2: resolve multiplicity (derives from resolved u when mode is Derived)
    let multiplicity = match &m_ast {
        ValueAst::Undetermined => match m_mode {
            MultiplicityMode::Required => Pattern::Any,
            MultiplicityMode::Derived => match unpaired_electrons {
                Pattern::Is(u) => {
                    let m_val = u
                        .checked_add(1)
                        .ok_or(LoweringError::InvalidMultiplicity(u8::MAX))?;
                    Pattern::Is(
                        SpinMultiplicity::from_multiplicity(m_val)
                            .ok_or(LoweringError::InvalidMultiplicity(m_val))?,
                    )
                }
                Pattern::Any => Pattern::Any,
            },
        },
        ValueAst::Lit(n) => {
            let m = u8::try_from(*n).map_err(|_| LoweringError::NonGround {
                field: "multiplicity",
            })?;
            Pattern::Is(
                SpinMultiplicity::from_multiplicity(m)
                    .ok_or(LoweringError::InvalidMultiplicity(m))?,
            )
        }
        _ => {
            return Err(LoweringError::NonGround {
                field: "multiplicity",
            })
        }
    };

    Ok((unpaired_electrons, multiplicity))
}

/// Coupled spin-state raising for ground types.
///
/// Returns `(u_ast, m_ast)` — the minimal AST fields that roundtrip through `lower_spin`.
///
/// `derived` = m == u + 1 (high-spin / Hund's rule relationship).
///
/// Suppression rules (Undetermined = field can be elided):
///   u side:
///     Zero    + u==0                        → Undetermined
///     Derived + derived + Required           → Undetermined
///     Derived + derived + Derived + u==0     → Undetermined
///     Derived + derived + Derived + u!=0     → Lit(u)
///     otherwise                              → Lit(u)
///   m side:
///     Derived + derived → Undetermined
///     otherwise         → Lit(m)
pub(crate) fn raise_spin_ground(
    u_value: u8,
    m_value: SpinMultiplicity,
    u_mode: &UnpairedElectronsMode,
    m_mode: &MultiplicityMode,
) -> (ValueAst, ValueAst) {
    let derived = m_value.multiplicity() == u_value + 1;

    let u_ast = match u_mode {
        UnpairedElectronsMode::Zero if u_value == 0 => ValueAst::Undetermined,
        UnpairedElectronsMode::Derived if derived => match m_mode {
            MultiplicityMode::Required => ValueAst::Undetermined,
            MultiplicityMode::Derived if u_value == 0 => ValueAst::Undetermined,
            MultiplicityMode::Derived => ValueAst::Lit(u_value as i64),
        },
        _ => ValueAst::Lit(u_value as i64),
    };

    let m_ast = match m_mode {
        MultiplicityMode::Derived if derived => ValueAst::Undetermined,
        _ => ValueAst::Lit(m_value.multiplicity() as i64),
    };

    (u_ast, m_ast)
}

/// Coupled spin-state raising for patterns.
///
/// Returns `(u_ast, m_ast)` — the minimal AST fields that roundtrip through `lower_spin`.
///
/// Suppression rules (Undetermined = field can be elided):
///   u pattern:
///     Is(0)  + Zero     → Undetermined
///     Any    + Required → Undetermined
///     Any    + Derived  + m=Any → Undetermined
///     Any    + Zero     → Undetermined
///     Any    + Derived  + m=Is(_) → Undetermined
///     Is(n)  + _        → Lit(n)
///   m pattern (given the u_pattern chosen above):
///     Any    + Required → Undetermined
///     Any    + Derived  + effective_u=Any → Undetermined
///     Any    + Derived  + effective_u=Is(_) → Undetermined
///     Is(m)  + Derived  + effective_u=Is(u) if m==u+1 → Undetermined
///     Is(m)  + _        → Lit(m)
pub(crate) fn raise_spin_pattern(
    u_pattern: Pattern<u8>,
    m_pattern: Pattern<SpinMultiplicity>,
    u_mode: &UnpairedElectronsMode,
    m_mode: &MultiplicityMode,
) -> (ValueAst, ValueAst) {
    let u_ast = match (u_pattern, u_mode) {
        (Pattern::Is(0), UnpairedElectronsMode::Zero) => ValueAst::Undetermined,
        (Pattern::Any, UnpairedElectronsMode::Required) => ValueAst::Undetermined,
        (Pattern::Any, UnpairedElectronsMode::Zero) => ValueAst::Undetermined,
        (Pattern::Any, UnpairedElectronsMode::Derived) => ValueAst::Undetermined,
        (Pattern::Is(n), _) => ValueAst::Lit(n as i64),
    };

    // Reconstruct what lower_spin will see after processing u_ast
    let effective_u = match &u_ast {
        ValueAst::Undetermined => match u_mode {
            UnpairedElectronsMode::Zero => Pattern::Is(0u8),
            UnpairedElectronsMode::Required | UnpairedElectronsMode::Derived => Pattern::Any,
        },
        ValueAst::Lit(n) => Pattern::Is(*n as u8),
        _ => unreachable!(),
    };

    let m_ast = match (m_pattern, m_mode) {
        (Pattern::Any, MultiplicityMode::Required) => ValueAst::Undetermined,
        (Pattern::Any, MultiplicityMode::Derived) => ValueAst::Undetermined,
        (Pattern::Is(m_val), MultiplicityMode::Derived) => match effective_u {
            Pattern::Is(u_val) if m_val.multiplicity() == u_val + 1 => ValueAst::Undetermined,
            _ => ValueAst::Lit(m_val.multiplicity() as i64),
        },
        (Pattern::Is(m_val), MultiplicityMode::Required) => {
            ValueAst::Lit(m_val.multiplicity() as i64)
        }
    };

    (u_ast, m_ast)
}
