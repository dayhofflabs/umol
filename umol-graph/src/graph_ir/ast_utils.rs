//! Utility functions for converting between AST and IR types.

use umol_data::SpinMultiplicity;

use super::atom_pattern::Pattern;
use crate::dsl::config::{MultiplicityMode, NumericMode, UnpairedElectronsMode};
use crate::dsl::error::LoweringError;
use crate::dsl::value::ValueAst;

pub(crate) fn raise_u8_ground(value: u8, mode: &NumericMode) -> Option<ValueAst> {
    match (value, mode) {
        (0, NumericMode::Zero) => None,
        (n, _) => Some(ValueAst::Lit(n as i32)),
    }
}

pub(crate) fn raise_i8_ground(value: i8, mode: &NumericMode) -> Option<ValueAst> {
    match (value, mode) {
        (0, NumericMode::Zero) => None,
        (n, _) => Some(ValueAst::Lit(n as i32)),
    }
}

pub(crate) fn raise_u8_pattern(pattern: Pattern<u8>, mode: &NumericMode) -> Option<ValueAst> {
    match (pattern, mode) {
        (Pattern::Is(0), NumericMode::Zero) => None,
        (Pattern::Any, NumericMode::Required) => None,
        (Pattern::Any, NumericMode::Zero) => Some(ValueAst::Wildcard),
        (Pattern::Is(n), _) => Some(ValueAst::Lit(n as i32)),
    }
}

pub(crate) fn raise_i8_pattern(pattern: Pattern<i8>, mode: &NumericMode) -> Option<ValueAst> {
    match (pattern, mode) {
        (Pattern::Is(0), NumericMode::Zero) => None,
        (Pattern::Any, NumericMode::Required) => None,
        (Pattern::Any, NumericMode::Zero) => Some(ValueAst::Wildcard),
        (Pattern::Is(n), _) => Some(ValueAst::Lit(n as i32)),
    }
}

/// Coupled spin-state lowering for unpaired electrons (u) and multiplicity (m).
///
/// Resolution is sequential: u resolves first, m second. Each absent field
/// falls back to its mode:
///
///   u_ast present  → use it directly (Lit → Is, Wildcard → Any)
///   u_ast absent   → Zero: Is(0), Required: Any, Derived: derive from m_ast (Lit(m) → Is(m-1), else Any)
///
///   m_ast present  → use it directly
///   m_ast absent   → Required: Any, Derived: derive from resolved u (Is(u) → Is(u+1), Any → Any)
///
/// The asymmetry matters: u in Derived mode reads raw m_ast, but m in Derived
/// mode reads the already-resolved u pattern. This lets `C#s3` (m=3, u absent)
/// derive u=2, then confirm m=3 matches u+1.
pub(crate) fn lower_spin(
    u_ast: Option<ValueAst>,
    m_ast: Option<ValueAst>,
    u_mode: &UnpairedElectronsMode,
    m_mode: &MultiplicityMode,
) -> Result<(Pattern<u8>, Pattern<SpinMultiplicity>), LoweringError> {
    // Step 1: resolve unpaired electrons
    let unpaired_electrons = match &u_ast {
        Some(ValueAst::Wildcard) => Pattern::Any,
        Some(ValueAst::Lit(n)) => {
            Pattern::Is(u8::try_from(*n).map_err(|_| LoweringError::NonGround {
                field: "unpaired_electrons",
            })?)
        }
        Some(_) => {
            return Err(LoweringError::NonGround {
                field: "unpaired_electrons",
            })
        }
        None => match u_mode {
            UnpairedElectronsMode::Zero => Pattern::Is(0),
            UnpairedElectronsMode::Required => Pattern::Any,
            UnpairedElectronsMode::Derived => match &m_ast {
                Some(ValueAst::Lit(m)) => {
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
    };

    // Step 2: resolve multiplicity (derives from resolved u when mode is Derived)
    let multiplicity = match &m_ast {
        Some(ValueAst::Wildcard) => Pattern::Any,
        Some(ValueAst::Lit(n)) => {
            let m = u8::try_from(*n).map_err(|_| LoweringError::NonGround {
                field: "multiplicity",
            })?;
            Pattern::Is(
                SpinMultiplicity::from_multiplicity(m)
                    .ok_or(LoweringError::InvalidMultiplicity(m))?,
            )
        }
        Some(_) => {
            return Err(LoweringError::NonGround {
                field: "multiplicity",
            })
        }
        None => match m_mode {
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
    };

    Ok((unpaired_electrons, multiplicity))
}

/// Coupled spin-state raising for ground types.
///
/// Returns `(u_ast, m_ast)` — the minimal AST fields that roundtrip through `lower_spin`.
///
/// `derived` = m == u + 1 (high-spin / Hund's rule relationship).
///
/// Suppression rules:
///   u side:
///     Zero    + u==0                        → None  (default produces 0)
///     Derived + derived + Required           → None  (m always emitted; u derivable from m)
///     Derived + derived + Derived + u==0     → None  (both absent → defaults match)
///     Derived + derived + Derived + u!=0     → Lit(u) (must fix u; m derived from u)
///     otherwise                              → Lit(u)
///   m side:
///     Derived + derived → None  (derivable from u or default)
///     otherwise         → Lit(m)
pub(crate) fn raise_spin_ground(
    u_value: u8,
    m_value: SpinMultiplicity,
    u_mode: &UnpairedElectronsMode,
    m_mode: &MultiplicityMode,
) -> (Option<ValueAst>, Option<ValueAst>) {
    let derived = m_value.multiplicity() == u_value + 1;

    let u_ast = match u_mode {
        UnpairedElectronsMode::Zero if u_value == 0 => None,
        UnpairedElectronsMode::Derived if derived => match m_mode {
            MultiplicityMode::Required => None,
            MultiplicityMode::Derived if u_value == 0 => None,
            MultiplicityMode::Derived => Some(ValueAst::Lit(u_value as i32)),
        },
        _ => Some(ValueAst::Lit(u_value as i32)),
    };

    let m_ast = match m_mode {
        MultiplicityMode::Derived if derived => None,
        _ => Some(ValueAst::Lit(m_value.multiplicity() as i32)),
    };

    (u_ast, m_ast)
}

/// Coupled spin-state raising for patterns.
///
/// Returns `(u_ast, m_ast)` — the minimal AST fields that roundtrip through `lower_spin`.
///
/// Suppression rules:
///   u pattern:
///     Is(0)  + Zero     → None          (default produces Is(0))
///     Any    + Required → None          (default produces Any)
///     Any    + Derived  + m=Any → None  (both absent → both Any)
///     Any    + Zero     → Wildcard      (need explicit marker)
///     Any    + Derived  + m=Is(_) → Wildcard  (m will be emitted; must mark u as non-ground)
///     Is(n)  + _        → Lit(n)
///   m pattern (given the u_pattern chosen above):
///     Any    + Required → None          (default produces Any)
///     Any    + Derived  + effective_u=Any → None   (derived from Any → Any)
///     Any    + Derived  + effective_u=Is(_) → Wildcard  (derived would give Is, not Any)
///     Is(m)  + Derived  + effective_u=Is(u) if m==u+1 → None  (derivable)
///     Is(m)  + _        → Lit(m)
pub(crate) fn raise_spin_pattern(
    u_pattern: Pattern<u8>,
    m_pattern: Pattern<SpinMultiplicity>,
    u_mode: &UnpairedElectronsMode,
    m_mode: &MultiplicityMode,
) -> (Option<ValueAst>, Option<ValueAst>) {
    let u_ast = match (u_pattern, u_mode) {
        (Pattern::Is(0), UnpairedElectronsMode::Zero) => None,
        (Pattern::Any, UnpairedElectronsMode::Required) => None,
        (Pattern::Any, UnpairedElectronsMode::Zero) => Some(ValueAst::Wildcard),
        (Pattern::Any, UnpairedElectronsMode::Derived) => match m_pattern {
            Pattern::Any => None,
            Pattern::Is(_) => Some(ValueAst::Wildcard),
        },
        (Pattern::Is(n), _) => Some(ValueAst::Lit(n as i32)),
    };

    // Reconstruct what lower_spin will see after processing u_ast
    let effective_u = match &u_ast {
        None => match u_mode {
            UnpairedElectronsMode::Zero => Pattern::Is(0u8),
            UnpairedElectronsMode::Required | UnpairedElectronsMode::Derived => Pattern::Any,
        },
        Some(ValueAst::Wildcard) => Pattern::Any,
        Some(ValueAst::Lit(n)) => Pattern::Is(*n as u8),
        _ => unreachable!(),
    };

    let m_ast = match (m_pattern, m_mode) {
        (Pattern::Any, MultiplicityMode::Required) => None,
        (Pattern::Any, MultiplicityMode::Derived) => match effective_u {
            Pattern::Any => None,
            Pattern::Is(_) => Some(ValueAst::Wildcard),
        },
        (Pattern::Is(m_val), MultiplicityMode::Derived) => match effective_u {
            Pattern::Is(u_val) if m_val.multiplicity() == u_val + 1 => None,
            _ => Some(ValueAst::Lit(m_val.multiplicity() as i32)),
        },
        (Pattern::Is(m_val), MultiplicityMode::Required) => {
            Some(ValueAst::Lit(m_val.multiplicity() as i32))
        }
    };

    (u_ast, m_ast)
}
