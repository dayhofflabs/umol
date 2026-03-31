//! Utility functions for converting between AST and IR types.

use umol_data::SpinMultiplicity;

use super::atom_pattern::Pattern;
use crate::dsl::config::{MultiplicityMode, NumericMode, UnpairedElectronsMode};
use crate::dsl::error::LoweringError;
use crate::dsl::value::ValueAst;

pub(crate) fn raise_pattern_u8(pat: Pattern<u8>, mode: &NumericMode) -> Option<ValueAst> {
    match (pat, mode) {
        (Pattern::Is(0), NumericMode::Zero) => None,
        (Pattern::Any, NumericMode::Required) => None,
        (Pattern::Any, NumericMode::Zero) => Some(ValueAst::Wildcard),
        (Pattern::Is(n), _) => Some(ValueAst::Lit(n as i32)),
    }
}

pub(crate) fn raise_pattern_i8(pat: Pattern<i8>, mode: &NumericMode) -> Option<ValueAst> {
    match (pat, mode) {
        (Pattern::Is(0), NumericMode::Zero) => None,
        (Pattern::Any, NumericMode::Required) => None,
        (Pattern::Any, NumericMode::Zero) => Some(ValueAst::Wildcard),
        (Pattern::Is(n), _) => Some(ValueAst::Lit(n as i32)),
    }
}

/// Coupled spin-state lowering: resolve u first (may derive from raw m), then m (may derive from resolved u).
pub(crate) fn lower_spin(
    raw_u: Option<ValueAst>,
    raw_m: Option<ValueAst>,
    u_mode: &UnpairedElectronsMode,
    m_mode: &MultiplicityMode,
) -> Result<(Pattern<u8>, Pattern<SpinMultiplicity>), LoweringError> {
    // Step 1: resolve unpaired electrons
    let unpaired_electrons = match &raw_u {
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
            UnpairedElectronsMode::Derived => match &raw_m {
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
    let multiplicity = match &raw_m {
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

pub(crate) fn raise_spin_u_pattern(
    u: Pattern<u8>,
    m: Pattern<SpinMultiplicity>,
    u_mode: &UnpairedElectronsMode,
) -> Option<ValueAst> {
    match (u, u_mode) {
        (Pattern::Is(0), UnpairedElectronsMode::Zero) => None,
        (Pattern::Any, UnpairedElectronsMode::Required) => None,
        (Pattern::Any, UnpairedElectronsMode::Zero) => Some(ValueAst::Wildcard),
        (Pattern::Any, UnpairedElectronsMode::Derived) => match m {
            Pattern::Any => None,                       // both absent → both Any
            Pattern::Is(_) => Some(ValueAst::Wildcard), // m emitted, must mark u as wildcard
        },
        (Pattern::Is(n), _) => Some(ValueAst::Lit(n as i32)),
    }
}

pub(crate) fn raise_spin_m_pattern(
    u: Pattern<u8>,
    m: Pattern<SpinMultiplicity>,
    u_mode: &UnpairedElectronsMode,
    m_mode: &MultiplicityMode,
) -> Option<ValueAst> {
    // Compute effective_u as seen by from_ast after processing u_ast
    let effective_u = match raise_spin_u_pattern(u, m, u_mode) {
        None => match u_mode {
            UnpairedElectronsMode::Zero => Pattern::Is(0u8),
            UnpairedElectronsMode::Required | UnpairedElectronsMode::Derived => Pattern::Any,
        },
        Some(ValueAst::Wildcard) => Pattern::Any,
        Some(ValueAst::Lit(n)) => Pattern::Is(n as u8),
        _ => unreachable!(),
    };

    match (m, m_mode) {
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
    }
}

pub(crate) fn raise_u8(value: u8, mode: &NumericMode) -> Option<ValueAst> {
    match (value, mode) {
        (0, NumericMode::Zero) => None,
        (n, _) => Some(ValueAst::Lit(n as i32)),
    }
}

pub(crate) fn raise_i8(value: i8, mode: &NumericMode) -> Option<ValueAst> {
    match (value, mode) {
        (0, NumericMode::Zero) => None,
        (n, _) => Some(ValueAst::Lit(n as i32)),
    }
}

pub(crate) fn raise_spin_u_ground(
    u: u8,
    m: SpinMultiplicity,
    u_mode: &UnpairedElectronsMode,
    m_mode: &MultiplicityMode,
) -> Option<ValueAst> {
    let derived = m.multiplicity() == u + 1;
    match u_mode {
        UnpairedElectronsMode::Zero if u == 0 => None,
        UnpairedElectronsMode::Derived if derived => match m_mode {
            MultiplicityMode::Required => None,
            MultiplicityMode::Derived if u == 0 => None,
            MultiplicityMode::Derived => Some(ValueAst::Lit(u as i32)),
        },
        _ => Some(ValueAst::Lit(u as i32)),
    }
}

pub(crate) fn raise_spin_m_ground(
    u: u8,
    m: SpinMultiplicity,
    _u_mode: &UnpairedElectronsMode,
    m_mode: &MultiplicityMode,
) -> Option<ValueAst> {
    let derived = m.multiplicity() == u + 1;
    match m_mode {
        MultiplicityMode::Derived if derived => None,
        _ => Some(ValueAst::Lit(m.multiplicity() as i32)),
    }
}
