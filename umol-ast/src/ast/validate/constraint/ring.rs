//! Constraints evaluated against the fixed Relevant ring projection through size 22.

use thiserror::Error;

use super::super::super::constraint::{AtomConstraintAst, BondConstraintAst};
use super::super::super::id::{AtomId, BondId};

/// Evaluates ring constraints with an explicit relevant-cycle algorithm selector.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RingConstraintValidator;

#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum RingConstraintContradiction {
    #[error("atom {atom:?} does not satisfy ring constraint {constraint:?}")]
    Atom {
        atom: AtomId,
        constraint: AtomConstraintAst,
    },
    #[error("bond {bond:?} does not satisfy ring constraint {constraint:?}")]
    Bond {
        bond: BondId,
        constraint: BondConstraintAst,
    },
}
