//! Model-independent constraints derived from entity fields and directly incident entities.

use thiserror::Error;

use super::super::super::constraint::{
    AromaticSystemConstraintAst, AtomConstraintAst, BondConstraintAst, DativeBondConstraintAst,
    MulticenterBondConstraintAst, NoncovalentBondConstraintAst,
};
use super::super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};

/// Evaluates model-independent incidence constraints without running a graph algorithm.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IncidenceConstraintValidator;

#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum IncidenceConstraintContradiction {
    #[error("atom {atom:?} does not satisfy incidence constraint {constraint:?}")]
    Atom {
        atom: AtomId,
        constraint: AtomConstraintAst,
    },
    #[error("bond {bond:?} does not satisfy incidence constraint {constraint:?}")]
    Bond {
        bond: BondId,
        constraint: BondConstraintAst,
    },
    #[error("dative bond {bond:?} does not satisfy incidence constraint {constraint:?}")]
    DativeBond {
        bond: DativeBondId,
        constraint: DativeBondConstraintAst,
    },
    #[error("aromatic system {system:?} does not satisfy incidence constraint {constraint:?}")]
    AromaticSystem {
        system: AromaticSystemId,
        constraint: AromaticSystemConstraintAst,
    },
    #[error("multicenter bond {bond:?} does not satisfy incidence constraint {constraint:?}")]
    MulticenterBond {
        bond: MulticenterBondId,
        constraint: MulticenterBondConstraintAst,
    },
    #[error("noncovalent bond {bond:?} does not satisfy incidence constraint {constraint:?}")]
    NoncovalentBond {
        bond: NoncovalentBondId,
        constraint: NoncovalentBondConstraintAst,
    },
}
