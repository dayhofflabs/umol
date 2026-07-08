//! AST constraints: per-scope predicates and their containers.
//!
//! Per-scope enums (`AtomConstraint`, `BondConstraint`, `DativeBondConstraint`,
//! `AromaticSystemConstraint`, `MulticenterBondConstraint`,
//! `NoncovalentBondConstraint`, `MoleculeConstraint`) each carry the predicates
//! admissible at that scope. `Constraint` is the tree node type admitting
//! per-entity leaves, a molecule-scope leaf, and `And`/`Or`/`Not` combinators.
//!
//! Per-entity constraints live inline on the entity AST via the typed
//! containers (`AtomConstraints`, `BondConstraints`, `DativeBondConstraints`,
//! `AromaticSystemConstraints`, `MulticenterBondConstraints`,
//! `NoncovalentBondConstraints`). Each exposes a uniform `new`/`len`/`iter`/
//! `add`/`retain`/`clear` surface; `add` enforces per-variant cardinality
//! (last-wins for unique-kind variants, append for multi-kind variants).
//! The molecule-level `Constraints` (under `molecule`) is a flat
//! `Vec<Constraint>` for molecule-scope predicates and combinator forms.

mod aromatic;
mod atom;
mod bond;
mod dative;
mod molecule;
mod multicenter;
mod noncovalent;
mod relational;
mod ring;
mod stereo;

pub use aromatic::{
    AromaticSystemConstraint, AromaticSystemConstraintKey, AromaticSystemConstraints,
};
pub use atom::{
    aromatic_increment, AromaticValenceAst, AtomConstraint, AtomConstraintKey, AtomConstraints,
    MulticenterValenceAst,
};
pub use bond::{BondConstraint, BondConstraintKey, BondConstraints};
pub use dative::{DativeBondConstraint, DativeBondConstraintKey, DativeBondConstraints};
pub use molecule::{Constraint, Constraints, MoleculeConstraint, SubPatternAnchor};
pub use multicenter::{
    MulticenterBondConstraint, MulticenterBondConstraintKey, MulticenterBondConstraints,
};
pub use noncovalent::{
    NoncovalentBondConstraint, NoncovalentBondConstraintKey, NoncovalentBondConstraints,
};
pub use relational::RelationalConstraint;
pub use ring::{RingMembershipAst, RingScope};
pub use stereo::{
    FluxionalityAst, LigandPermutation, LigandSymmetryAst, OrientedLigandPermutation,
    StereoAtomConstraint, StereoAtomConstraintKey, StereoAtomConstraints, StereoBondConstraint,
    StereoBondConstraintKey, StereoBondConstraints, StereoLigandPair, StereogenicityAst,
    TopicityAst, TopicityRelationAst,
};
