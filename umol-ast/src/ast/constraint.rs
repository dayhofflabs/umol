//! AST constraints: per-scope predicates and their containers.
//!
//! Per-scope enums (`AtomConstraintAst`, `BondConstraintAst`, `DativeBondConstraintAst`,
//! `AromaticSystemConstraintAst`, `MulticenterBondConstraintAst`,
//! `NoncovalentBondConstraintAst`, `MoleculeConstraint`) each carry the predicates
//! admissible at that scope. `Constraint` is the tree node type admitting
//! per-entity leaves, a molecule-scope leaf, and `And`/`Or`/`Not` combinators.
//!
//! Per-entity constraints live inline on the entity AST via the typed
//! containers (`AtomConstraintsAst`, `BondConstraintsAst`, `DativeBondConstraintsAst`,
//! `AromaticSystemConstraintsAst`, `MulticenterBondConstraintsAst`,
//! `NoncovalentBondConstraintsAst`). Each exposes a uniform `new`/`len`/`iter`/
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
    AromaticSystemConstraintAst, AromaticSystemConstraintKey, AromaticSystemConstraintsAst,
};
pub use atom::{
    aromatic_covalence, AromaticValence, AromaticValenceAst, AtomConstraintAst, AtomConstraintKey,
    AtomConstraintsAst, MulticenterValence, MulticenterValenceAst,
};
pub use bond::{BondConstraintAst, BondConstraintKey, BondConstraintsAst};
pub use dative::{DativeBondConstraintAst, DativeBondConstraintKey, DativeBondConstraintsAst};
pub use molecule::{Constraint, Constraints, MoleculeConstraint, SubPatternAnchor};
pub use multicenter::{
    MulticenterBondConstraintAst, MulticenterBondConstraintKey, MulticenterBondConstraintsAst,
};
pub use noncovalent::{
    NoncovalentBondConstraintAst, NoncovalentBondConstraintKey, NoncovalentBondConstraintsAst,
};
pub use relational::RelationalConstraint;
pub use ring::{RingMembershipAst, RingScope};
pub use stereo::{
    FluxionalityAst, LigandPermutation, LigandSymmetryAst, OrientedLigandPermutation,
    StereoAtomConstraintAst, StereoAtomConstraintKey, StereoAtomConstraintsAst,
    StereoBondConstraintAst, StereoBondConstraintKey, StereoBondConstraintsAst, StereoLigandPair,
    StereogenicityAst, TopicityAst, TopicityRelationAst,
};
