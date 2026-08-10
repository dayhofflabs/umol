//! AST constraints: per-scope predicates and their containers.
//!
//! Per-scope enums (`AtomConstraintForm`, `BondConstraintForm`, `DativeBondConstraintForm`,
//! `AromaticSystemConstraintForm`, `MulticenterBondConstraintForm`,
//! `NoncovalentBondConstraintForm`, `MoleculeConstraint`) each carry the predicates
//! admissible at that scope. `Constraint` is the tree node type admitting
//! per-entity leaves, a molecule-scope leaf, and `And`/`Or`/`Not` combinators.
//!
//! Per-entity constraints live inline on the entity AST via the typed
//! containers (`AtomConstraintsForm`, `BondConstraintsForm`, `DativeBondConstraintsForm`,
//! `AromaticSystemConstraintsForm`, `MulticenterBondConstraintsForm`,
//! `NoncovalentBondConstraintsForm`). Each exposes a uniform `new`/`len`/`iter`/
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
    AromaticSystemConstraintForm, AromaticSystemConstraintKey, AromaticSystemConstraintsForm,
};
pub use atom::{
    aromatic_covalence, AromaticValence, AromaticValenceForm, AtomConstraintForm,
    AtomConstraintKey, AtomConstraintsForm, MulticenterValence, MulticenterValenceForm,
};
pub use bond::{BondConstraintForm, BondConstraintKey, BondConstraintsForm};
pub use dative::{DativeBondConstraintForm, DativeBondConstraintKey, DativeBondConstraintsForm};
pub use molecule::{Constraint, Constraints, MoleculeConstraint};
pub use multicenter::{
    MulticenterBondConstraintForm, MulticenterBondConstraintKey, MulticenterBondConstraintsForm,
};
pub use noncovalent::{
    NoncovalentBondConstraintForm, NoncovalentBondConstraintKey, NoncovalentBondConstraintsForm,
};
pub use relational::RelationalConstraint;
pub use ring::{RingMembershipForm, RingScope};
pub use stereo::{
    FluxionalityForm, LigandPermutation, LigandSymmetryForm, OrientedLigandPermutation,
    StereoAtomConstraintForm, StereoAtomConstraintKey, StereoAtomConstraintsForm,
    StereoBondConstraintForm, StereoBondConstraintKey, StereoBondConstraintsForm, StereoLigandPair,
    StereogenicityForm, TopicityForm, TopicityRelationForm,
};
