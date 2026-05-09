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

pub use aromatic::{
    AromaticSystemConstraint, AromaticSystemConstraintKind, AromaticSystemConstraints,
};
pub use atom::{
    AromaticValenceAst, AtomConstraint, AtomConstraintKind, AtomConstraints, MulticenterValenceAst,
};
pub use bond::{BondConstraint, BondConstraintKind, BondConstraints};
pub use dative::{DativeBondConstraint, DativeBondConstraintKind, DativeBondConstraints};
pub use molecule::{Constraint, Constraints, MoleculeConstraint, SubPatternAnchor};
pub use multicenter::{MulticenterBondConstraint, MulticenterBondConstraints};
pub use noncovalent::{NoncovalentBondConstraint, NoncovalentBondConstraints};
pub use relational::RelationalConstraint;
