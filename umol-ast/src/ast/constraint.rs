//! AST constraints: per-scope predicates and their routing.
//!
//! Per-scope enums (`AtomConstraint`, `BondConstraint`, `DativeBondConstraint`,
//! `AromaticSystemConstraint`, `MulticenterBondConstraint`,
//! `NoncovalentBondConstraint`, `MoleculeConstraint`) each carry the predicates
//! admissible at that scope. `Constraint` is the tree node type admitting
//! per-entity leaves, a molecule-scope leaf, and `And`/`Or`/`Not` combinators.
//!
//! Storage is dual. Each atom carries an `AtomConstraints` slotmap keyed by
//! `AtomConstraintKind`; other entity ASTs carry per-kind inline vecs. The
//! molecule-level `Constraints` (under `molecule`) carries per-scope
//! `IndexMap` buckets plus a flat `molecule` vec for molecule-scope and
//! combinator forms. Consumers read the union of inline and molecule-level
//! entries for any given idx; there is no invariant between the two stores.

pub mod aromatic;
pub mod atom;
pub mod bond;
pub mod dative;
pub mod molecule;
pub mod multicenter;
pub mod noncovalent;

pub use aromatic::AromaticSystemConstraint;
pub use atom::{
    AromaticValenceConstraint, AtomConstraint, AtomConstraintKind, AtomConstraints,
    MulticenterValenceConstraint,
};
pub use bond::{BondConstraint, BondConstraintKind};
pub use dative::{DativeBondConstraint, DativeBondConstraintKind};
pub use molecule::{Constraint, Constraints, MoleculeConstraint, SubPatternAnchor};
pub use multicenter::MulticenterBondConstraint;
pub use noncovalent::NoncovalentBondConstraint;
