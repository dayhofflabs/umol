//! Semantic AST layer.
//!
//! `umol_ast::ast` is the canonical import location for every AST type.
//! Sub-modules are `pub(crate)` implementation detail; external users access
//! everything through this facade.

pub(crate) mod aromatic;
pub(crate) mod atom;
pub(crate) mod automorphism;
pub(crate) mod bond;
pub(crate) mod constraint;
pub(crate) mod dative;
pub(crate) mod error;
pub(crate) mod idx;
pub(crate) mod matching;
pub(crate) mod molecule;
pub(crate) mod multicenter;
pub(crate) mod noncovalent;
pub(crate) mod reaction;
pub(crate) mod remap;
pub(crate) mod rings;
pub(crate) mod spin;
pub(crate) mod subgraph;
pub(crate) mod traits;
pub(crate) mod value;
pub(crate) mod views;

pub use aromatic::AromaticSystemAst;
pub use atom::{AtomAst, ElementAst, ImplicitHydrogensAst, IsotopeAst};
pub use automorphism::AtomAutomorphism;
pub use bond::BondAst;
pub use constraint::{
    AromaticSystemConstraint, AromaticSystemConstraints, AromaticValenceAst, AtomConstraint,
    AtomConstraintKind, AtomConstraints, BondConstraint, BondConstraintKind, BondConstraints,
    Constraint, Constraints, DativeBondConstraint, DativeBondConstraintKind, DativeBondConstraints,
    MoleculeConstraint, MulticenterBondConstraint, MulticenterBondConstraints,
    MulticenterValenceAst, NoncovalentBondConstraint, NoncovalentBondConstraints, SubPatternAnchor,
};
pub use dative::{DativeBondAst, DativeDirection};
pub use error::{EvaluationError, RewriteError};
pub use idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
pub use matching::BondMatching;
pub use molecule::{MoleculeAst, MoleculeBuilder};
pub use multicenter::MulticenterBondAst;
pub use noncovalent::{NoncovalentBondAst, NoncovalentKind, NoncovalentKindAst};
pub use reaction::{Assignment, ReactionRuleAst};
pub use remap::IdxRemapping;
pub use rings::{RingFamily, RingGraph, RingGraphEdge, RingIdx, RingRelation, RingSet, RingView};
pub use spin::SpinStateAst;
pub use subgraph::MoleculeSubgraph;
pub use traits::{FromAst, IntoAst};
pub use value::{ArithOp, Bindings, Expr, RelOp, ValueAst};
pub use views::{
    AromaticSystemView, AromaticSystemViews, AtomView, AtomViewMut, AtomViews, BondView,
    BondViewMut, BondViews, DativeBondView, DativeBondViews, MulticenterBondView,
    MulticenterBondViews, NeighborView, NoncovalentBondView, NoncovalentBondViews,
};
