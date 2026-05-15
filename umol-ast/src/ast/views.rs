//! Read-only views over `MoleculeAst` topology and relations.
//!
//! View records bundle an index with the underlying data so consumers
//! never assemble (id, data, participants) tuples by hand. Namespace
//! types group per-relation accessors (`count`, `ids`, `iter`, `get`,
//! and `Index`) without burying them on `MoleculeAst` itself.

mod aromatic_system;
mod atom;
mod bond;
mod dative_bond;
mod graph;
mod multicenter_bond;
mod neighbor;
mod noncovalent_bond;

pub use aromatic_system::{
    AromaticSystemBuilderView, AromaticSystemBuilderViewMut, AromaticSystemView,
    AromaticSystemViews,
};
pub use atom::{AtomBuilderView, AtomBuilderViewMut, AtomView, AtomViewMut, AtomViews};
pub use bond::{BondBuilderView, BondBuilderViewMut, BondView, BondViewMut, BondViews};
pub use dative_bond::{
    DativeBondBuilderView, DativeBondBuilderViewMut, DativeBondView, DativeBondViews,
};
pub use graph::GraphView;
pub use multicenter_bond::{
    MulticenterBondBuilderView, MulticenterBondBuilderViewMut, MulticenterBondView,
    MulticenterBondViews,
};
pub use neighbor::NeighborView;
pub use noncovalent_bond::{
    NoncovalentBondBuilderView, NoncovalentBondBuilderViewMut, NoncovalentBondView,
    NoncovalentBondViews,
};
