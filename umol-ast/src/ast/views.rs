//! Read-only views over `MoleculeAst` topology and relations.
//!
//! View records bundle an index with the underlying data so consumers
//! never assemble (id, data, participants) tuples by hand. Namespace
//! types group per-relation accessors (`count`, `ids`, `iter`, `get`,
//! and `Index`) without burying them on `MoleculeAst` itself.

mod aromatic;
mod atom;
mod bond;
mod dative;
mod graph;
mod ligand;
mod multicenter;
mod neighbor;
mod noncovalent;
mod stereo;

pub use aromatic::{
    AromaticSystemBuilderView, AromaticSystemBuilderViewMut, AromaticSystemView,
    AromaticSystemViews,
};
pub use atom::{AtomBuilderView, AtomBuilderViewMut, AtomView, AtomViewMut, AtomViews};
pub use bond::{BondBuilderView, BondBuilderViewMut, BondView, BondViewMut, BondViews};
pub use dative::{
    DativeBondBuilderView, DativeBondBuilderViewMut, DativeBondView, DativeBondViews,
};
pub use graph::GraphView;
pub use ligand::StereoLigandView;
pub use multicenter::{
    MulticenterBondBuilderView, MulticenterBondBuilderViewMut, MulticenterBondView,
    MulticenterBondViews,
};
pub use neighbor::NeighborView;
pub use noncovalent::{
    NoncovalentBondBuilderView, NoncovalentBondBuilderViewMut, NoncovalentBondView,
    NoncovalentBondViews,
};
pub use stereo::{
    StereoAtomBuilderView, StereoAtomBuilderViewMut, StereoAtomView, StereoAtomViews,
    StereoBondBuilderView, StereoBondBuilderViewMut, StereoBondView, StereoBondViews,
};
