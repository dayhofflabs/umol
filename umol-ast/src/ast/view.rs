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
mod ring;
mod stereo;

pub use aromatic::{
    AromaticSystemEditorView, AromaticSystemEditorViewMut, AromaticSystemView,
    AromaticSystemViewMut, AromaticSystemViews,
};
pub use atom::{AtomEditorView, AtomEditorViewMut, AtomView, AtomViewMut, AtomViews};
pub use bond::{BondEditorView, BondEditorViewMut, BondView, BondViewMut, BondViews};
pub use dative::{
    DativeBondEditorView, DativeBondEditorViewMut, DativeBondView, DativeBondViewMut,
    DativeBondViews,
};
pub use graph::{AtomAutomorphism, GraphView};
pub use ligand::StereoLigandView;
pub use multicenter::{
    MulticenterBondEditorView, MulticenterBondEditorViewMut, MulticenterBondView,
    MulticenterBondViewMut, MulticenterBondViews,
};
pub use neighbor::NeighborView;
pub use noncovalent::{
    NoncovalentBondEditorView, NoncovalentBondEditorViewMut, NoncovalentBondView,
    NoncovalentBondViewMut, NoncovalentBondViews,
};
pub use ring::{RingAtomView, RingBondView, RingView, RingViews};
pub use stereo::{
    StereoAtomEditorView, StereoAtomEditorViewMut, StereoAtomView, StereoAtomViewMut,
    StereoAtomViews, StereoBondEditorView, StereoBondEditorViewMut, StereoBondView,
    StereoBondViewMut, StereoBondViews,
};
