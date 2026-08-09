//! Read-only views over `Molecule` topology and relations.
//!
//! View records bundle an index with the underlying data so consumers
//! never assemble (id, data, participants) tuples by hand. Namespace
//! types group per-relation accessors (`count`, `ids`, `iter`, `get`,
//! and `Index`) without burying them on `Molecule` itself.

#[cfg(test)]
use std::fmt::Debug;

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

#[cfg(test)]
fn assert_exact_size_by<I, T, U>(mut iterator: I, expected: Vec<U>, mut project: impl FnMut(T) -> U)
where
    I: ExactSizeIterator<Item = T>,
    U: Debug + PartialEq,
{
    assert_eq!(iterator.len(), expected.len());
    assert_eq!(iterator.size_hint(), (expected.len(), Some(expected.len())));
    while let Some(expected_item) = expected.get(expected.len() - iterator.len()) {
        let previous = iterator.len();
        assert_eq!(
            iterator.next().map(&mut project).as_ref(),
            Some(expected_item),
        );
        let remaining = iterator.len();
        assert_eq!(remaining, previous - 1);
        assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
    }
    assert_eq!(iterator.next().map(project), None);
    assert_eq!(iterator.len(), 0);
    assert_eq!(iterator.size_hint(), (0, Some(0)));
}
