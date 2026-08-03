//! Python edit values and symbolic creation handles.

use pyo3::prelude::*;
use umol_ast::ast::{
    AromaticSystemHandle as AstAromaticSystemHandle, AromaticSystemId as AstAromaticSystemId,
    AtomHandle as AstAtomHandle, AtomId as AstAtomId, BondHandle as AstBondHandle,
    BondId as AstBondId, DativeBondHandle as AstDativeBondHandle, DativeBondId as AstDativeBondId,
    MulticenterBondHandle as AstMulticenterBondHandle, MulticenterBondId as AstMulticenterBondId,
    NoncovalentBondHandle as AstNoncovalentBondHandle, NoncovalentBondId as AstNoncovalentBondId,
    StereoAtomHandle as AstStereoAtomHandle, StereoAtomId as AstStereoAtomId,
    StereoBondHandle as AstStereoBondHandle, StereoBondId as AstStereoBondId,
};

/// A same-kind creation ordinal in an edit sequence.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct New {
    index: usize,
}

#[pymethods]
impl New {
    #[new]
    fn from_index(index: usize) -> Self {
        Self { index }
    }

    /// Zero-based creation ordinal within the surrounding entity kind.
    #[getter]
    fn index(&self) -> usize {
        self.index
    }

    fn __repr__(&self) -> String {
        format!("New({})", self.index)
    }
}

/// A Python host id or same-kind creation ordinal. The argument position supplies the entity kind.
#[derive(FromPyObject)]
#[allow(
    dead_code,
    reason = "Python-to-Rust handle input for the Edit and Edits bindings"
)]
pub(crate) enum HandleLike {
    New(New),
    Id(u32),
}

#[allow(
    dead_code,
    reason = "Python-to-Rust handle conversions for the Edit and Edits bindings"
)]
impl HandleLike {
    pub(crate) fn to_atom_handle(&self) -> AstAtomHandle {
        match self {
            Self::Id(index) => AstAtomHandle::Id(AstAtomId(*index)),
            Self::New(new) => AstAtomHandle::New(new.index),
        }
    }

    pub(crate) fn to_bond_handle(&self) -> AstBondHandle {
        match self {
            Self::Id(index) => AstBondHandle::Id(AstBondId(*index)),
            Self::New(new) => AstBondHandle::New(new.index),
        }
    }

    pub(crate) fn to_dative_bond_handle(&self) -> AstDativeBondHandle {
        match self {
            Self::Id(index) => AstDativeBondHandle::Id(AstDativeBondId(*index)),
            Self::New(new) => AstDativeBondHandle::New(new.index),
        }
    }

    pub(crate) fn to_aromatic_system_handle(&self) -> AstAromaticSystemHandle {
        match self {
            Self::Id(index) => AstAromaticSystemHandle::Id(AstAromaticSystemId(*index)),
            Self::New(new) => AstAromaticSystemHandle::New(new.index),
        }
    }

    pub(crate) fn to_multicenter_bond_handle(&self) -> AstMulticenterBondHandle {
        match self {
            Self::Id(index) => AstMulticenterBondHandle::Id(AstMulticenterBondId(*index)),
            Self::New(new) => AstMulticenterBondHandle::New(new.index),
        }
    }

    pub(crate) fn to_noncovalent_bond_handle(&self) -> AstNoncovalentBondHandle {
        match self {
            Self::Id(index) => AstNoncovalentBondHandle::Id(AstNoncovalentBondId(*index)),
            Self::New(new) => AstNoncovalentBondHandle::New(new.index),
        }
    }

    pub(crate) fn to_stereo_atom_handle(&self) -> AstStereoAtomHandle {
        match self {
            Self::Id(index) => AstStereoAtomHandle::Id(AstStereoAtomId(*index)),
            Self::New(new) => AstStereoAtomHandle::New(new.index),
        }
    }

    pub(crate) fn to_stereo_bond_handle(&self) -> AstStereoBondHandle {
        match self {
            Self::Id(index) => AstStereoBondHandle::Id(AstStereoBondId(*index)),
            Self::New(new) => AstStereoBondHandle::New(new.index),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::id(HandleLike::Id(7), AstAtomHandle::Id(AstAtomId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), AstAtomHandle::New(7))]
    fn test_handle_like_to_atom_handle(#[case] input: HandleLike, #[case] expected: AstAtomHandle) {
        assert_eq!(input.to_atom_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), AstBondHandle::Id(AstBondId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), AstBondHandle::New(7))]
    fn test_handle_like_to_bond_handle(#[case] input: HandleLike, #[case] expected: AstBondHandle) {
        assert_eq!(input.to_bond_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), AstDativeBondHandle::Id(AstDativeBondId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), AstDativeBondHandle::New(7))]
    fn test_handle_like_to_dative_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: AstDativeBondHandle,
    ) {
        assert_eq!(input.to_dative_bond_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), AstAromaticSystemHandle::Id(AstAromaticSystemId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), AstAromaticSystemHandle::New(7))]
    fn test_handle_like_to_aromatic_system_handle(
        #[case] input: HandleLike,
        #[case] expected: AstAromaticSystemHandle,
    ) {
        assert_eq!(input.to_aromatic_system_handle(), expected);
    }

    #[rstest]
    #[case::id(
        HandleLike::Id(7),
        AstMulticenterBondHandle::Id(AstMulticenterBondId(7))
    )]
    #[case::new(HandleLike::New(New { index: 7 }), AstMulticenterBondHandle::New(7))]
    fn test_handle_like_to_multicenter_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: AstMulticenterBondHandle,
    ) {
        assert_eq!(input.to_multicenter_bond_handle(), expected);
    }

    #[rstest]
    #[case::id(
        HandleLike::Id(7),
        AstNoncovalentBondHandle::Id(AstNoncovalentBondId(7))
    )]
    #[case::new(HandleLike::New(New { index: 7 }), AstNoncovalentBondHandle::New(7))]
    fn test_handle_like_to_noncovalent_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: AstNoncovalentBondHandle,
    ) {
        assert_eq!(input.to_noncovalent_bond_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), AstStereoAtomHandle::Id(AstStereoAtomId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), AstStereoAtomHandle::New(7))]
    fn test_handle_like_to_stereo_atom_handle(
        #[case] input: HandleLike,
        #[case] expected: AstStereoAtomHandle,
    ) {
        assert_eq!(input.to_stereo_atom_handle(), expected);
    }

    #[rstest]
    #[case::id(HandleLike::Id(7), AstStereoBondHandle::Id(AstStereoBondId(7)))]
    #[case::new(HandleLike::New(New { index: 7 }), AstStereoBondHandle::New(7))]
    fn test_handle_like_to_stereo_bond_handle(
        #[case] input: HandleLike,
        #[case] expected: AstStereoBondHandle,
    ) {
        assert_eq!(input.to_stereo_bond_handle(), expected);
    }
}
