//! Structural AST: a generalization that can express both ground molecules and
//! pattern queries. AST types are independent of any particular text format;
//! the [`crate::dsl`] module provides one EDN-based serialization.

pub mod atom;
pub mod bond;
pub mod config;
pub mod constraint;
pub mod error;
pub mod matcher;
pub mod molecule;
pub mod morgan;
pub mod rings;
pub mod symmetry;

use umol_edn::{FromEdn, ToEdn};
use umol_graph_core::{EdgeId, NodeId};
use umol_graph_core::relation::RelationId;

macro_rules! define_idx {
    ($($(#[doc = $doc:literal])* $name:ident),+ $(,)?) => {
        $(
            $(#[doc = $doc])*
            #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, FromEdn, ToEdn)]
            #[edn(transparent)]
            pub struct $name(pub u32);

            impl From<usize> for $name {
                fn from(v: usize) -> Self { Self(v as u32) }
            }

            impl $name {
                pub fn from_usize(v: usize) -> Self { Self(v as u32) }
                pub fn index(self) -> usize { self.0 as usize }
            }
        )+
    };
}

define_idx!(
    /// Atom index — maps directly to `NodeId` in the underlying graph.
    AtomIdx,
    /// Bond index — maps directly to `EdgeId` in the underlying graph.
    BondIdx,
    /// Dative bond index — maps to `RelationId` in the dative bonds relation set.
    DativeBondIdx,
    /// Aromatic system index — maps to `RelationId` in the aromatic systems relation set.
    AromaticSystemIdx,
    /// Multicenter bond index — maps to `RelationId` in the multicenter bonds relation set.
    MulticenterBondIdx,
    /// Noncovalent bond index — maps to `RelationId` in the noncovalent bonds relation set.
    NoncovalentBondIdx,
);

impl From<NodeId> for AtomIdx {
    fn from(id: NodeId) -> Self { Self(id.0) }
}
impl From<AtomIdx> for NodeId {
    fn from(idx: AtomIdx) -> Self { Self(idx.0) }
}

impl From<EdgeId> for BondIdx {
    fn from(id: EdgeId) -> Self { Self(id.0) }
}
impl From<BondIdx> for EdgeId {
    fn from(idx: BondIdx) -> Self { Self(idx.0) }
}

macro_rules! relation_idx_from {
    ($name:ident) => {
        impl From<RelationId> for $name {
            fn from(id: RelationId) -> Self { Self(id.0) }
        }
        impl From<$name> for RelationId {
            fn from(idx: $name) -> Self { Self(idx.0) }
        }
    };
}

relation_idx_from!(DativeBondIdx);
relation_idx_from!(AromaticSystemIdx);
relation_idx_from!(MulticenterBondIdx);
relation_idx_from!(NoncovalentBondIdx);
