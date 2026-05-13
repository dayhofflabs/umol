//! AST indices into atom, bond, and relation tables.

use std::fmt;

use umol_edn::{FromEdn, ToEdn};
use umol_graph_core::{EdgeId, NodeId, RelationId};

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
                pub fn index(self) -> usize { self.0 as usize }
            }

            impl fmt::Display for $name {
                fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    write!(f, "{}", self.0)
                }
            }
        )+
    };
}

define_idx!(
    /// Atom index — maps directly to `NodeId` in the underlying graph.
    AtomId,
    /// Bond index — maps directly to `EdgeId` in the underlying graph.
    BondId,
    /// Dative bond index — maps to `RelationId` in the dative bonds relation set.
    DativeBondId,
    /// Aromatic system index — maps to `RelationId` in the aromatic systems relation set.
    AromaticSystemId,
    /// Multicenter bond index — maps to `RelationId` in the multicenter bonds relation set.
    MulticenterBondId,
    /// Noncovalent bond index — maps to `RelationId` in the noncovalent bonds relation set.
    NoncovalentBondId,
);

impl From<NodeId> for AtomId {
    fn from(id: NodeId) -> Self {
        Self(id.0)
    }
}
impl From<AtomId> for NodeId {
    fn from(idx: AtomId) -> Self {
        Self(idx.0)
    }
}

impl From<EdgeId> for BondId {
    fn from(id: EdgeId) -> Self {
        Self(id.0)
    }
}
impl From<BondId> for EdgeId {
    fn from(idx: BondId) -> Self {
        Self(idx.0)
    }
}

macro_rules! relation_idx_from {
    ($name:ident) => {
        impl From<RelationId> for $name {
            fn from(id: RelationId) -> Self {
                Self(id.0)
            }
        }
        impl From<$name> for RelationId {
            fn from(idx: $name) -> Self {
                Self(idx.0)
            }
        }
    };
}

relation_idx_from!(DativeBondId);
relation_idx_from!(AromaticSystemId);
relation_idx_from!(MulticenterBondId);
relation_idx_from!(NoncovalentBondId);
