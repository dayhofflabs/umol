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

use error::LoweringError;
use index_vec::Idx;
use umol_edn::{FromEdn, ToEdn};

macro_rules! define_idx {
    ($($(#[doc = $doc:literal])* $name:ident),+ $(,)?) => {
        $(
            $(#[doc = $doc])*
            #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, FromEdn, ToEdn)]
            #[edn(transparent)]
            pub struct $name(pub usize);

            impl From<usize> for $name {
                fn from(v: usize) -> Self { Self(v) }
            }

            impl Idx for $name {
                fn from_usize(v: usize) -> Self { Self(v) }
                fn index(self) -> usize { self.0 }
            }
        )+
    };
}

define_idx!(
    /// Index into `MoleculeAst::atoms`.
    AtomIdx,
    /// Index into `MoleculeAst::bonds`.
    BondIdx,
    /// Index into `MoleculeAst::dative_bonds`.
    DativeBondIdx,
    /// Index into `MoleculeAst::aromatic_systems`.
    AromaticSystemIdx,
    /// Index into `MoleculeAst::multicenter_bonds`.
    MulticenterBondIdx,
    /// Index into `MoleculeAst::noncovalent_bonds`.
    NoncovalentBondIdx,
);

/// Construct an `IndexVec<AtomIdx, AtomAst>` from a list of atoms.
#[macro_export]
macro_rules! atoms {
    ($($atom:expr),* $(,)?) => {
        ::index_vec::IndexVec::<$crate::ast::AtomIdx, _>::from(vec![$($atom),*])
    };
}

/// AST marker trait carrying the lowering/raising config type.
pub trait Ast {
    type Config;
}

/// Lowering targets that can be built from an AST.
pub trait FromAst<A: Ast>: Sized {
    fn from_ast(ast: &A, cfg: &A::Config) -> Result<Self, LoweringError>;
}

/// Raise a ground or pattern type to its AST representation.
pub trait ToAst<A: Ast> {
    fn to_ast(&self, cfg: &A::Config) -> A;
}
