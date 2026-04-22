//! Per-dative-bond constraints.

use strum::EnumDiscriminants;

use super::atom::AtomConstraint;
use crate::ast::idx::{AtomIdx, BondIdx};
use crate::ast::remap::IdxRemapping;
use crate::ast::value::ValueAst;

#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(DativeBondConstraintKind), derive(Hash))]
pub enum DativeBondConstraint {
    RingCount(ValueAst),
    RingSize(ValueAst),
    Donor(AtomIdx),
    Acceptor(AtomIdx),
    DonorSatisfies(Box<AtomConstraint>),
    AcceptorSatisfies(Box<AtomConstraint>),
    Parallels(BondIdx),
}

impl DativeBondConstraint {
    pub fn kind(&self) -> DativeBondConstraintKind {
        self.into()
    }

    pub fn remap(self, remap: &IdxRemapping) -> Option<Self> {
        match self {
            Self::RingCount(v) => Some(Self::RingCount(v)),
            Self::RingSize(v) => Some(Self::RingSize(v)),
            Self::Donor(a) => remap.atom(a).map(Self::Donor),
            Self::Acceptor(a) => remap.atom(a).map(Self::Acceptor),
            Self::DonorSatisfies(c) => Some(Self::DonorSatisfies(c)),
            Self::AcceptorSatisfies(c) => Some(Self::AcceptorSatisfies(c)),
            Self::Parallels(b) => remap.bond(b).map(Self::Parallels),
        }
    }
}
