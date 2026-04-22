//! Per-aromatic-system constraints.

use super::atom::AtomConstraint;
use crate::ast::idx::AtomIdx;
use crate::ast::remap::IdxRemapping;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticSystemConstraint {
    Atoms(Vec<AtomIdx>),
    Contains(AtomIdx),
    ContainsAll(Vec<AtomIdx>),
    AllAtoms(Box<AtomConstraint>),
    AnyAtom(Box<AtomConstraint>),
}

impl AromaticSystemConstraint {
    pub fn remap(self, remap: &IdxRemapping) -> Option<Self> {
        match self {
            Self::Atoms(atoms) => {
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                atoms.map(Self::Atoms)
            }
            Self::Contains(a) => remap.atom(a).map(Self::Contains),
            Self::ContainsAll(atoms) => {
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                atoms.map(Self::ContainsAll)
            }
            Self::AllAtoms(c) => Some(Self::AllAtoms(c)),
            Self::AnyAtom(c) => Some(Self::AnyAtom(c)),
        }
    }
}
