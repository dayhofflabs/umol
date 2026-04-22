//! Per-noncovalent-bond constraints.

use super::atom::AtomConstraint;
use crate::ast::idx::AtomIdx;
use crate::ast::remap::IdxRemapping;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondConstraint {
    Ends([AtomIdx; 2]),
    Contains(AtomIdx),
    EndsSatisfy([Box<AtomConstraint>; 2]),
}

impl NoncovalentBondConstraint {
    pub fn remap(self, remap: &IdxRemapping) -> Option<Self> {
        match self {
            Self::Ends([a, b]) => {
                let a = remap.atom(a)?;
                let b = remap.atom(b)?;
                Some(Self::Ends([a, b]))
            }
            Self::Contains(a) => remap.atom(a).map(Self::Contains),
            Self::EndsSatisfy(cs) => Some(Self::EndsSatisfy(cs)),
        }
    }
}
