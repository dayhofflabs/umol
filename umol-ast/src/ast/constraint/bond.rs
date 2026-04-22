//! Per-covalent-bond constraints.

use strum::EnumDiscriminants;

use crate::ast::value::ValueAst;

#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(BondConstraintKind), derive(Hash))]
pub enum BondConstraint {
    Aromatic,
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl BondConstraint {
    pub fn kind(&self) -> BondConstraintKind {
        self.into()
    }
}
