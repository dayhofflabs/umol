// GraphBond implementation

use super::atom::GraphAtom;
use crate::link::AtomLink;
use super::types::AtomIndex;

#[derive(Debug, Clone, Copy)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
}

#[derive(Debug, Clone)]
pub struct GraphBond {
    order: BondOrder,
    between: (AtomIndex, AtomIndex),
}

impl GraphBond {
    pub fn new(from: AtomIndex, to: AtomIndex, order: BondOrder) -> Self {
        GraphBond {
            order,
            between: if from < to { (from, to) } else { (to, from) },
        }
    }

    pub fn single(from: AtomIndex, to: AtomIndex) -> Self {
        Self::new(from, to, BondOrder::Single)
    }

    pub fn double(from: AtomIndex, to: AtomIndex) -> Self {
        Self::new(from, to, BondOrder::Double)
    }

    pub fn triple(from: AtomIndex, to: AtomIndex) -> Self {
        Self::new(from, to, BondOrder::Triple)
    }

    pub fn order(&self) -> BondOrder {
        self.order
    }
}

impl AtomLink<GraphAtom> for GraphBond {
    type SiteRef = AtomIndex;
    fn between(&self) -> (Self::SiteRef, Self::SiteRef) {
        self.between
    }
}
