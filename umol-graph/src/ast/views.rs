//! Read-only views over `MoleculeAst` topology and relations.
//!
//! View records bundle an index with the underlying data so consumers
//! never assemble (idx, data, participants) tuples by hand. Namespace
//! types group per-relation accessors (`count`, `ids`, `iter`, `get`,
//! and `Index`) without burying them on `MoleculeAst` itself.

use std::ops::Index;

use umol_graph_core::relation::RelationId;
use umol_graph_core::{EdgeId, FixedRelationSet, Graph, NodeId, VarRelationSet};

use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::molecule::{AromaticSystemAst, MulticenterBondAst};
use crate::ast::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};

#[derive(Clone, Copy, Debug)]
pub struct AtomView<'a> {
    pub idx: AtomIdx,
    pub data: &'a AtomAst,
}

#[derive(Debug)]
pub struct AtomViewMut<'a> {
    pub idx: AtomIdx,
    pub data: &'a mut AtomAst,
}

#[derive(Clone, Copy, Debug)]
pub struct BondView<'a> {
    pub idx: BondIdx,
    pub src: AtomIdx,
    pub tgt: AtomIdx,
    pub data: &'a BondAst,
}

#[derive(Debug)]
pub struct BondViewMut<'a> {
    pub idx: BondIdx,
    pub src: AtomIdx,
    pub tgt: AtomIdx,
    pub data: &'a mut BondAst,
}

#[derive(Clone, Copy, Debug)]
pub struct NeighborView<'a> {
    pub atom: AtomIdx,
    pub bond: BondIdx,
    pub data: &'a BondAst,
}

#[derive(Clone, Copy, Debug)]
pub struct DativeBondView<'a> {
    pub idx: DativeBondIdx,
    pub donor: AtomIdx,
    pub acceptor: AtomIdx,
    pub data: &'a BondAst,
}

#[derive(Clone, Copy, Debug)]
pub struct NoncovalentBondView<'a> {
    pub idx: NoncovalentBondIdx,
    pub atoms: [AtomIdx; 2],
    pub data: &'a BondAst,
}

#[derive(Clone, Copy, Debug)]
pub struct AromaticSystemView<'a> {
    pub idx: AromaticSystemIdx,
    pub data: &'a AromaticSystemAst,
    atoms: &'a [NodeId],
}

impl<'a> AromaticSystemView<'a> {
    pub fn atoms(&self) -> impl Iterator<Item = AtomIdx> + '_ {
        self.atoms.iter().map(|&n| AtomIdx::from(n))
    }
}

#[derive(Clone, Copy, Debug)]
pub struct MulticenterBondView<'a> {
    pub idx: MulticenterBondIdx,
    pub data: &'a MulticenterBondAst,
    atoms: &'a [NodeId],
}

impl<'a> MulticenterBondView<'a> {
    pub fn atoms(&self) -> impl Iterator<Item = AtomIdx> + '_ {
        self.atoms.iter().map(|&n| AtomIdx::from(n))
    }
}

#[derive(Clone, Copy)]
pub struct AtomViews<'a> {
    pub(crate) atoms: &'a [AtomAst],
}

impl<'a> AtomViews<'a> {
    pub fn count(&self) -> usize {
        self.atoms.len()
    }

    pub fn ids(&self) -> impl Iterator<Item = AtomIdx> {
        (0..self.atoms.len() as u32).map(AtomIdx)
    }

    pub fn iter(&self) -> impl Iterator<Item = AtomView<'a>> {
        self.atoms
            .iter()
            .enumerate()
            .map(|(i, data)| AtomView { idx: AtomIdx(i as u32), data })
    }

    pub fn get(&self, idx: AtomIdx) -> AtomView<'a> {
        AtomView { idx, data: &self.atoms[idx.index()] }
    }
}

impl<'a> Index<AtomIdx> for AtomViews<'a> {
    type Output = AtomAst;
    fn index(&self, idx: AtomIdx) -> &AtomAst {
        &self.atoms[idx.index()]
    }
}

#[derive(Clone, Copy)]
pub struct BondViews<'a> {
    pub(crate) bonds: &'a [BondAst],
    pub(crate) graph: &'a Graph,
}

impl<'a> BondViews<'a> {
    pub fn count(&self) -> usize {
        self.bonds.len()
    }

    pub fn ids(&self) -> impl Iterator<Item = BondIdx> {
        (0..self.bonds.len() as u32).map(BondIdx)
    }

    pub fn iter(&self) -> impl Iterator<Item = BondView<'a>> {
        let bonds = self.bonds;
        let graph = self.graph;
        graph.edge_ids().map(move |id| {
            let [s, t] = graph.edge_endpoints(id);
            BondView {
                idx: BondIdx::from(id),
                src: AtomIdx::from(s),
                tgt: AtomIdx::from(t),
                data: &bonds[id.index()],
            }
        })
    }

    pub fn get(&self, idx: BondIdx) -> BondView<'a> {
        let [s, t] = self.graph.edge_endpoints(EdgeId::from(idx));
        BondView {
            idx,
            src: AtomIdx::from(s),
            tgt: AtomIdx::from(t),
            data: &self.bonds[idx.index()],
        }
    }
}

impl<'a> Index<BondIdx> for BondViews<'a> {
    type Output = BondAst;
    fn index(&self, idx: BondIdx) -> &BondAst {
        &self.bonds[idx.index()]
    }
}

#[derive(Clone, Copy)]
pub struct DativeBondViews<'a> {
    pub(crate) set: &'a FixedRelationSet<BondAst, 2>,
}

impl<'a> DativeBondViews<'a> {
    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = DativeBondIdx> {
        self.set.relation_ids().map(DativeBondIdx::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = DativeBondView<'a>> {
        let set = self.set;
        set.relation_ids().map(move |rid| {
            let parts = set.participants(rid);
            DativeBondView {
                idx: DativeBondIdx::from(rid),
                donor: AtomIdx::from(parts[0]),
                acceptor: AtomIdx::from(parts[1]),
                data: set.data(rid),
            }
        })
    }

    pub fn get(&self, idx: DativeBondIdx) -> DativeBondView<'a> {
        let rid = RelationId::from(idx);
        let parts = self.set.participants(rid);
        DativeBondView {
            idx,
            donor: AtomIdx::from(parts[0]),
            acceptor: AtomIdx::from(parts[1]),
            data: self.set.data(rid),
        }
    }
}

impl<'a> Index<DativeBondIdx> for DativeBondViews<'a> {
    type Output = BondAst;
    fn index(&self, idx: DativeBondIdx) -> &BondAst {
        self.set.data(RelationId::from(idx))
    }
}

#[derive(Clone, Copy)]
pub struct NoncovalentBondViews<'a> {
    pub(crate) set: &'a FixedRelationSet<BondAst, 2>,
}

impl<'a> NoncovalentBondViews<'a> {
    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = NoncovalentBondIdx> {
        self.set.relation_ids().map(NoncovalentBondIdx::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = NoncovalentBondView<'a>> {
        let set = self.set;
        set.relation_ids().map(move |rid| {
            let parts = set.participants(rid);
            NoncovalentBondView {
                idx: NoncovalentBondIdx::from(rid),
                atoms: [AtomIdx::from(parts[0]), AtomIdx::from(parts[1])],
                data: set.data(rid),
            }
        })
    }

    pub fn get(&self, idx: NoncovalentBondIdx) -> NoncovalentBondView<'a> {
        let rid = RelationId::from(idx);
        let parts = self.set.participants(rid);
        NoncovalentBondView {
            idx,
            atoms: [AtomIdx::from(parts[0]), AtomIdx::from(parts[1])],
            data: self.set.data(rid),
        }
    }
}

impl<'a> Index<NoncovalentBondIdx> for NoncovalentBondViews<'a> {
    type Output = BondAst;
    fn index(&self, idx: NoncovalentBondIdx) -> &BondAst {
        self.set.data(RelationId::from(idx))
    }
}

#[derive(Clone, Copy)]
pub struct AromaticSystemViews<'a> {
    pub(crate) set: &'a VarRelationSet<AromaticSystemAst>,
}

impl<'a> AromaticSystemViews<'a> {
    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = AromaticSystemIdx> {
        self.set.relation_ids().map(AromaticSystemIdx::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = AromaticSystemView<'a>> {
        let set = self.set;
        set.relation_ids().map(move |rid| AromaticSystemView {
            idx: AromaticSystemIdx::from(rid),
            data: set.data(rid),
            atoms: set.participants(rid),
        })
    }

    pub fn get(&self, idx: AromaticSystemIdx) -> AromaticSystemView<'a> {
        let rid = RelationId::from(idx);
        AromaticSystemView {
            idx,
            data: self.set.data(rid),
            atoms: self.set.participants(rid),
        }
    }
}

impl<'a> Index<AromaticSystemIdx> for AromaticSystemViews<'a> {
    type Output = AromaticSystemAst;
    fn index(&self, idx: AromaticSystemIdx) -> &AromaticSystemAst {
        self.set.data(RelationId::from(idx))
    }
}

#[derive(Clone, Copy)]
pub struct MulticenterBondViews<'a> {
    pub(crate) set: &'a VarRelationSet<MulticenterBondAst>,
}

impl<'a> MulticenterBondViews<'a> {
    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = MulticenterBondIdx> {
        self.set.relation_ids().map(MulticenterBondIdx::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = MulticenterBondView<'a>> {
        let set = self.set;
        set.relation_ids().map(move |rid| MulticenterBondView {
            idx: MulticenterBondIdx::from(rid),
            data: set.data(rid),
            atoms: set.participants(rid),
        })
    }

    pub fn get(&self, idx: MulticenterBondIdx) -> MulticenterBondView<'a> {
        let rid = RelationId::from(idx);
        MulticenterBondView {
            idx,
            data: self.set.data(rid),
            atoms: self.set.participants(rid),
        }
    }
}

impl<'a> Index<MulticenterBondIdx> for MulticenterBondViews<'a> {
    type Output = MulticenterBondAst;
    fn index(&self, idx: MulticenterBondIdx) -> &MulticenterBondAst {
        self.set.data(RelationId::from(idx))
    }
}
