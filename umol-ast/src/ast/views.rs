//! Read-only views over `MoleculeAst` topology and relations.
//!
//! View records bundle an index with the underlying data so consumers
//! never assemble (idx, data, participants) tuples by hand. Namespace
//! types group per-relation accessors (`count`, `ids`, `iter`, `get`,
//! and `Index`) without burying them on `MoleculeAst` itself.

use std::ops::Index;

use umol_graph_core::relation::RelationId;
use umol_graph_core::{EdgeId, FixedRelationSet, Graph, NodeId, VarRelationSet};

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::dative::{DativeBondAst, DativeBondDirection};
use super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;

/// Borrowed view of an atom: its index and the underlying `AtomAst`.
#[derive(Clone, Copy, Debug)]
pub struct AtomView<'a> {
    pub idx: AtomIdx,
    pub data: &'a AtomAst,
}

/// Mutable borrowed view of an atom.
#[derive(Debug)]
pub struct AtomViewMut<'a> {
    pub idx: AtomIdx,
    pub data: &'a mut AtomAst,
}

/// Borrowed view of a bond: its index, endpoints (`src`, `tgt`), and data.
#[derive(Clone, Copy, Debug)]
pub struct BondView<'a> {
    pub idx: BondIdx,
    pub src: AtomIdx,
    pub tgt: AtomIdx,
    pub data: &'a BondAst,
}

/// Mutable borrowed view of a bond.
#[derive(Debug)]
pub struct BondViewMut<'a> {
    pub idx: BondIdx,
    pub src: AtomIdx,
    pub tgt: AtomIdx,
    pub data: &'a mut BondAst,
}

/// Neighbor-side view of a bond: the atom on the other end (`atom`), the
/// bond index, and the bond data. Yielded by `MoleculeAst::neighbors`.
#[derive(Clone, Copy, Debug)]
pub struct NeighborView<'a> {
    pub atom: AtomIdx,
    pub bond: BondIdx,
    pub data: &'a BondAst,
}

/// Borrowed view of a dative bond: donor and acceptor atom indices plus data.
#[derive(Clone, Copy, Debug)]
pub struct DativeBondView<'a> {
    pub idx: DativeBondIdx,
    pub donor: AtomIdx,
    pub acceptor: AtomIdx,
    pub data: &'a DativeBondAst,
}

/// Borrowed view of a noncovalent bond: the two participating atoms plus data.
#[derive(Clone, Copy, Debug)]
pub struct NoncovalentBondView<'a> {
    pub idx: NoncovalentBondIdx,
    pub atoms: [AtomIdx; 2],
    pub data: &'a NoncovalentBondAst,
}

/// Borrowed view of an aromatic system: its index, the `AromaticSystemAst`,
/// and accessors for member atoms and induced ring bonds via `atoms()` and
/// `bonds()`.
#[derive(Clone, Copy, Debug)]
pub struct AromaticSystemView<'a> {
    pub idx: AromaticSystemIdx,
    pub data: &'a AromaticSystemAst,
    atoms: &'a [NodeId],
    graph: &'a Graph,
}

impl<'a> AromaticSystemView<'a> {
    pub fn atoms(&self) -> impl Iterator<Item = AtomIdx> + '_ {
        self.atoms.iter().map(|&n| AtomIdx::from(n))
    }

    pub fn bonds(&self) -> impl Iterator<Item = BondIdx> + '_ {
        self.graph.induced_edges(self.atoms).map(BondIdx::from)
    }
}

/// Borrowed view of a multicenter bond: its index, data, and member atoms
/// via `atoms()`.
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

/// Namespace accessor for atom views on a `MoleculeAst`. Provides `count`,
/// `ids`, `iter`, `get`, and `Index` without burying them on `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct AtomViews<'a> {
    atoms: &'a [AtomAst],
}

impl<'a> AtomViews<'a> {
    pub(super) fn new(atoms: &'a [AtomAst]) -> Self {
        Self { atoms }
    }

    pub fn count(&self) -> usize {
        self.atoms.len()
    }

    pub fn ids(&self) -> impl Iterator<Item = AtomIdx> {
        (0..self.atoms.len() as u32).map(AtomIdx)
    }

    pub fn iter(&self) -> impl Iterator<Item = AtomView<'a>> {
        self.atoms.iter().enumerate().map(|(i, data)| AtomView {
            idx: AtomIdx(i as u32),
            data,
        })
    }

    pub fn get(&self, idx: AtomIdx) -> AtomView<'a> {
        AtomView {
            idx,
            data: &self.atoms[idx.index()],
        }
    }
}

impl<'a> Index<AtomIdx> for AtomViews<'a> {
    type Output = AtomAst;
    fn index(&self, idx: AtomIdx) -> &AtomAst {
        &self.atoms[idx.index()]
    }
}

/// Namespace accessor for bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct BondViews<'a> {
    bonds: &'a [BondAst],
    graph: &'a Graph,
}

impl<'a> BondViews<'a> {
    pub(super) fn new(bonds: &'a [BondAst], graph: &'a Graph) -> Self {
        Self { bonds, graph }
    }

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

/// Namespace accessor for dative-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct DativeBondViews<'a> {
    set: &'a FixedRelationSet<DativeBondAst, 2>,
}

impl<'a> DativeBondViews<'a> {
    pub(super) fn new(set: &'a FixedRelationSet<DativeBondAst, 2>) -> Self {
        Self { set }
    }

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
            let data = set.data(rid);
            let (donor, acceptor) = directed_endpoints(parts, data.direction);
            DativeBondView {
                idx: DativeBondIdx::from(rid),
                donor,
                acceptor,
                data,
            }
        })
    }

    pub fn get(&self, idx: DativeBondIdx) -> DativeBondView<'a> {
        let rid = RelationId::from(idx);
        let parts = self.set.participants(rid);
        let data = self.set.data(rid);
        let (donor, acceptor) = directed_endpoints(parts, data.direction);
        DativeBondView {
            idx,
            donor,
            acceptor,
            data,
        }
    }
}

fn directed_endpoints(parts: &[NodeId; 2], direction: DativeBondDirection) -> (AtomIdx, AtomIdx) {
    match direction {
        DativeBondDirection::Forward => (AtomIdx::from(parts[0]), AtomIdx::from(parts[1])),
        DativeBondDirection::Reverse => (AtomIdx::from(parts[1]), AtomIdx::from(parts[0])),
    }
}

impl<'a> Index<DativeBondIdx> for DativeBondViews<'a> {
    type Output = DativeBondAst;
    fn index(&self, idx: DativeBondIdx) -> &DativeBondAst {
        self.set.data(RelationId::from(idx))
    }
}

/// Namespace accessor for noncovalent-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct NoncovalentBondViews<'a> {
    set: &'a FixedRelationSet<NoncovalentBondAst, 2>,
}

impl<'a> NoncovalentBondViews<'a> {
    pub(super) fn new(set: &'a FixedRelationSet<NoncovalentBondAst, 2>) -> Self {
        Self { set }
    }

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
    type Output = NoncovalentBondAst;
    fn index(&self, idx: NoncovalentBondIdx) -> &NoncovalentBondAst {
        self.set.data(RelationId::from(idx))
    }
}

/// Namespace accessor for aromatic-system views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct AromaticSystemViews<'a> {
    set: &'a VarRelationSet<AromaticSystemAst>,
    graph: &'a Graph,
}

impl<'a> AromaticSystemViews<'a> {
    pub(super) fn new(set: &'a VarRelationSet<AromaticSystemAst>, graph: &'a Graph) -> Self {
        Self { set, graph }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = AromaticSystemIdx> {
        self.set.relation_ids().map(AromaticSystemIdx::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = AromaticSystemView<'a>> {
        let set = self.set;
        let graph = self.graph;
        set.relation_ids().map(move |rid| AromaticSystemView {
            idx: AromaticSystemIdx::from(rid),
            data: set.data(rid),
            atoms: set.participants(rid),
            graph,
        })
    }

    pub fn get(&self, idx: AromaticSystemIdx) -> AromaticSystemView<'a> {
        let rid = RelationId::from(idx);
        AromaticSystemView {
            idx,
            data: self.set.data(rid),
            atoms: self.set.participants(rid),
            graph: self.graph,
        }
    }
}

impl<'a> Index<AromaticSystemIdx> for AromaticSystemViews<'a> {
    type Output = AromaticSystemAst;
    fn index(&self, idx: AromaticSystemIdx) -> &AromaticSystemAst {
        self.set.data(RelationId::from(idx))
    }
}

/// Namespace accessor for multicenter-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct MulticenterBondViews<'a> {
    set: &'a VarRelationSet<MulticenterBondAst>,
}

impl<'a> MulticenterBondViews<'a> {
    pub(super) fn new(set: &'a VarRelationSet<MulticenterBondAst>) -> Self {
        Self { set }
    }

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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::constraint::Constraints;
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};

    #[fixture]
    fn rich() -> MoleculeAst {
        MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            vec![
                (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
                (AtomIdx(1), AtomIdx(2), BondAst::from_order(2)),
                (AtomIdx(2), AtomIdx(3), BondAst::from_order(1)),
            ],
            vec![(AtomIdx(2), AtomIdx(3), DativeBondAst::new())],
            vec![(
                vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
                AromaticSystemAst::default(),
            )],
            vec![(
                vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
                MulticenterBondAst::default(),
            )],
            vec![(
                AtomIdx(0),
                AtomIdx(3),
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_atom_views_count_and_ids(rich: MoleculeAst) {
        let views = rich.atoms();
        assert_eq!(views.count(), 4);
        assert_eq!(
            views.ids().collect::<Vec<_>>(),
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)],
        );
    }

    #[rstest]
    fn test_atom_views_index_trait(rich: MoleculeAst) {
        let views = rich.atoms();
        let atom: &AtomAst = &views[AtomIdx(2)];
        assert_eq!(*atom, AtomAst::from_element(Element::N));
    }

    #[rstest]
    fn test_bond_views_count_and_ids(rich: MoleculeAst) {
        let views = rich.bonds();
        assert_eq!(views.count(), 3);
        assert_eq!(
            views.ids().collect::<Vec<_>>(),
            vec![BondIdx(0), BondIdx(1), BondIdx(2)],
        );
    }

    #[rstest]
    fn test_bond_views_index_trait(rich: MoleculeAst) {
        let views = rich.bonds();
        let bond: &BondAst = &views[BondIdx(1)];
        assert_eq!(*bond, BondAst::from_order(2));
    }

    #[rstest]
    fn test_dative_bond_views_count_ids_and_index(rich: MoleculeAst) {
        let views = rich.dative_bonds();
        assert_eq!(views.count(), 1);
        assert_eq!(views.ids().collect::<Vec<_>>(), vec![DativeBondIdx(0)]);
        let _: &DativeBondAst = &views[DativeBondIdx(0)];
    }

    #[rstest]
    fn test_aromatic_system_views_count_ids_and_index(rich: MoleculeAst) {
        let views = rich.aromatic_systems();
        assert_eq!(views.count(), 1);
        assert_eq!(views.ids().collect::<Vec<_>>(), vec![AromaticSystemIdx(0)],);
        let _: &AromaticSystemAst = &views[AromaticSystemIdx(0)];
    }

    #[rstest]
    fn test_multicenter_bond_views_count_ids_and_index(rich: MoleculeAst) {
        let views = rich.multicenter_bonds();
        assert_eq!(views.count(), 1);
        assert_eq!(views.ids().collect::<Vec<_>>(), vec![MulticenterBondIdx(0)],);
        let _: &MulticenterBondAst = &views[MulticenterBondIdx(0)];
    }

    #[rstest]
    fn test_noncovalent_bond_views_count_ids_and_index(rich: MoleculeAst) {
        let views = rich.noncovalent_bonds();
        assert_eq!(views.count(), 1);
        assert_eq!(views.ids().collect::<Vec<_>>(), vec![NoncovalentBondIdx(0)],);
        let _: &NoncovalentBondAst = &views[NoncovalentBondIdx(0)];
    }
}
