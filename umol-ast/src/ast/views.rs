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
use super::constraint::{
    AromaticValenceAst, AtomConstraint, AtomConstraintKind, MulticenterValenceAst,
};
use super::dative::{DativeBondAst, DativeBondDirection};
use super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use super::molecule::MoleculeAst;
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;
use super::value::ValueAst;

/// Borrowed view of an atom: index, underlying `AtomAst`, and the parent
/// `MoleculeAst` for cross-relation chemistry methods.
///
/// Chemistry methods come in pairs: the topology-derived value (summed from
/// incident bonds / dative bonds / aromatic system / multicenter bonds) and
/// the matching local-constraint value carried in `data.constraints`. The
/// validator cross-checks the two when both are ground.
#[derive(Clone, Copy, Debug)]
pub struct AtomView<'a> {
    pub idx: AtomIdx,
    pub data: &'a AtomAst,
    ast: &'a MoleculeAst,
}

impl<'a> AtomView<'a> {
    /// σ-valence summed from incident bond orders. `None` if any incident
    /// bond's order is not a non-negative `Lit`.
    pub fn bond_order_sum(&self) -> Option<u32> {
        let mut sum: u32 = 0;
        for n in self.ast.neighbors(self.idx) {
            match n.data.order {
                ValueAst::Lit(v) if v >= 0 => sum = sum.checked_add(v as u32)?,
                _ => return None,
            }
        }
        Some(sum)
    }

    /// Local valence constraint, if asserted.
    pub fn valence_constraint(&self) -> Option<&'a ValueAst> {
        atom_constraint_value(self.data, AtomConstraintKind::Valence, |c| match c {
            AtomConstraint::Valence(v) => Some(v),
            _ => None,
        })
    }

    /// Number of incident dative bonds where this atom is the donor.
    pub fn donated_pairs(&self) -> u32 {
        self.ast
            .dative_bonds_incident(self.idx)
            .filter(|&id| self.ast.dative_bond(id).donor == self.idx)
            .count() as u32
    }

    pub fn donated_pairs_constraint(&self) -> Option<&'a ValueAst> {
        atom_constraint_value(self.data, AtomConstraintKind::DonatedPairs, |c| match c {
            AtomConstraint::DonatedPairs(v) => Some(v),
            _ => None,
        })
    }

    /// Number of incident dative bonds where this atom is the acceptor.
    pub fn accepted_pairs(&self) -> u32 {
        self.ast
            .dative_bonds_incident(self.idx)
            .filter(|&id| self.ast.dative_bond(id).acceptor == self.idx)
            .count() as u32
    }

    pub fn accepted_pairs_constraint(&self) -> Option<&'a ValueAst> {
        atom_constraint_value(self.data, AtomConstraintKind::AcceptedPairs, |c| match c {
            AtomConstraint::AcceptedPairs(v) => Some(v),
            _ => None,
        })
    }

    /// π contribution from the aromatic system this atom belongs to.
    /// `Some(0)` if the atom is not in any aromatic system. `None` if the
    /// recorded contribution is not a non-negative `Lit`.
    ///
    /// An atom belongs to at most one aromatic system; the first incident
    /// system is consulted.
    pub fn aromatic_contribution(&self) -> Option<u32> {
        let Some(sys_id) = self.ast.aromatic_systems_incident(self.idx).next() else {
            return Some(0);
        };
        let view = self.ast.aromatic_system(sys_id);
        let pos = view.atoms().position(|a| a == self.idx)?;
        match view.data.electrons.get(pos)? {
            ValueAst::Lit(v) if *v >= 0 => Some(*v as u32),
            _ => None,
        }
    }

    pub fn is_in_aromatic_system(&self) -> bool {
        self.ast.aromatic_systems_incident(self.idx).next().is_some()
    }

    pub fn aromatic_valence_constraint(&self) -> Option<&'a AromaticValenceAst> {
        atom_constraint_value(self.data, AtomConstraintKind::AromaticValence, |c| match c {
            AtomConstraint::AromaticValence(v) => Some(v),
            _ => None,
        })
    }

    /// Sum of per-atom contributions across incident multicenter bonds.
    /// `None` if any contribution is not a non-negative `Lit`.
    pub fn multicenter_contribution(&self) -> Option<u32> {
        let mut sum: u32 = 0;
        for mc_id in self.ast.multicenter_bonds_incident(self.idx) {
            let view = self.ast.multicenter_bond(mc_id);
            let pos = view.atoms().position(|a| a == self.idx)?;
            match view.data.electrons.get(pos)? {
                ValueAst::Lit(v) if *v >= 0 => sum = sum.checked_add(*v as u32)?,
                _ => return None,
            }
        }
        Some(sum)
    }

    pub fn multicenter_valence_constraint(&self) -> Option<&'a MulticenterValenceAst> {
        atom_constraint_value(
            self.data,
            AtomConstraintKind::MulticenterValence,
            |c| match c {
                AtomConstraint::MulticenterValence(v) => Some(v),
                _ => None,
            },
        )
    }
}

fn atom_constraint_value<'a, T>(
    atom: &'a AtomAst,
    kind: AtomConstraintKind,
    extract: impl FnOnce(&'a AtomConstraint) -> Option<&'a T>,
) -> Option<&'a T> {
    extract(atom.constraints.get(kind)?)
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
    ast: &'a MoleculeAst,
    atoms: &'a [AtomAst],
}

impl<'a> AtomViews<'a> {
    pub(super) fn new(ast: &'a MoleculeAst, atoms: &'a [AtomAst]) -> Self {
        Self { ast, atoms }
    }

    pub fn count(&self) -> usize {
        self.atoms.len()
    }

    pub fn ids(&self) -> impl Iterator<Item = AtomIdx> {
        (0..self.atoms.len() as u32).map(AtomIdx)
    }

    pub fn iter(&self) -> impl Iterator<Item = AtomView<'a>> {
        let ast = self.ast;
        self.atoms.iter().enumerate().map(move |(i, data)| AtomView {
            idx: AtomIdx(i as u32),
            data,
            ast,
        })
    }

    pub fn get(&self, idx: AtomIdx) -> AtomView<'a> {
        AtomView {
            idx,
            data: &self.atoms[idx.index()],
            ast: self.ast,
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

    use crate::ast::bond::BondAst;
    use crate::ast::constraint::AtomConstraint;
    use crate::ast::dative::{DativeBondAst, DativeBondDirection};
    use crate::ast::spin::SpinStateAst;
    use crate::ast::value::ValueAst;

    fn atom_with_constraints(element: Element, cs: Vec<AtomConstraint>) -> AtomAst {
        let mut atom = AtomAst::from_element(element);
        for c in cs {
            atom.constraints.add(c);
        }
        atom
    }

    fn dative_with_direction(direction: DativeBondDirection) -> DativeBondAst {
        DativeBondAst {
            direction,
            constraints: Default::default(),
        }
    }

    fn aromatic_with_electrons(electrons: Vec<ValueAst>) -> AromaticSystemAst {
        AromaticSystemAst::new(electrons, ValueAst::Lit(0), SpinStateAst::default())
    }

    fn multicenter_with_electrons(electrons: Vec<ValueAst>) -> MulticenterBondAst {
        MulticenterBondAst::new(electrons, ValueAst::Lit(0), SpinStateAst::default())
    }

    #[rstest]
    #[case::no_bonds(AtomIdx(3), Some(0))]
    #[case::single(AtomIdx(0), Some(1))]
    #[case::two_incident(AtomIdx(1), Some(3))]
    #[case::double(AtomIdx(2), Some(2))]
    fn test_atom_view_bond_order_sum_ground(
        #[case] center: AtomIdx,
        #[case] expected: Option<u32>,
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ];
        let bonds = vec![
            (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
            (AtomIdx(1), AtomIdx(2), BondAst::from_order(2)),
        ];
        let ast = MoleculeAst::new(
            atoms,
            bonds,
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(ast.atom(center).bond_order_sum(), expected);
    }

    #[rstest]
    fn test_atom_view_bond_order_sum_undetermined() {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ];
        let mut undetermined = BondAst::from_order(1);
        undetermined.order = ValueAst::Undetermined;
        let bonds = vec![(AtomIdx(0), AtomIdx(1), undetermined)];
        let ast = MoleculeAst::new(
            atoms,
            bonds,
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(ast.atom(AtomIdx(0)).bond_order_sum(), None);
    }

    #[rstest]
    #[case::donor_forward(AtomIdx(0), 1, 0)]
    #[case::acceptor_forward(AtomIdx(1), 0, 1)]
    fn test_atom_view_dative_pair_counts(
        #[case] atom: AtomIdx,
        #[case] expected_donated: u32,
        #[case] expected_accepted: u32,
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::C),
        ];
        let dative = vec![(
            AtomIdx(0),
            AtomIdx(1),
            dative_with_direction(DativeBondDirection::Forward),
        )];
        let ast = MoleculeAst::new(
            atoms,
            vec![],
            dative,
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        let view = ast.atom(atom);
        assert_eq!(view.donated_pairs(), expected_donated);
        assert_eq!(view.accepted_pairs(), expected_accepted);
    }

    #[rstest]
    #[case::lit(ValueAst::Lit(2), Some(2))]
    #[case::undetermined(ValueAst::Undetermined, None)]
    fn test_atom_view_aromatic_contribution(
        #[case] entry: ValueAst,
        #[case] expected: Option<u32>,
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ];
        let bonds = vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))];
        let aromatic = vec![(
            vec![AtomIdx(0), AtomIdx(1)],
            aromatic_with_electrons(vec![entry, ValueAst::Lit(1)]),
        )];
        let ast = MoleculeAst::new(
            atoms,
            bonds,
            vec![],
            aromatic,
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(ast.atom(AtomIdx(0)).aromatic_contribution(), expected);
    }

    #[rstest]
    fn test_atom_view_aromatic_contribution_not_in_system() {
        let atoms = vec![AtomAst::from_element(Element::C)];
        let ast = MoleculeAst::new(
            atoms,
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        let view = ast.atom(AtomIdx(0));
        assert_eq!(view.aromatic_contribution(), Some(0));
        assert!(!view.is_in_aromatic_system());
    }

    #[rstest]
    #[case::single_bond(vec![(vec![AtomIdx(0), AtomIdx(1)], vec![ValueAst::Lit(2), ValueAst::Lit(2)])], Some(2))]
    #[case::two_bonds(
        vec![
            (vec![AtomIdx(0), AtomIdx(1)], vec![ValueAst::Lit(2), ValueAst::Lit(2)]),
            (vec![AtomIdx(0), AtomIdx(2)], vec![ValueAst::Lit(1), ValueAst::Lit(1)]),
        ],
        Some(3),
    )]
    #[case::undetermined_aborts(
        vec![(vec![AtomIdx(0), AtomIdx(1)], vec![ValueAst::Undetermined, ValueAst::Lit(2)])],
        None,
    )]
    fn test_atom_view_multicenter_contribution(
        #[case] bonds: Vec<(Vec<AtomIdx>, Vec<ValueAst>)>,
        #[case] expected: Option<u32>,
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ];
        let multicenter: Vec<_> = bonds
            .into_iter()
            .map(|(parts, electrons)| (parts, multicenter_with_electrons(electrons)))
            .collect();
        let ast = MoleculeAst::new(
            atoms,
            vec![],
            vec![],
            vec![],
            multicenter,
            vec![],
            Constraints::default(),
        );
        assert_eq!(ast.atom(AtomIdx(0)).multicenter_contribution(), expected);
    }

    #[rstest]
    fn test_atom_view_constraint_accessors_present() {
        let atoms = vec![atom_with_constraints(
            Element::C,
            vec![
                AtomConstraint::Valence(ValueAst::Lit(4)),
                AtomConstraint::DonatedPairs(ValueAst::Lit(0)),
                AtomConstraint::AcceptedPairs(ValueAst::Lit(0)),
                AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
                AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(
                    ValueAst::Lit(2),
                )),
            ],
        )];
        let ast = MoleculeAst::new(
            atoms,
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        let view = ast.atom(AtomIdx(0));
        assert_eq!(view.valence_constraint(), Some(&ValueAst::Lit(4)));
        assert_eq!(view.donated_pairs_constraint(), Some(&ValueAst::Lit(0)));
        assert_eq!(view.accepted_pairs_constraint(), Some(&ValueAst::Lit(0)));
        assert_eq!(
            view.aromatic_valence_constraint(),
            Some(&AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
        );
        assert_eq!(
            view.multicenter_valence_constraint(),
            Some(&MulticenterValenceAst::Multicenter(ValueAst::Lit(2))),
        );
    }

    #[rstest]
    fn test_atom_view_constraint_accessors_absent() {
        let atoms = vec![AtomAst::from_element(Element::C)];
        let ast = MoleculeAst::new(
            atoms,
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        let view = ast.atom(AtomIdx(0));
        assert!(view.valence_constraint().is_none());
        assert!(view.donated_pairs_constraint().is_none());
        assert!(view.accepted_pairs_constraint().is_none());
        assert!(view.aromatic_valence_constraint().is_none());
        assert!(view.multicenter_valence_constraint().is_none());
    }
}
