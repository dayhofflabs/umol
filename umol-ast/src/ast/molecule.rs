//! Molecule structural AST.

mod builder;

use std::collections::HashSet;
use std::ops::Index;
use std::sync::Arc;

pub use builder::MoleculeBuilder;
use umol_graph_core::relation::RelationId;
use umol_graph_core::{FixedRelationSet, Graph, NodeId, VarRelationSet};

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::constraint::Constraints;
use super::dative::DativeBondAst;
use super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;
use super::views::{
    AromaticSystemView, AromaticSystemViews, AtomView, AtomViewMut, AtomViews, BondView,
    BondViewMut, BondViews, DativeBondView, DativeBondViews, MulticenterBondView,
    MulticenterBondViews, NeighborView, NoncovalentBondView, NoncovalentBondViews,
};

/// Molecule AST: structural representation of a molecule (ground or pattern).
///
/// Topology and per-atom/bond data are `Arc`-shared (copy-on-write). The AST
/// itself only allows attribute mutation (`atom_mut`, `bond_mut`); structural
/// edits go through `MoleculeBuilder` via [`MoleculeAst::edit`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeAst {
    graph: Graph,
    atoms: Arc<Vec<AtomAst>>,
    bonds: Arc<Vec<BondAst>>,
    dative_bonds: Arc<FixedRelationSet<DativeBondAst, 2>>,
    aromatic_systems: Arc<VarRelationSet<AromaticSystemAst>>,
    multicenter_bonds: Arc<VarRelationSet<MulticenterBondAst>>,
    noncovalent_bonds: Arc<FixedRelationSet<NoncovalentBondAst, 2>>,
    constraints: Constraints,
}

impl Default for MoleculeAst {
    fn default() -> Self {
        Self {
            graph: Graph::default(),
            atoms: Arc::new(Vec::new()),
            bonds: Arc::new(Vec::new()),
            dative_bonds: Arc::new(FixedRelationSet::default()),
            aromatic_systems: Arc::new(VarRelationSet::default()),
            multicenter_bonds: Arc::new(VarRelationSet::default()),
            noncovalent_bonds: Arc::new(FixedRelationSet::default()),
            constraints: Constraints::new(),
        }
    }
}

impl MoleculeAst {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        atoms: Vec<AtomAst>,
        bonds: Vec<(AtomIdx, AtomIdx, BondAst)>,
        dative: Vec<(AtomIdx, AtomIdx, DativeBondAst)>,
        aromatic: Vec<(Vec<AtomIdx>, AromaticSystemAst)>,
        multicenter: Vec<(Vec<AtomIdx>, MulticenterBondAst)>,
        noncovalent: Vec<(AtomIdx, AtomIdx, NoncovalentBondAst)>,
        constraints: Constraints,
    ) -> Self {
        let node_count = atoms.len();
        let edges: Vec<[u32; 2]> = bonds.iter().map(|(s, t, _)| [s.0, t.0]).collect();
        let bond_data: Vec<BondAst> = bonds.into_iter().map(|(_, _, d)| d).collect();
        let graph = Graph::new(node_count, &edges);

        let dative_bonds = FixedRelationSet::new(
            dative
                .into_iter()
                .map(|(a, b, d)| ([NodeId::from(a), NodeId::from(b)], d))
                .collect(),
        );

        let aromatic_systems = VarRelationSet::new(
            aromatic
                .into_iter()
                .map(|(atoms, d)| (atoms.into_iter().map(NodeId::from).collect(), d))
                .collect(),
        );

        let multicenter_bonds = VarRelationSet::new(
            multicenter
                .into_iter()
                .map(|(atoms, d)| (atoms.into_iter().map(NodeId::from).collect(), d))
                .collect(),
        );

        let noncovalent_bonds = FixedRelationSet::new(
            noncovalent
                .into_iter()
                .map(|(a, b, d)| ([NodeId::from(a), NodeId::from(b)], d))
                .collect(),
        );

        Self {
            graph,
            atoms: Arc::new(atoms),
            bonds: Arc::new(bond_data),
            dative_bonds: Arc::new(dative_bonds),
            aromatic_systems: Arc::new(aromatic_systems),
            multicenter_bonds: Arc::new(multicenter_bonds),
            noncovalent_bonds: Arc::new(noncovalent_bonds),
            constraints,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn from_arcs(
        graph: Graph,
        atoms: Arc<Vec<AtomAst>>,
        bonds: Arc<Vec<BondAst>>,
        dative_bonds: Arc<FixedRelationSet<DativeBondAst, 2>>,
        aromatic_systems: Arc<VarRelationSet<AromaticSystemAst>>,
        multicenter_bonds: Arc<VarRelationSet<MulticenterBondAst>>,
        noncovalent_bonds: Arc<FixedRelationSet<NoncovalentBondAst, 2>>,
        constraints: Constraints,
    ) -> Self {
        Self {
            graph,
            atoms,
            bonds,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            constraints,
        }
    }

    // -- Read: topology ---------------------------------------------------

    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    pub fn neighbors(&self, atom: AtomIdx) -> impl Iterator<Item = NeighborView<'_>> {
        let bonds = &self.bonds;
        self.graph
            .neighbors(NodeId::from(atom))
            .iter()
            .map(move |n| NeighborView {
                atom: AtomIdx::from(n.node),
                bond: BondIdx::from(n.edge),
                data: &bonds[n.edge.index()],
            })
    }

    // -- Read: atoms ------------------------------------------------------

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn atoms(&self) -> AtomViews<'_> {
        AtomViews::new(&self.atoms)
    }

    pub fn atom(&self, idx: AtomIdx) -> AtomView<'_> {
        self.atoms().get(idx)
    }

    // -- Read: bonds ------------------------------------------------------

    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    pub fn bonds(&self) -> BondViews<'_> {
        BondViews::new(&self.bonds, &self.graph)
    }

    pub fn bond(&self, idx: BondIdx) -> BondView<'_> {
        self.bonds().get(idx)
    }

    // -- Read: dative bonds -----------------------------------------------

    pub fn dative_bond_count(&self) -> usize {
        self.dative_bonds.relation_count()
    }

    pub fn dative_bonds(&self) -> DativeBondViews<'_> {
        DativeBondViews::new(&self.dative_bonds)
    }

    pub fn dative_bond(&self, idx: DativeBondIdx) -> DativeBondView<'_> {
        self.dative_bonds().get(idx)
    }

    // -- Read: aromatic systems -------------------------------------------

    pub fn aromatic_system_count(&self) -> usize {
        self.aromatic_systems.relation_count()
    }

    pub fn aromatic_systems(&self) -> AromaticSystemViews<'_> {
        AromaticSystemViews::new(&self.aromatic_systems, &self.graph)
    }

    pub fn aromatic_system(&self, idx: AromaticSystemIdx) -> AromaticSystemView<'_> {
        self.aromatic_systems().get(idx)
    }

    // -- Read: multicenter bonds ------------------------------------------

    pub fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.relation_count()
    }

    pub fn multicenter_bonds(&self) -> MulticenterBondViews<'_> {
        MulticenterBondViews::new(&self.multicenter_bonds)
    }

    pub fn multicenter_bond(&self, idx: MulticenterBondIdx) -> MulticenterBondView<'_> {
        self.multicenter_bonds().get(idx)
    }

    // -- Read: noncovalent bonds ------------------------------------------

    pub fn noncovalent_bond_count(&self) -> usize {
        self.noncovalent_bonds.relation_count()
    }

    pub fn noncovalent_bonds(&self) -> NoncovalentBondViews<'_> {
        NoncovalentBondViews::new(&self.noncovalent_bonds)
    }

    pub fn noncovalent_bond(&self, idx: NoncovalentBondIdx) -> NoncovalentBondView<'_> {
        self.noncovalent_bonds().get(idx)
    }

    // -- Read: incidence ----------------------------------------------------

    pub fn connecting_bond(&self, a: AtomIdx, b: AtomIdx) -> Option<BondIdx> {
        self.graph
            .find_edge(NodeId::from(a), NodeId::from(b))
            .map(BondIdx::from)
    }

    pub fn dative_bonds_incident(
        &self,
        atom: AtomIdx,
    ) -> impl Iterator<Item = DativeBondIdx> + '_ {
        self.dative_bonds
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| DativeBondIdx::from(rid))
    }

    pub fn aromatic_systems_incident(
        &self,
        atom: AtomIdx,
    ) -> impl Iterator<Item = AromaticSystemIdx> + '_ {
        self.aromatic_systems
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| AromaticSystemIdx::from(rid))
    }

    pub fn multicenter_bonds_incident(
        &self,
        atom: AtomIdx,
    ) -> impl Iterator<Item = MulticenterBondIdx> + '_ {
        self.multicenter_bonds
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| MulticenterBondIdx::from(rid))
    }

    pub fn noncovalent_bonds_incident(
        &self,
        atom: AtomIdx,
    ) -> impl Iterator<Item = NoncovalentBondIdx> + '_ {
        self.noncovalent_bonds
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| NoncovalentBondIdx::from(rid))
    }

    // -- Read: induced subsets ----------------------------------------------

    pub fn induced_bonds(&self, atoms: &[AtomIdx]) -> Vec<BondIdx> {
        let mut nodes: Vec<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        nodes.sort_unstable();
        self.graph.induced_edges(&nodes).map(BondIdx::from).collect()
    }

    pub fn induced_dative_bonds(&self, atoms: &[AtomIdx]) -> Vec<DativeBondIdx> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.dative_bonds
            .relation_ids()
            .filter(|&rid| {
                self.dative_bonds
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(DativeBondIdx::from)
            .collect()
    }

    pub fn induced_aromatic_systems(&self, atoms: &[AtomIdx]) -> Vec<AromaticSystemIdx> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.aromatic_systems
            .relation_ids()
            .filter(|&rid| {
                self.aromatic_systems
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(AromaticSystemIdx::from)
            .collect()
    }

    pub fn induced_multicenter_bonds(&self, atoms: &[AtomIdx]) -> Vec<MulticenterBondIdx> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.multicenter_bonds
            .relation_ids()
            .filter(|&rid| {
                self.multicenter_bonds
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(MulticenterBondIdx::from)
            .collect()
    }

    pub fn induced_noncovalent_bonds(&self, atoms: &[AtomIdx]) -> Vec<NoncovalentBondIdx> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.noncovalent_bonds
            .relation_ids()
            .filter(|&rid| {
                self.noncovalent_bonds
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(NoncovalentBondIdx::from)
            .collect()
    }

    pub fn is_ground(&self) -> bool {
        self.atoms.iter().all(|a| a.is_ground())
            && self.bonds.iter().all(|b| b.is_ground())
            && self
                .dative_bonds
                .relation_ids()
                .all(|id| self.dative_bonds.data(id).is_ground())
            && self
                .aromatic_systems
                .relation_ids()
                .all(|id| self.aromatic_systems.data(id).is_ground())
            && self
                .multicenter_bonds
                .relation_ids()
                .all(|id| self.multicenter_bonds.data(id).is_ground())
            && self
                .noncovalent_bonds
                .relation_ids()
                .all(|id| self.noncovalent_bonds.data(id).is_ground())
    }

    // -- Entity mutation: atoms -------------------------------------------

    pub fn atom_mut(&mut self, idx: AtomIdx) -> AtomViewMut<'_> {
        let data = &mut Arc::make_mut(&mut self.atoms)[idx.index()];
        AtomViewMut { idx, data }
    }

    pub fn atoms_mut(&mut self) -> impl Iterator<Item = &mut AtomAst> {
        Arc::make_mut(&mut self.atoms).iter_mut()
    }

    // -- Entity mutation: bonds -------------------------------------------

    pub fn bond_mut(&mut self, idx: BondIdx) -> BondViewMut<'_> {
        let [s, t] = self.graph.edge_endpoints(idx.into());
        let data = &mut Arc::make_mut(&mut self.bonds)[idx.index()];
        BondViewMut {
            idx,
            src: AtomIdx::from(s),
            tgt: AtomIdx::from(t),
            data,
        }
    }

    pub fn bonds_mut(&mut self) -> impl Iterator<Item = &mut BondAst> {
        Arc::make_mut(&mut self.bonds).iter_mut()
    }

    // -- Entity mutation: dative bonds ------------------------------------

    pub fn dative_bond_mut(&mut self, idx: DativeBondIdx) -> &mut DativeBondAst {
        Arc::make_mut(&mut self.dative_bonds).data_mut(RelationId::from(idx))
    }

    // -- Entity mutation: aromatic systems --------------------------------

    pub fn aromatic_system_mut(&mut self, idx: AromaticSystemIdx) -> &mut AromaticSystemAst {
        Arc::make_mut(&mut self.aromatic_systems).data_mut(RelationId::from(idx))
    }

    pub fn aromatic_systems_mut(&mut self) -> impl Iterator<Item = &mut AromaticSystemAst> {
        Arc::make_mut(&mut self.aromatic_systems).data_iter_mut()
    }

    // -- Entity mutation: multicenter bonds -------------------------------

    pub fn multicenter_bond_mut(&mut self, idx: MulticenterBondIdx) -> &mut MulticenterBondAst {
        Arc::make_mut(&mut self.multicenter_bonds).data_mut(RelationId::from(idx))
    }

    pub fn multicenter_bonds_mut(&mut self) -> impl Iterator<Item = &mut MulticenterBondAst> {
        Arc::make_mut(&mut self.multicenter_bonds).data_iter_mut()
    }

    // -- Entity mutation: noncovalent bonds -------------------------------

    pub fn noncovalent_bond_mut(&mut self, idx: NoncovalentBondIdx) -> &mut NoncovalentBondAst {
        Arc::make_mut(&mut self.noncovalent_bonds).data_mut(RelationId::from(idx))
    }

    // -- Constraints ------------------------------------------------------

    pub fn constraints(&self) -> &Constraints {
        &self.constraints
    }

    pub fn constraints_mut(&mut self) -> &mut Constraints {
        &mut self.constraints
    }

    // -- Topological mutation (via builder) --------------------------------

    pub fn edit(&self) -> MoleculeBuilder {
        MoleculeBuilder::from_parts(
            self.graph.clone(),
            Arc::clone(&self.atoms),
            Arc::clone(&self.bonds),
            Arc::clone(&self.dative_bonds),
            Arc::clone(&self.aromatic_systems),
            Arc::clone(&self.multicenter_bonds),
            Arc::clone(&self.noncovalent_bonds),
            self.constraints.clone(),
        )
    }
}

impl Index<AtomIdx> for MoleculeAst {
    type Output = AtomAst;
    fn index(&self, idx: AtomIdx) -> &AtomAst {
        &self.atoms[idx.index()]
    }
}

impl Index<BondIdx> for MoleculeAst {
    type Output = BondAst;
    fn index(&self, idx: BondIdx) -> &BondAst {
        &self.bonds[idx.index()]
    }
}

impl Index<DativeBondIdx> for MoleculeAst {
    type Output = DativeBondAst;
    fn index(&self, idx: DativeBondIdx) -> &DativeBondAst {
        self.dative_bonds.data(RelationId::from(idx))
    }
}

impl Index<AromaticSystemIdx> for MoleculeAst {
    type Output = AromaticSystemAst;
    fn index(&self, idx: AromaticSystemIdx) -> &AromaticSystemAst {
        self.aromatic_systems.data(RelationId::from(idx))
    }
}

impl Index<MulticenterBondIdx> for MoleculeAst {
    type Output = MulticenterBondAst;
    fn index(&self, idx: MulticenterBondIdx) -> &MulticenterBondAst {
        self.multicenter_bonds.data(RelationId::from(idx))
    }
}

impl Index<NoncovalentBondIdx> for MoleculeAst {
    type Output = NoncovalentBondAst;
    fn index(&self, idx: NoncovalentBondIdx) -> &NoncovalentBondAst {
        self.noncovalent_bonds.data(RelationId::from(idx))
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use umol_shared::element::Element;

    use super::super::atom::ElementAst;
    use super::super::constraint::{Constraint, MoleculeConstraint};
    use super::super::dative::DativeBondAst;
    use super::super::multicenter::MulticenterBondAst;
    use super::super::noncovalent::{NoncovalentBondAst, NoncovalentKind};
    use super::super::value::ValueAst;
    use super::*;

    fn ground_atom() -> AtomAst {
        let mut a = AtomAst::from_element(Element::C);
        a.isotope_mass = super::super::atom::IsotopeAst::Natural;
        a.charge = ValueAst::Lit(0);
        a.implicit_hydrogens = super::super::atom::ImplicitHydrogensAst::Value(ValueAst::Lit(4));
        a.lone_pairs = ValueAst::Lit(0);
        a.spin = super::super::spin::SpinStateAst::new(0, 1);
        a
    }

    #[test]
    fn test_molecule_ast_is_ground_empty() {
        assert!(MoleculeAst::default().is_ground());
    }

    #[test]
    fn test_molecule_ast_is_ground_atom() {
        let ast = MoleculeAst::new(
            vec![ground_atom()],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert!(ast.is_ground());
    }

    #[test]
    fn test_molecule_ast_is_ground_wildcard_element() {
        let ast = MoleculeAst::new(
            vec![AtomAst::new(ElementAst::Undetermined)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert!(!ast.is_ground());
    }

    #[test]
    fn test_molecule_ast_is_ground_wildcard_bond() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::new(ValueAst::Undetermined))],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert!(!ast.is_ground());
    }

    #[test]
    fn test_molecule_ast_is_ground_ignores_constraints() {
        let mut ast = MoleculeAst::new(
            vec![ground_atom()],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        ast.constraints
            .push_molecule(Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: vec![],
                sum: ValueAst::Undetermined,
            }));
        assert!(ast.is_ground());
    }

    #[test]
    fn test_molecule_ast_neighbors() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
                AtomAst::from_element(Element::N),
            ],
            vec![
                (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
                (AtomIdx(0), AtomIdx(2), BondAst::from_order(2)),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(ast.neighbors(AtomIdx(0)).count(), 2);
        assert_eq!(ast.neighbors(AtomIdx(1)).count(), 1);
        assert_eq!(ast.neighbors(AtomIdx(2)).count(), 1);
    }

    #[test]
    fn test_molecule_ast_edit_add_aromatic_system() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        let mut b = ast.edit();
        let id = b.add_aromatic_system(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst::default());
        let new_ast = b.build();
        assert_eq!(id, AromaticSystemIdx(0));
        assert_eq!(new_ast.aromatic_systems().count(), 1);
        assert_eq!(ast.aromatic_systems().count(), 0);
    }

    #[test]
    fn test_molecule_ast_counts() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(2))],
            vec![],
            vec![(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst::default())],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(ast.atom_count(), 2);
        assert_eq!(ast.bond_count(), 1);
        assert_eq!(ast.aromatic_system_count(), 1);
        assert_eq!(ast.dative_bond_count(), 0);
        assert_eq!(ast.multicenter_bond_count(), 0);
        assert_eq!(ast.noncovalent_bond_count(), 0);
    }

    fn rich_molecule() -> MoleculeAst {
        // 4 atoms: C(0)—C(1)—N(2)—O(3)
        // 3 covalent bonds: 0–1 (E0), 1–2 (E1), 2–3 (E2)
        // dative: 2→3 (donor=N, acceptor=O)
        // aromatic system: {0,1,2}
        // multicenter bond: {0,1,2}
        // noncovalent: 0↔3
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
                NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond),
            )],
            Constraints::default(),
        )
    }

    #[test]
    fn test_molecule_ast_bond_view() {
        let ast = rich_molecule();
        let bv = ast.bond(BondIdx(0));
        assert_eq!(bv.idx, BondIdx(0));
        assert_eq!(bv.src, AtomIdx(0));
        assert_eq!(bv.tgt, AtomIdx(1));
        assert_eq!(bv.data.order, ValueAst::Lit(1));

        let bv2 = ast.bond(BondIdx(2));
        assert_eq!(bv2.src, AtomIdx(2));
        assert_eq!(bv2.tgt, AtomIdx(3));
    }

    #[test]
    fn test_molecule_ast_bond_views_iter() {
        let ast = rich_molecule();
        let views: Vec<_> = ast.bonds().iter().collect();
        assert_eq!(views.len(), 3);
        assert_eq!(views[0].src, AtomIdx(0));
        assert_eq!(views[1].src, AtomIdx(1));
        assert_eq!(views[2].src, AtomIdx(2));
    }

    #[test]
    fn test_molecule_ast_dative_bond_view() {
        let ast = rich_molecule();
        let dv = ast.dative_bond(DativeBondIdx(0));
        assert_eq!(dv.idx, DativeBondIdx(0));
        assert_eq!(dv.donor, AtomIdx(2));
        assert_eq!(dv.acceptor, AtomIdx(3));
    }

    #[test]
    fn test_molecule_ast_dative_bond_views_iter() {
        let ast = rich_molecule();
        let views: Vec<_> = ast.dative_bonds().iter().collect();
        assert_eq!(views.len(), 1);
        assert_eq!(views[0].donor, AtomIdx(2));
        assert_eq!(views[0].acceptor, AtomIdx(3));
    }

    #[test]
    fn test_molecule_ast_aromatic_system_view() {
        let ast = rich_molecule();
        let av = ast.aromatic_system(AromaticSystemIdx(0));
        assert_eq!(av.idx, AromaticSystemIdx(0));
        let mut atoms: Vec<_> = av.atoms().collect();
        atoms.sort_unstable();
        assert_eq!(atoms, vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
        let mut bonds: Vec<_> = av.bonds().collect();
        bonds.sort_unstable();
        assert_eq!(bonds, vec![BondIdx(0), BondIdx(1)]);
    }

    #[test]
    fn test_molecule_ast_aromatic_system_views_iter() {
        let ast = rich_molecule();
        let views: Vec<_> = ast.aromatic_systems().iter().collect();
        assert_eq!(views.len(), 1);
        assert_eq!(views[0].atoms().count(), 3);
        assert_eq!(views[0].bonds().count(), 2);
    }

    #[test]
    fn test_molecule_ast_multicenter_bond_view() {
        let ast = rich_molecule();
        let mv = ast.multicenter_bond(MulticenterBondIdx(0));
        assert_eq!(mv.idx, MulticenterBondIdx(0));
        let mut atoms: Vec<_> = mv.atoms().collect();
        atoms.sort_unstable();
        assert_eq!(atoms, vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
    }

    #[test]
    fn test_molecule_ast_multicenter_bond_views_iter() {
        let ast = rich_molecule();
        let views: Vec<_> = ast.multicenter_bonds().iter().collect();
        assert_eq!(views.len(), 1);
        assert_eq!(views[0].atoms().count(), 3);
    }

    #[test]
    fn test_molecule_ast_noncovalent_bond_view() {
        let ast = rich_molecule();
        let nv = ast.noncovalent_bond(NoncovalentBondIdx(0));
        assert_eq!(nv.idx, NoncovalentBondIdx(0));
        let mut atoms = nv.atoms;
        atoms.sort_unstable();
        assert_eq!(atoms, [AtomIdx(0), AtomIdx(3)]);
    }

    #[test]
    fn test_molecule_ast_noncovalent_bond_views_iter() {
        let ast = rich_molecule();
        let views: Vec<_> = ast.noncovalent_bonds().iter().collect();
        assert_eq!(views.len(), 1);
    }

    #[test]
    fn test_molecule_ast_connecting_bond() {
        let ast = rich_molecule();
        assert_eq!(ast.connecting_bond(AtomIdx(0), AtomIdx(1)), Some(BondIdx(0)));
        assert_eq!(ast.connecting_bond(AtomIdx(1), AtomIdx(0)), Some(BondIdx(0)));
        assert_eq!(ast.connecting_bond(AtomIdx(0), AtomIdx(3)), None);
    }

    #[test]
    fn test_molecule_ast_dative_bonds_incident() {
        let ast = rich_molecule();
        let inc: Vec<_> = ast.dative_bonds_incident(AtomIdx(2)).collect();
        assert_eq!(inc, vec![DativeBondIdx(0)]);
        let inc: Vec<_> = ast.dative_bonds_incident(AtomIdx(3)).collect();
        assert_eq!(inc, vec![DativeBondIdx(0)]);
        let inc: Vec<_> = ast.dative_bonds_incident(AtomIdx(0)).collect();
        assert!(inc.is_empty());
    }

    #[test]
    fn test_molecule_ast_aromatic_systems_incident() {
        let ast = rich_molecule();
        let inc: Vec<_> = ast.aromatic_systems_incident(AtomIdx(1)).collect();
        assert_eq!(inc, vec![AromaticSystemIdx(0)]);
        let inc: Vec<_> = ast.aromatic_systems_incident(AtomIdx(3)).collect();
        assert!(inc.is_empty());
    }

    #[test]
    fn test_molecule_ast_multicenter_bonds_incident() {
        let ast = rich_molecule();
        let inc: Vec<_> = ast.multicenter_bonds_incident(AtomIdx(0)).collect();
        assert_eq!(inc, vec![MulticenterBondIdx(0)]);
        let inc: Vec<_> = ast.multicenter_bonds_incident(AtomIdx(3)).collect();
        assert!(inc.is_empty());
    }

    #[test]
    fn test_molecule_ast_noncovalent_bonds_incident() {
        let ast = rich_molecule();
        let inc: Vec<_> = ast.noncovalent_bonds_incident(AtomIdx(0)).collect();
        assert_eq!(inc, vec![NoncovalentBondIdx(0)]);
        let inc: Vec<_> = ast.noncovalent_bonds_incident(AtomIdx(3)).collect();
        assert_eq!(inc, vec![NoncovalentBondIdx(0)]);
        let inc: Vec<_> = ast.noncovalent_bonds_incident(AtomIdx(1)).collect();
        assert!(inc.is_empty());
    }

    #[test]
    fn test_molecule_ast_induced_dative_bonds() {
        let ast = rich_molecule();
        assert_eq!(
            ast.induced_dative_bonds(&[AtomIdx(2), AtomIdx(3)]),
            vec![DativeBondIdx(0)]
        );
        assert!(ast.induced_dative_bonds(&[AtomIdx(0), AtomIdx(2)]).is_empty());
    }

    #[test]
    fn test_molecule_ast_induced_aromatic_systems() {
        let ast = rich_molecule();
        assert_eq!(
            ast.induced_aromatic_systems(&[AtomIdx(0), AtomIdx(1), AtomIdx(2)]),
            vec![AromaticSystemIdx(0)]
        );
        assert!(ast.induced_aromatic_systems(&[AtomIdx(0), AtomIdx(1)]).is_empty());
    }

    #[test]
    fn test_molecule_ast_induced_multicenter_bonds() {
        let ast = rich_molecule();
        assert_eq!(
            ast.induced_multicenter_bonds(&[AtomIdx(0), AtomIdx(1), AtomIdx(2)]),
            vec![MulticenterBondIdx(0)]
        );
        assert!(ast.induced_multicenter_bonds(&[AtomIdx(0), AtomIdx(1)]).is_empty());
    }

    #[test]
    fn test_molecule_ast_induced_noncovalent_bonds() {
        let ast = rich_molecule();
        assert_eq!(
            ast.induced_noncovalent_bonds(&[AtomIdx(0), AtomIdx(3)]),
            vec![NoncovalentBondIdx(0)]
        );
        assert!(ast.induced_noncovalent_bonds(&[AtomIdx(0), AtomIdx(1)]).is_empty());
    }

    #[test]
    fn test_molecule_ast_neighbor_view() {
        let ast = rich_molecule();
        let nbrs: Vec<_> = ast.neighbors(AtomIdx(1)).collect();
        assert_eq!(nbrs.len(), 2);
        assert!(nbrs.iter().any(|n| n.atom == AtomIdx(0) && n.bond == BondIdx(0)));
        assert!(nbrs.iter().any(|n| n.atom == AtomIdx(2) && n.bond == BondIdx(1)));
    }

    #[test]
    fn test_molecule_ast_atom_view() {
        let ast = rich_molecule();
        let av = ast.atom(AtomIdx(2));
        assert_eq!(av.idx, AtomIdx(2));
        assert_eq!(av.data.element, ElementAst::Lit(Element::N));
    }

    #[test]
    fn test_molecule_ast_atom_views_iter() {
        let ast = rich_molecule();
        let views: Vec<_> = ast.atoms().iter().collect();
        assert_eq!(views.len(), 4);
        assert_eq!(views[0].idx, AtomIdx(0));
        assert_eq!(views[3].idx, AtomIdx(3));
    }

    #[test]
    fn test_molecule_ast_induced_bonds() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![
                (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
                (AtomIdx(1), AtomIdx(2), BondAst::from_order(1)),
                (AtomIdx(0), AtomIdx(2), BondAst::from_order(1)),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        let bonds = ast.induced_bonds(&[AtomIdx(0), AtomIdx(1)]);
        assert_eq!(bonds, vec![BondIdx(0)]);

        let mut all = ast.induced_bonds(&[AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
        all.sort_unstable();
        assert_eq!(all, vec![BondIdx(0), BondIdx(1), BondIdx(2)]);
    }
}
