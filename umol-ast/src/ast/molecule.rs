//! Molecule structural AST.

use std::ops::Index;
use std::sync::Arc;

use umol_graph_core::relation::RelationId;
use umol_graph_core::{FixedRelationSet, Graph, NodeId, VarRelationSet};

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::builder::MoleculeBuilder;
use super::constraint::Constraints;
use super::dative::DativeBondAst;
use super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;
use super::views::{
    AromaticSystemViews, AtomView, AtomViewMut, AtomViews, BondView, BondViewMut, BondViews,
    DativeBondViews, MulticenterBondViews, NeighborView, NoncovalentBondViews,
};

/// Molecule AST: structural representation of a molecule (ground or pattern).
///
/// Topology and per-atom/bond data are `Arc`-shared (copy-on-write). The AST
/// itself only allows attribute mutation (`atom_mut`, `bond_mut`); structural
/// edits go through `MoleculeBuilder` via [`MoleculeAst::edit`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeAst {
    pub(super) graph: Graph,
    pub(super) atoms: Arc<Vec<AtomAst>>,
    pub(super) bonds: Arc<Vec<BondAst>>,
    pub(super) dative_bonds: Arc<FixedRelationSet<DativeBondAst, 2>>,
    pub(super) aromatic_systems: Arc<VarRelationSet<AromaticSystemAst>>,
    pub(super) multicenter_bonds: Arc<VarRelationSet<MulticenterBondAst>>,
    pub(super) noncovalent_bonds: Arc<FixedRelationSet<NoncovalentBondAst, 2>>,
    pub(super) constraints: Constraints,
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
    pub(super) fn from_arcs(
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

    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    pub fn atoms(&self) -> AtomViews<'_> {
        AtomViews::new(&self.atoms)
    }

    pub fn bonds(&self) -> BondViews<'_> {
        BondViews::new(&self.bonds, &self.graph)
    }

    pub fn dative_bonds(&self) -> DativeBondViews<'_> {
        DativeBondViews::new(&self.dative_bonds)
    }

    pub fn aromatic_systems(&self) -> AromaticSystemViews<'_> {
        AromaticSystemViews::new(&self.aromatic_systems)
    }

    pub fn multicenter_bonds(&self) -> MulticenterBondViews<'_> {
        MulticenterBondViews::new(&self.multicenter_bonds)
    }

    pub fn noncovalent_bonds(&self) -> NoncovalentBondViews<'_> {
        NoncovalentBondViews::new(&self.noncovalent_bonds)
    }

    pub fn atom(&self, idx: AtomIdx) -> AtomView<'_> {
        self.atoms().get(idx)
    }

    pub fn atom_mut(&mut self, idx: AtomIdx) -> AtomViewMut<'_> {
        let data = &mut Arc::make_mut(&mut self.atoms)[idx.index()];
        AtomViewMut { idx, data }
    }

    pub fn atoms_mut(&mut self) -> impl Iterator<Item = &mut AtomAst> {
        Arc::make_mut(&mut self.atoms).iter_mut()
    }

    pub fn bond(&self, idx: BondIdx) -> BondView<'_> {
        self.bonds().get(idx)
    }

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

    pub fn constraints(&self) -> &Constraints {
        &self.constraints
    }

    pub fn constraints_mut(&mut self) -> &mut Constraints {
        &mut self.constraints
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
            vec![(
                AtomIdx(0),
                AtomIdx(1),
                BondAst::new(ValueAst::Undetermined),
            )],
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
        let id =
            b.add_aromatic_system(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst::default());
        let new_ast = b.build();
        assert_eq!(id, AromaticSystemIdx(0));
        assert_eq!(new_ast.aromatic_systems().count(), 1);
        assert_eq!(ast.aromatic_systems().count(), 0);
    }
}
