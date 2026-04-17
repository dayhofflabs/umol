//! Molecule structural AST.

use std::ops::Index;
use std::sync::Arc;

use umol_graph_core::relation::RelationId;
use umol_graph_core::{
    EdgeId, FixedRelationSet, Graph, NodeId, VarRelationSet,
};
use umol_shared::value_ast::ValueAst;

use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::builder::MoleculeBuilder;
use crate::ast::config::MoleculeAstConfig;
use crate::ast::constraint::MoleculeConstraint;
use crate::ast::views::{
    AromaticSystemViews, AtomView, AtomViewMut, AtomViews, BondView, BondViewMut, BondViews,
    DativeBondViews, MulticenterBondViews, NeighborView, NoncovalentBondViews,
};
use crate::ast::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx,
    NoncovalentBondIdx,
};
use crate::table_ir::Molecule as TableMolecule;
use crate::table_ir::bond::BondDonation;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AromaticSystemAst {}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MulticenterBondAst {}

/// Molecule AST: structural representation of a molecule (ground or pattern).
///
/// Topology and per-atom/bond data are `Arc`-shared (copy-on-write) so
/// derived molecules in a reaction network share storage. The AST itself
/// only allows attribute mutation (`atom_mut`, `bond_mut`); structural
/// edits go through `MoleculeBuilder` via [`MoleculeAst::edit`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeAst {
    graph: Graph,
    atoms: Arc<Vec<AtomAst>>,
    bonds: Arc<Vec<BondAst>>,
    dative_bonds: Arc<FixedRelationSet<BondAst, 2>>,
    noncovalent_bonds: Arc<FixedRelationSet<BondAst, 2>>,
    aromatic_systems: Arc<VarRelationSet<AromaticSystemAst>>,
    multicenter_bonds: Arc<VarRelationSet<MulticenterBondAst>>,
    pub constraints: Vec<MoleculeConstraint>,
}

impl Default for MoleculeAst {
    fn default() -> Self {
        Self {
            graph: Graph::default(),
            atoms: Arc::new(Vec::new()),
            bonds: Arc::new(Vec::new()),
            dative_bonds: Arc::new(FixedRelationSet::default()),
            noncovalent_bonds: Arc::new(FixedRelationSet::default()),
            aromatic_systems: Arc::new(VarRelationSet::default()),
            multicenter_bonds: Arc::new(VarRelationSet::default()),
            constraints: Vec::new(),
        }
    }
}

impl MoleculeAst {
    pub fn new(
        atoms: Vec<AtomAst>,
        bonds: Vec<(AtomIdx, AtomIdx, BondAst)>,
        dative: Vec<(AtomIdx, AtomIdx, BondAst)>,
        noncovalent: Vec<(AtomIdx, AtomIdx, BondAst)>,
        aromatic: Vec<(Vec<AtomIdx>, AromaticSystemAst)>,
        multicenter: Vec<(Vec<AtomIdx>, MulticenterBondAst)>,
        constraints: Vec<MoleculeConstraint>,
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

        let noncovalent_bonds = FixedRelationSet::new(
            noncovalent
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

        Self {
            graph,
            atoms: Arc::new(atoms),
            bonds: Arc::new(bond_data),
            dative_bonds: Arc::new(dative_bonds),
            noncovalent_bonds: Arc::new(noncovalent_bonds),
            aromatic_systems: Arc::new(aromatic_systems),
            multicenter_bonds: Arc::new(multicenter_bonds),
            constraints,
        }
    }

    pub(crate) fn from_arcs(
        graph: Graph,
        atoms: Arc<Vec<AtomAst>>,
        bonds: Arc<Vec<BondAst>>,
        dative_bonds: Arc<FixedRelationSet<BondAst, 2>>,
        noncovalent_bonds: Arc<FixedRelationSet<BondAst, 2>>,
        aromatic_systems: Arc<VarRelationSet<AromaticSystemAst>>,
        multicenter_bonds: Arc<VarRelationSet<MulticenterBondAst>>,
        constraints: Vec<MoleculeConstraint>,
    ) -> Self {
        Self {
            graph,
            atoms,
            bonds,
            dative_bonds,
            noncovalent_bonds,
            aromatic_systems,
            multicenter_bonds,
            constraints,
        }
    }

    /// Lift a `table_ir::Molecule` to a `MoleculeAst` by lifting atoms and
    /// bonds individually (`AtomAst::from_table_atom`, `BondAst::from_table_bond`)
    /// and splitting bonds into regular/dative/noncovalent by their table-level
    /// tags. Aromatic systems and constraints are not derived here.
    pub fn from_table_molecule(mol: &TableMolecule) -> Self {
        let atoms: Vec<AtomAst> = mol.atoms.iter().map(AtomAst::from_table_atom).collect();

        let mut regular = Vec::new();
        let mut dative = Vec::new();
        let mut noncovalent = Vec::new();
        for b in &mol.bonds {
            let a_idx = AtomIdx(b.atoms.first());
            let b_idx = AtomIdx(b.atoms.second());
            let bond_ast = BondAst::from_table_bond(b);
            if b.noncovalent.is_some() {
                noncovalent.push((a_idx, b_idx, bond_ast));
            } else if matches!(
                b.donation,
                Some(BondDonation::Donating | BondDonation::Accepting)
            ) {
                dative.push((a_idx, b_idx, bond_ast));
            } else {
                regular.push((a_idx, b_idx, bond_ast));
            }
        }

        let multicenter: Vec<(Vec<AtomIdx>, MulticenterBondAst)> = mol
            .multicenter_bonds
            .iter()
            .map(|mc| {
                let mut seen = std::collections::HashSet::new();
                let atoms: Vec<AtomIdx> = mc
                    .all_atoms()
                    .into_iter()
                    .filter(|a| seen.insert(*a))
                    .map(AtomIdx)
                    .collect();
                (atoms, MulticenterBondAst {})
            })
            .collect();

        Self::new(atoms, regular, dative, noncovalent, vec![], multicenter, vec![])
    }
}

impl MoleculeAst {
    pub fn atoms(&self) -> AtomViews<'_> {
        AtomViews { atoms: &self.atoms }
    }

    pub fn bonds(&self) -> BondViews<'_> {
        BondViews { bonds: &self.bonds, graph: &self.graph }
    }

    pub fn dative_bonds(&self) -> DativeBondViews<'_> {
        DativeBondViews { set: &self.dative_bonds }
    }

    pub fn noncovalent_bonds(&self) -> NoncovalentBondViews<'_> {
        NoncovalentBondViews { set: &self.noncovalent_bonds }
    }

    pub fn aromatic_systems(&self) -> AromaticSystemViews<'_> {
        AromaticSystemViews { set: &self.aromatic_systems }
    }

    pub fn multicenter_bonds(&self) -> MulticenterBondViews<'_> {
        MulticenterBondViews { set: &self.multicenter_bonds }
    }

    pub fn atom(&self, idx: AtomIdx) -> AtomView<'_> {
        self.atoms().get(idx)
    }

    pub fn atom_mut(&mut self, idx: AtomIdx) -> AtomViewMut<'_> {
        let data = &mut Arc::make_mut(&mut self.atoms)[idx.index()];
        AtomViewMut { idx, data }
    }

    pub fn bond(&self, idx: BondIdx) -> BondView<'_> {
        self.bonds().get(idx)
    }

    pub fn bond_mut(&mut self, idx: BondIdx) -> BondViewMut<'_> {
        let [s, t] = self.graph.edge_endpoints(EdgeId::from(idx));
        let data = &mut Arc::make_mut(&mut self.bonds)[idx.index()];
        BondViewMut {
            idx,
            src: AtomIdx::from(s),
            tgt: AtomIdx::from(t),
            data,
        }
    }

    pub fn neighbors(&self, atom: AtomIdx) -> impl Iterator<Item = NeighborView<'_>> {
        let bonds = &self.bonds;
        self.graph.neighbors(NodeId::from(atom)).iter().map(move |n| NeighborView {
            atom: AtomIdx::from(n.node),
            bond: BondIdx::from(n.edge),
            data: &bonds[n.edge.index()],
        })
    }

    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    pub fn edit(&self) -> MoleculeBuilder {
        MoleculeBuilder::from_parts(
            self.graph.clone(),
            Arc::clone(&self.atoms),
            Arc::clone(&self.bonds),
            Arc::clone(&self.dative_bonds),
            Arc::clone(&self.noncovalent_bonds),
            Arc::clone(&self.aromatic_systems),
            Arc::clone(&self.multicenter_bonds),
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
    type Output = BondAst;
    fn index(&self, idx: DativeBondIdx) -> &BondAst {
        self.dative_bonds.data(RelationId::from(idx))
    }
}

impl Index<NoncovalentBondIdx> for MoleculeAst {
    type Output = BondAst;
    fn index(&self, idx: NoncovalentBondIdx) -> &BondAst {
        self.noncovalent_bonds.data(RelationId::from(idx))
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

impl MoleculeAst {
    pub fn coerce(&mut self, config: &MoleculeAstConfig) {
        for atom in Arc::make_mut(&mut self.atoms) {
            atom.coerce(&config.atom);
        }
    }

    pub fn release(&mut self, config: &MoleculeAstConfig) {
        for atom in Arc::make_mut(&mut self.atoms) {
            atom.release(&config.atom);
        }
    }

    pub fn is_ground(&self) -> bool {
        self.atoms.iter().all(|a| a.is_ground())
            && self.bonds.iter().all(|b| b.is_ground())
            && self
                .dative_bonds
                .relation_ids()
                .all(|id| self.dative_bonds.data(id).is_ground())
            && self
                .noncovalent_bonds
                .relation_ids()
                .all(|id| self.noncovalent_bonds.data(id).is_ground())
            && self.constraints.iter().all(|c| c.is_ground_assertion())
    }

    pub fn bond_order_sum(&self, atom: AtomIdx) -> Option<u8> {
        let mut sum: u8 = 0;
        for n in self.graph.neighbors(NodeId::from(atom)) {
            match self.bonds[n.edge.index()].order {
                ValueAst::Lit(v) => sum += v as u8,
                _ => return None,
            }
        }
        Some(sum)
    }

    pub fn dative_bond_order_sums(&self, atom: AtomIdx) -> (u8, u8) {
        let node = NodeId::from(atom);
        let mut donated: u8 = 0;
        let mut accepted: u8 = 0;
        for &rel_id in self.dative_bonds.incident(node) {
            let data = self.dative_bonds.data(rel_id);
            let order = match data.order {
                ValueAst::Lit(v) => v as u8,
                _ => continue,
            };
            let participants = self.dative_bonds.participants(rel_id);
            if participants[0] == node {
                donated += order;
            } else {
                accepted += order;
            }
        }
        (donated, accepted)
    }

    pub fn is_in_aromatic_system(&self, atom: AtomIdx) -> bool {
        self.aromatic_systems.has_incident(NodeId::from(atom))
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use umol_shared::atom_ast::{AromaticValenceAst, ElementAst, HydrogenAst, IsotopeAst};
    use umol_shared::element::Element;
    use umol_shared::spin::SpinState;
    use umol_shared::spin_ast::SpinStateAst;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use crate::ast::constraint::{DerivedPred, RelationRefs};

    fn ground_atom() -> AtomAst {
        AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: HydrogenAst::Value(ValueAst::Lit(4)),
            lone_pairs: ValueAst::Lit(0),
            spin: SpinStateAst::Lit(SpinState::closed_shell()),
            valence: ValueAst::Lit(4),
            donated_pairs: ValueAst::Lit(0),
            accepted_pairs: ValueAst::Lit(0),
            aromatic_valence: AromaticValenceAst::NotAromatic,
            multicenter_valence: ValueAst::Lit(0),
        }
    }

    fn ground_ast() -> MoleculeAst {
        MoleculeAst::new(vec![ground_atom()], vec![], vec![], vec![], vec![], vec![], vec![])
    }

    #[test]
    fn test_molecule_ast_is_ground_empty() {
        assert!(MoleculeAst::default().is_ground());
    }

    #[test]
    fn test_molecule_ast_is_ground_atom() {
        assert!(ground_ast().is_ground());
    }

    #[test]
    fn test_molecule_ast_is_ground_with_constraint() {
        let mut ast = ground_ast();
        ast.constraints.push(MoleculeConstraint::Derived {
            predicate: DerivedPred::TotalCharge(ValueAst::Lit(-1)),
            refs: RelationRefs::default(),
        });
        assert!(ast.is_ground());
    }

    #[test]
    fn test_molecule_ast_is_ground_wildcard_element() {
        let ast = MoleculeAst::new(
            vec![AtomAst::new(ElementAst::Undetermined)],
            vec![], vec![], vec![], vec![], vec![], vec![],
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
            vec![], vec![], vec![], vec![], vec![],
        );
        assert!(!ast.is_ground());
    }

    #[test]
    fn test_molecule_ast_is_ground_non_ground_constraint() {
        let mut ast = ground_ast();
        ast.constraints.push(MoleculeConstraint::Derived {
            predicate: DerivedPred::TotalSpin(SpinStateAst::default()),
            refs: RelationRefs::default(),
        });
        assert!(!ast.is_ground());
    }

    #[test]
    fn test_molecule_ast_is_ground_sub_pattern() {
        let mut ast = ground_ast();
        ast.constraints.push(MoleculeConstraint::SubPattern {
            anchor: AtomIdx(0),
            pattern: Box::new(MoleculeAst::default()),
        });
        assert!(!ast.is_ground());
    }

    #[test]
    fn test_molecule_ast_bond_order_sum() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(2))],
            vec![], vec![], vec![], vec![], vec![],
        );
        assert_eq!(ast.bond_order_sum(AtomIdx(0)), Some(2));
        assert_eq!(ast.bond_order_sum(AtomIdx(1)), Some(2));
    }

    #[test]
    fn test_molecule_ast_bond_order_sum_no_bonds() {
        let ast = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C)],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        assert_eq!(ast.bond_order_sum(AtomIdx(0)), Some(0));
    }

    #[test]
    fn test_molecule_ast_bond_order_sum_wildcard() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::new(ValueAst::Undetermined))],
            vec![], vec![], vec![], vec![], vec![],
        );
        assert_eq!(ast.bond_order_sum(AtomIdx(0)), None);
    }

    #[test]
    fn test_molecule_ast_is_in_aromatic_system() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![], vec![],
            vec![(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst {})],
            vec![], vec![],
        );
        assert!(ast.is_in_aromatic_system(AtomIdx(0)));
        assert!(ast.is_in_aromatic_system(AtomIdx(1)));
        assert!(!ast.is_in_aromatic_system(AtomIdx(2)));
    }

    #[test]
    fn test_molecule_ast_dative_bond_order_sums() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::B),
            ],
            vec![],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![], vec![], vec![], vec![],
        );
        assert_eq!(ast.dative_bond_order_sums(AtomIdx(0)), (1, 0));
        assert_eq!(ast.dative_bond_order_sums(AtomIdx(1)), (0, 1));
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
            vec![], vec![], vec![], vec![], vec![],
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
            vec![], vec![], vec![], vec![], vec![],
        );
        let mut b = ast.edit();
        let id = b.add_aromatic_system(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst {});
        let new_ast = b.build();
        assert_eq!(id, AromaticSystemIdx(0));
        assert_eq!(new_ast.aromatic_systems().count(), 1);
        assert!(new_ast.is_in_aromatic_system(AtomIdx(0)));
        assert_eq!(ast.aromatic_systems().count(), 0);
    }
}
