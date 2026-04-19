//! Molecule structural AST.

use std::collections::HashSet;
use std::ops::Index;
use std::sync::Arc;

use umol_graph_core::relation::RelationId;
use umol_graph_core::{EdgeId, FixedRelationSet, Graph, NodeId, VarRelationSet};
use umol_shared::value_ast::ValueAst;

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::builder::MoleculeBuilder;
use super::config::MoleculeAstConfig;
use super::constraint::{
    AromaticValenceConstraint, AtomConstraint, AtomConstraintKind, MoleculeConstraint,
    MoleculeConstraints,
};
use super::multicenter::MulticenterBondAst;
use super::views::{
    AromaticSystemViews, AtomView, AtomViewMut, AtomViews, BondView, BondViewMut, BondViews,
    DativeBondViews, MulticenterBondViews, NeighborView, NoncovalentBondViews,
};
use super::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use crate::api::pattern::{coerce_atom_constraints, release_atom_constraints};
use crate::table_ir::bond::BondDonation;
use crate::table_ir::Molecule as TableMolecule;

/// Molecule AST: structural representation of a molecule (ground or pattern).
///
/// Topology and per-atom/bond data are `Arc`-shared (copy-on-write).
/// AST itself only allows attribute mutation (`atom_mut`, `bond_mut`);
/// structural edits go through `MoleculeBuilder` via [`MoleculeAst::edit`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeAst {
    graph: Graph,
    atoms: Arc<Vec<AtomAst>>,
    bonds: Arc<Vec<BondAst>>,
    dative_bonds: Arc<FixedRelationSet<BondAst, 2>>,
    noncovalent_bonds: Arc<FixedRelationSet<BondAst, 2>>,
    aromatic_systems: Arc<VarRelationSet<AromaticSystemAst>>,
    multicenter_bonds: Arc<VarRelationSet<MulticenterBondAst>>,
    constraints: MoleculeConstraints,
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
            constraints: MoleculeConstraints::new(),
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
        constraints: impl Into<MoleculeConstraints>,
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
            constraints: constraints.into(),
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn from_arcs(
        graph: Graph,
        atoms: Arc<Vec<AtomAst>>,
        bonds: Arc<Vec<BondAst>>,
        dative_bonds: Arc<FixedRelationSet<BondAst, 2>>,
        noncovalent_bonds: Arc<FixedRelationSet<BondAst, 2>>,
        aromatic_systems: Arc<VarRelationSet<AromaticSystemAst>>,
        multicenter_bonds: Arc<VarRelationSet<MulticenterBondAst>>,
        constraints: MoleculeConstraints,
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
    /// tags. Per-atom aromatic hints and CTAB valence overrides lift to the
    /// molecule constraint vec.
    pub fn from_table_molecule(mol: &TableMolecule) -> Self {
        let atoms: Vec<AtomAst> = mol.atoms.iter().map(AtomAst::from_table_atom).collect();

        let mut constraints: MoleculeConstraints = MoleculeConstraints::new();
        for (i, atom) in mol.atoms.iter().enumerate() {
            let idx = AtomIdx::from(i);
            if atom.aromatic == Some(true) {
                constraints.insert(MoleculeConstraint::AtomPred(
                    idx,
                    AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(
                        ValueAst::Undetermined,
                    )),
                ));
            }
            if let Some(v) = atom.valence {
                constraints.insert(MoleculeConstraint::AtomPred(
                    idx,
                    AtomConstraint::Valence(ValueAst::Lit(v as i64)),
                ));
            }
        }

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
                let mut seen = HashSet::new();
                let atoms: Vec<AtomIdx> = mc
                    .all_atoms()
                    .into_iter()
                    .filter(|a| seen.insert(*a))
                    .map(AtomIdx)
                    .collect();
                (atoms, MulticenterBondAst {})
            })
            .collect();

        Self::new(
            atoms,
            regular,
            dative,
            noncovalent,
            vec![],
            multicenter,
            constraints,
        )
    }
}

impl MoleculeAst {
    pub fn atoms(&self) -> AtomViews<'_> {
        AtomViews { atoms: &self.atoms }
    }

    pub fn bonds(&self) -> BondViews<'_> {
        BondViews {
            bonds: &self.bonds,
            graph: &self.graph,
        }
    }

    pub fn dative_bonds(&self) -> DativeBondViews<'_> {
        DativeBondViews {
            set: &self.dative_bonds,
        }
    }

    pub fn noncovalent_bonds(&self) -> NoncovalentBondViews<'_> {
        NoncovalentBondViews {
            set: &self.noncovalent_bonds,
        }
    }

    pub fn aromatic_systems(&self) -> AromaticSystemViews<'_> {
        AromaticSystemViews {
            set: &self.aromatic_systems,
        }
    }

    pub fn multicenter_bonds(&self) -> MulticenterBondViews<'_> {
        MulticenterBondViews {
            set: &self.multicenter_bonds,
        }
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
        let [s, t] = self.graph.edge_endpoints(EdgeId::from(idx));
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

    pub fn aromatic_systems_mut(&mut self) -> impl Iterator<Item = &mut AromaticSystemAst> {
        Arc::make_mut(&mut self.aromatic_systems).data_iter_mut()
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

    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    pub fn constraints(&self) -> &MoleculeConstraints {
        &self.constraints
    }

    pub fn constraints_mut(&mut self) -> &mut MoleculeConstraints {
        &mut self.constraints
    }

    /// Per-atom aromatic-valence constraint hint, if any. Returns `None` when no
    /// `AtomPred(_, AromaticValence(_))` entry is attached to this atom.
    pub fn atom_aromatic_hint(&self, idx: AtomIdx) -> Option<&AromaticValenceConstraint> {
        match self
            .constraints
            .atoms()
            .get(&idx)?
            .get(AtomConstraintKind::AromaticValence)?
        {
            AtomConstraint::AromaticValence(c) => Some(c),
            _ => None,
        }
    }

    /// True if the atom is hinted aromatic via constraint or already placed in an aromatic system.
    pub fn atom_is_aromatic(&self, idx: AtomIdx) -> bool {
        if self.is_in_aromatic_system(idx) {
            return true;
        }
        matches!(
            self.atom_aromatic_hint(idx),
            Some(AromaticValenceConstraint::Value(_))
        )
    }

    /// Aromatic pi-electron count for this atom, if pinned to a literal in the
    /// constraint vec. `None` when absent, set to `Undetermined`, or `NotAromatic`.
    pub fn atom_aromatic_valence(&self, idx: AtomIdx) -> Option<u8> {
        match self.atom_aromatic_hint(idx)? {
            AromaticValenceConstraint::Value(ValueAst::Lit(n)) => Some(*n as u8),
            _ => None,
        }
    }

    /// Multicenter-bond electron contribution for this atom, if pinned to a
    /// literal in the constraint vec. `None` when absent or non-literal.
    pub fn atom_multicenter_valence(&self, idx: AtomIdx) -> Option<u8> {
        match self
            .constraints
            .atoms()
            .get(&idx)?
            .get(AtomConstraintKind::MulticenterValence)?
        {
            AtomConstraint::MulticenterValence(ValueAst::Lit(n)) => Some(*n as u8),
            _ => None,
        }
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
        for bond in Arc::make_mut(&mut self.bonds) {
            bond.coerce(&config.bond);
        }
        let dative = Arc::make_mut(&mut self.dative_bonds);
        for rid in dative.relation_ids().collect::<Vec<_>>() {
            dative.data_mut(rid).coerce(&config.bond);
        }
        let noncov = Arc::make_mut(&mut self.noncovalent_bonds);
        for rid in noncov.relation_ids().collect::<Vec<_>>() {
            noncov.data_mut(rid).coerce(&config.bond);
        }
        let aromatic = Arc::make_mut(&mut self.aromatic_systems);
        for rid in aromatic.relation_ids().collect::<Vec<_>>() {
            aromatic.data_mut(rid).coerce(&config.bond);
        }
        for i in 0..self.atoms.len() {
            let idx = AtomIdx::from(i);
            let set = self.constraints.atoms_mut().entry(idx).or_default();
            let mut bucket: Vec<AtomConstraint> = set.iter().cloned().collect();
            coerce_atom_constraints(&mut bucket, &config.atom);
            *set = bucket.into_iter().collect();
        }
    }

    pub fn release(&mut self, config: &MoleculeAstConfig) {
        for atom in Arc::make_mut(&mut self.atoms) {
            atom.release(&config.atom);
        }
        for bond in Arc::make_mut(&mut self.bonds) {
            bond.release(&config.bond);
        }
        let dative = Arc::make_mut(&mut self.dative_bonds);
        for rid in dative.relation_ids().collect::<Vec<_>>() {
            dative.data_mut(rid).release(&config.bond);
        }
        let noncov = Arc::make_mut(&mut self.noncovalent_bonds);
        for rid in noncov.relation_ids().collect::<Vec<_>>() {
            noncov.data_mut(rid).release(&config.bond);
        }
        let aromatic = Arc::make_mut(&mut self.aromatic_systems);
        for rid in aromatic.relation_ids().collect::<Vec<_>>() {
            aromatic.data_mut(rid).release(&config.bond);
        }
        for set in self.constraints.atoms_mut().values_mut() {
            set.retain(|ac| {
                let mut tmp = vec![ac.clone()];
                release_atom_constraints(&mut tmp, &config.atom);
                !tmp.is_empty()
            });
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
            && self
                .aromatic_systems
                .relation_ids()
                .all(|id| self.aromatic_systems.data(id).is_ground())
            && self
                .multicenter_bonds
                .relation_ids()
                .all(|id| self.multicenter_bonds.data(id).is_ground())
            && self.constraints.iter().all(|c| c.is_ground_assertion())
    }

    pub fn valence(&self, atom: AtomIdx) -> Option<u8> {
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

    /// Bonds of the induced subgraph over `atoms`: edges whose both endpoints
    /// are in the set. Result is sorted by `BondIdx`.
    pub fn induced_bonds(&self, atoms: &[AtomIdx]) -> Vec<BondIdx> {
        let set: HashSet<AtomIdx> = atoms.iter().copied().collect();
        let mut seen: HashSet<BondIdx> = HashSet::new();
        let mut bonds: Vec<BondIdx> = Vec::new();
        for &a in atoms {
            for n in self.graph.neighbors(NodeId::from(a)) {
                let other = AtomIdx::from(n.node);
                if set.contains(&other) {
                    let bond = BondIdx::from(n.edge);
                    if seen.insert(bond) {
                        bonds.push(bond);
                    }
                }
            }
        }
        bonds.sort_unstable();
        bonds
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use umol_shared::atom_ast::{ElementAst, HydrogenAst, IsotopeAst};
    use umol_shared::element::Element;
    use umol_shared::spin::SpinState;
    use umol_shared::spin_ast::SpinStateAst;
    use umol_shared::value_ast::ValueAst;

    use super::*;

    fn ground_atom() -> AtomAst {
        AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: HydrogenAst::Value(ValueAst::Lit(4)),
            lone_pairs: ValueAst::Lit(0),
            spin: SpinStateAst::Lit(SpinState::closed_shell()),
        }
    }

    fn ground_ast() -> MoleculeAst {
        MoleculeAst::new(
            vec![ground_atom()],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        )
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
        ast.constraints
            .insert(MoleculeConstraint::TotalCharge(ValueAst::Lit(-1)));
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
            vec![],
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
            vec![],
        );
        assert!(!ast.is_ground());
    }

    #[test]
    fn test_molecule_ast_is_ground_non_ground_constraint() {
        let mut ast = ground_ast();
        ast.constraints
            .insert(MoleculeConstraint::TotalSpin(SpinStateAst::default()));
        assert!(!ast.is_ground());
    }

    #[test]
    fn test_molecule_ast_is_ground_sub_pattern() {
        let mut ast = ground_ast();
        ast.constraints.insert(MoleculeConstraint::SubPattern {
            target_anchor: AtomIdx(0),
            pattern_anchor: AtomIdx(0),
            pattern: Box::new(MoleculeAst::default()),
        });
        assert!(!ast.is_ground());
    }

    #[test]
    fn test_molecule_ast_valence() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(2))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(ast.valence(AtomIdx(0)), Some(2));
        assert_eq!(ast.valence(AtomIdx(1)), Some(2));
    }

    #[test]
    fn test_molecule_ast_valence_no_bonds() {
        let ast = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(ast.valence(AtomIdx(0)), Some(0));
    }

    #[test]
    fn test_molecule_ast_valence_wildcard() {
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
            vec![],
        );
        assert_eq!(ast.valence(AtomIdx(0)), None);
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
            vec![],
            vec![],
            vec![(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst::default())],
            vec![],
            vec![],
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
            vec![],
            vec![],
            vec![],
            vec![],
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
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
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
            vec![],
        );
        let mut b = ast.edit();
        let id = b.add_aromatic_system(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst::default());
        let new_ast = b.build();
        assert_eq!(id, AromaticSystemIdx(0));
        assert_eq!(new_ast.aromatic_systems().count(), 1);
        assert!(new_ast.is_in_aromatic_system(AtomIdx(0)));
        assert_eq!(ast.aromatic_systems().count(), 0);
    }
}
