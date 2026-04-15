//! Molecule structural AST.

use index_vec::Idx;
use umol_graph_core::relation::RelationId;
use umol_graph_core::{
    EdgeId, FixedRelationSet, Graph, Neighbor, NodeId, Remapping, VarRelationSet,
};
use umol_shared::value_ast::ValueAst;

use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::constraint::MoleculeConstraint;
use crate::ast::error::GroundError;
use crate::ast::{
    AromaticSystemIdx, Ast, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx,
    NoncovalentBondIdx,
};
use crate::ast::config::MoleculeAstConfig;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AromaticSystemAst {}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MulticenterBondAst {}

/// Molecule AST: structural representation of a molecule (ground or pattern).
///
/// Topology lives in a CSR `Graph` (Arc-shared, copy-on-write).
/// Atom and bond data are stored in flat arrays indexed by graph position.
/// Secondary relations (dative, noncovalent, aromatic, multicenter) use
/// flat relation sets sharing the same `NodeId` space.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MoleculeAst {
    graph: Graph,
    atoms: Vec<AtomAst>,
    bonds: Vec<BondAst>,
    dative_bonds: FixedRelationSet<BondAst, 2>,
    noncovalent_bonds: FixedRelationSet<BondAst, 2>,
    aromatic_systems: VarRelationSet<AromaticSystemAst>,
    multicenter_bonds: VarRelationSet<MulticenterBondAst>,
    pub constraints: Vec<MoleculeConstraint>,
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

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
            atoms,
            bonds: bond_data,
            dative_bonds,
            noncovalent_bonds,
            aromatic_systems,
            multicenter_bonds,
            constraints,
        }
    }
}

// ---------------------------------------------------------------------------
// Mutations
// ---------------------------------------------------------------------------

impl MoleculeAst {
    pub fn add_atom(&mut self, atom: AtomAst) -> AtomIdx {
        let node_id = self.graph.add_node();
        self.atoms.push(atom);
        AtomIdx::from(node_id)
    }

    pub fn add_bond(&mut self, source: AtomIdx, target: AtomIdx, bond: BondAst) -> BondIdx {
        let edge_id = self.graph.add_edge(NodeId::from(source), NodeId::from(target));
        self.bonds.push(bond);
        BondIdx::from(edge_id)
    }

    pub fn remove_atom(&mut self, idx: AtomIdx) -> Remapping {
        let remap = self.graph.remove_node(NodeId::from(idx));
        self.apply_remapping(&remap);
        remap
    }

    pub fn remove_bond(&mut self, idx: BondIdx) -> Remapping {
        let remap = self.graph.remove_edge(EdgeId::from(idx));
        self.apply_remapping(&remap);
        remap
    }

    pub fn remove(&mut self, atoms: &[AtomIdx], bonds: &[BondIdx]) -> Remapping {
        let nodes: Vec<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        let edges: Vec<EdgeId> = bonds.iter().map(|&b| EdgeId::from(b)).collect();
        let remap = self.graph.remove(&nodes, &edges);
        self.apply_remapping(&remap);
        remap
    }

    pub fn set_aromatic_systems(&mut self, systems: Vec<(Vec<AtomIdx>, AromaticSystemAst)>) {
        self.aromatic_systems = VarRelationSet::new(
            systems
                .into_iter()
                .map(|(atoms, d)| (atoms.into_iter().map(NodeId::from).collect(), d))
                .collect(),
        );
    }

    fn apply_remapping(&mut self, remap: &Remapping) {
        self.atoms = remap.apply_to_node_vec(&self.atoms);
        self.bonds = remap.apply_to_edge_vec(&self.bonds);
        self.dative_bonds = remap.apply_to_fixed_relation_set(&self.dative_bonds);
        self.noncovalent_bonds = remap.apply_to_fixed_relation_set(&self.noncovalent_bonds);
        self.aromatic_systems = remap.apply_to_var_relation_set(&self.aromatic_systems);
        self.multicenter_bonds = remap.apply_to_var_relation_set(&self.multicenter_bonds);
    }
}

// ---------------------------------------------------------------------------
// Atom accessors
// ---------------------------------------------------------------------------

impl MoleculeAst {
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn atom(&self, idx: AtomIdx) -> &AtomAst {
        &self.atoms[idx.index()]
    }

    pub fn atom_mut(&mut self, idx: AtomIdx) -> &mut AtomAst {
        &mut self.atoms[idx.index()]
    }

    pub fn atoms(&self) -> impl Iterator<Item = (AtomIdx, &AtomAst)> {
        self.atoms
            .iter()
            .enumerate()
            .map(|(i, a)| (AtomIdx(i as u32), a))
    }

    pub fn neighbors(&self, idx: AtomIdx) -> &[Neighbor] {
        self.graph.neighbors(NodeId::from(idx))
    }
}

// ---------------------------------------------------------------------------
// Bond accessors
// ---------------------------------------------------------------------------

impl MoleculeAst {
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    pub fn bond(&self, idx: BondIdx) -> &BondAst {
        &self.bonds[idx.index()]
    }

    pub fn bond_endpoints(&self, idx: BondIdx) -> (AtomIdx, AtomIdx) {
        let [a, b] = self.graph.edge_endpoints(EdgeId::from(idx));
        (AtomIdx::from(a), AtomIdx::from(b))
    }

    pub fn bonds(&self) -> impl Iterator<Item = (BondIdx, AtomIdx, AtomIdx, &BondAst)> {
        self.graph.edge_ids().map(|id| {
            let [a, b] = self.graph.edge_endpoints(id);
            (
                BondIdx::from(id),
                AtomIdx::from(a),
                AtomIdx::from(b),
                &self.bonds[id.index()],
            )
        })
    }
}

// ---------------------------------------------------------------------------
// Dative bond accessors
// ---------------------------------------------------------------------------

impl MoleculeAst {
    pub fn dative_bond(&self, idx: DativeBondIdx) -> &BondAst {
        self.dative_bonds.data(RelationId::from(idx))
    }

    pub fn dative_bond_participants(&self, idx: DativeBondIdx) -> &[NodeId] {
        self.dative_bonds.participants(RelationId::from(idx))
    }

    pub fn dative_bond_ids(&self) -> impl Iterator<Item = DativeBondIdx> {
        self.dative_bonds
            .relation_ids()
            .map(DativeBondIdx::from)
    }

    pub fn dative_bond_count(&self) -> usize {
        self.dative_bonds.relation_count()
    }
}

// ---------------------------------------------------------------------------
// Noncovalent bond accessors
// ---------------------------------------------------------------------------

impl MoleculeAst {
    pub fn noncovalent_bond(&self, idx: NoncovalentBondIdx) -> &BondAst {
        self.noncovalent_bonds.data(RelationId::from(idx))
    }

    pub fn noncovalent_bond_participants(&self, idx: NoncovalentBondIdx) -> &[NodeId] {
        self.noncovalent_bonds.participants(RelationId::from(idx))
    }

    pub fn noncovalent_bond_ids(&self) -> impl Iterator<Item = NoncovalentBondIdx> {
        self.noncovalent_bonds
            .relation_ids()
            .map(NoncovalentBondIdx::from)
    }

    pub fn noncovalent_bond_count(&self) -> usize {
        self.noncovalent_bonds.relation_count()
    }
}

// ---------------------------------------------------------------------------
// Aromatic system accessors
// ---------------------------------------------------------------------------

impl MoleculeAst {
    pub fn aromatic_system(&self, idx: AromaticSystemIdx) -> &AromaticSystemAst {
        self.aromatic_systems.data(RelationId::from(idx))
    }

    pub fn aromatic_system_participants(&self, idx: AromaticSystemIdx) -> &[NodeId] {
        self.aromatic_systems.participants(RelationId::from(idx))
    }

    pub fn aromatic_system_ids(&self) -> impl Iterator<Item = AromaticSystemIdx> {
        self.aromatic_systems
            .relation_ids()
            .map(AromaticSystemIdx::from)
    }

    pub fn aromatic_system_count(&self) -> usize {
        self.aromatic_systems.relation_count()
    }
}

// ---------------------------------------------------------------------------
// Multicenter bond accessors
// ---------------------------------------------------------------------------

impl MoleculeAst {
    pub fn multicenter_bond(&self, idx: MulticenterBondIdx) -> &MulticenterBondAst {
        self.multicenter_bonds.data(RelationId::from(idx))
    }

    pub fn multicenter_bond_participants(&self, idx: MulticenterBondIdx) -> &[NodeId] {
        self.multicenter_bonds.participants(RelationId::from(idx))
    }

    pub fn multicenter_bond_ids(&self) -> impl Iterator<Item = MulticenterBondIdx> {
        self.multicenter_bonds
            .relation_ids()
            .map(MulticenterBondIdx::from)
    }

    pub fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.relation_count()
    }
}

// ---------------------------------------------------------------------------
// Graph access (read-only)
// ---------------------------------------------------------------------------

impl MoleculeAst {
    pub fn graph(&self) -> &Graph {
        &self.graph
    }
}

// ---------------------------------------------------------------------------
// Derived queries
// ---------------------------------------------------------------------------

impl MoleculeAst {
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
        let node = NodeId::from(atom);
        let mut sum: u8 = 0;
        for neighbor in self.graph.neighbors(node) {
            match self.bonds[neighbor.edge.index()].order {
                ValueAst::Lit(n) => sum += n as u8,
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
                ValueAst::Lit(n) => n as u8,
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

impl Ast for MoleculeAst {
    type Config = MoleculeAstConfig;
}


/// A `MoleculeAst` whose fields are all concrete and whose constraints are
/// all ground assertions. The invariant is checked once at construction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GroundMolecule(MoleculeAst);

impl GroundMolecule {
    pub fn new(ast: MoleculeAst) -> Result<Self, GroundError> {
        if ast.is_ground() {
            Ok(Self(ast))
        } else {
            Err(GroundError)
        }
    }

    pub fn as_ast(&self) -> &MoleculeAst {
        &self.0
    }

    pub fn into_ast(self) -> MoleculeAst {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use umol_shared::atom_ast::{AromaticValenceAst, ElementAst, HydrogenAst, IsotopeAst};
    use umol_shared::spin::SpinState;
    use umol_shared::spin_ast::SpinStateAst;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use crate::ast::constraint::{DerivedPred, RelationRefs};

    fn ground_atom() -> AtomAst {
        AtomAst {
            element: ElementAst::Lit(umol_shared::element::Element::C),
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
                AtomAst::from_element(umol_shared::element::Element::C),
                AtomAst::from_element(umol_shared::element::Element::O),
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
    fn test_ground_molecule_new() {
        let ast = ground_ast();
        assert!(GroundMolecule::new(ast).is_ok());
    }

    #[test]
    fn test_ground_molecule_new_error() {
        let ast = MoleculeAst::new(
            vec![AtomAst::new(ElementAst::Undetermined)],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        assert_eq!(GroundMolecule::new(ast), Err(GroundError));
    }

    #[test]
    fn test_molecule_ast_bond_order_sum() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(umol_shared::element::Element::C),
                AtomAst::from_element(umol_shared::element::Element::O),
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
            vec![AtomAst::from_element(umol_shared::element::Element::C)],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        assert_eq!(ast.bond_order_sum(AtomIdx(0)), Some(0));
    }

    #[test]
    fn test_molecule_ast_bond_order_sum_wildcard() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(umol_shared::element::Element::C),
                AtomAst::from_element(umol_shared::element::Element::O),
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
                AtomAst::from_element(umol_shared::element::Element::C),
                AtomAst::from_element(umol_shared::element::Element::C),
                AtomAst::from_element(umol_shared::element::Element::C),
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
                AtomAst::from_element(umol_shared::element::Element::N),
                AtomAst::from_element(umol_shared::element::Element::B),
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
                AtomAst::from_element(umol_shared::element::Element::C),
                AtomAst::from_element(umol_shared::element::Element::O),
                AtomAst::from_element(umol_shared::element::Element::N),
            ],
            vec![
                (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
                (AtomIdx(0), AtomIdx(2), BondAst::from_order(2)),
            ],
            vec![], vec![], vec![], vec![], vec![],
        );
        assert_eq!(ast.neighbors(AtomIdx(0)).len(), 2);
        assert_eq!(ast.neighbors(AtomIdx(1)).len(), 1);
        assert_eq!(ast.neighbors(AtomIdx(2)).len(), 1);
    }
}
