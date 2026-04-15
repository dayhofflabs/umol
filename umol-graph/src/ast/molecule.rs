//! Molecule structural AST.

use umol_graph_core::relation::RelationId;
use umol_graph_core::{EdgeId, Graph, Neighbor, NodeId, RelationSet};
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
/// Primary topology (atoms + localized bonds) lives in `Graph<AtomAst, BondAst>`.
/// Secondary relations (dative bonds, noncovalent bonds, aromatic systems,
/// multicenter bonds) live in typed `RelationSet`s sharing the same `NodeId` space.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MoleculeAst {
    graph: Graph<AtomAst, BondAst>,
    dative_bonds: RelationSet<BondAst>,
    noncovalent_bonds: RelationSet<BondAst>,
    aromatic_systems: RelationSet<AromaticSystemAst>,
    multicenter_bonds: RelationSet<MulticenterBondAst>,
    pub constraints: Vec<MoleculeConstraint>,
}

// ---------------------------------------------------------------------------
// Atom accessors
// ---------------------------------------------------------------------------

impl MoleculeAst {
    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    pub fn atom(&self, idx: AtomIdx) -> &AtomAst {
        &self.graph[NodeId::from(idx)]
    }

    pub fn atom_mut(&mut self, idx: AtomIdx) -> &mut AtomAst {
        &mut self.graph[NodeId::from(idx)]
    }

    pub fn atoms(&self) -> impl Iterator<Item = (AtomIdx, &AtomAst)> {
        self.graph
            .node_ids()
            .map(|id| (AtomIdx::from(id), &self.graph[id]))
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
        self.graph.edge_count()
    }

    pub fn bond(&self, idx: BondIdx) -> &BondAst {
        &self.graph[EdgeId::from(idx)]
    }

    pub fn bond_endpoints(&self, idx: BondIdx) -> (AtomIdx, AtomIdx) {
        let [a, b] = self.graph.edge_endpoints(EdgeId::from(idx)).unwrap();
        (AtomIdx::from(a), AtomIdx::from(b))
    }

    pub fn bonds(&self) -> impl Iterator<Item = (BondIdx, AtomIdx, AtomIdx, &BondAst)> {
        self.graph.edge_ids().map(|id| {
            let [a, b] = self.graph.edge_endpoints(id).unwrap();
            (
                BondIdx::from(id),
                AtomIdx::from(a),
                AtomIdx::from(b),
                &self.graph[id],
            )
        })
    }
}

// ---------------------------------------------------------------------------
// Dative bond accessors
// ---------------------------------------------------------------------------

impl MoleculeAst {
    pub fn dative_bond(&self, idx: DativeBondIdx) -> &BondAst {
        self.dative_bonds
            .data(RelationId::from(idx))
            .expect("invalid DativeBondIdx")
    }

    pub fn dative_bond_participants(&self, idx: DativeBondIdx) -> &[NodeId] {
        self.dative_bonds
            .participants(RelationId::from(idx))
            .expect("invalid DativeBondIdx")
    }

    pub fn dative_bond_ids(&self) -> impl Iterator<Item = DativeBondIdx> + '_ {
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
        self.noncovalent_bonds
            .data(RelationId::from(idx))
            .expect("invalid NoncovalentBondIdx")
    }

    pub fn noncovalent_bond_participants(&self, idx: NoncovalentBondIdx) -> &[NodeId] {
        self.noncovalent_bonds
            .participants(RelationId::from(idx))
            .expect("invalid NoncovalentBondIdx")
    }

    pub fn noncovalent_bond_ids(&self) -> impl Iterator<Item = NoncovalentBondIdx> + '_ {
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
        self.aromatic_systems
            .data(RelationId::from(idx))
            .expect("invalid AromaticSystemIdx")
    }

    pub fn aromatic_system_participants(&self, idx: AromaticSystemIdx) -> &[NodeId] {
        self.aromatic_systems
            .participants(RelationId::from(idx))
            .expect("invalid AromaticSystemIdx")
    }

    pub fn aromatic_system_ids(&self) -> impl Iterator<Item = AromaticSystemIdx> + '_ {
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
        self.multicenter_bonds
            .data(RelationId::from(idx))
            .expect("invalid MulticenterBondIdx")
    }

    pub fn multicenter_bond_participants(&self, idx: MulticenterBondIdx) -> &[NodeId] {
        self.multicenter_bonds
            .participants(RelationId::from(idx))
            .expect("invalid MulticenterBondIdx")
    }

    pub fn multicenter_bond_ids(&self) -> impl Iterator<Item = MulticenterBondIdx> + '_ {
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
    pub fn graph(&self) -> &Graph<AtomAst, BondAst> {
        &self.graph
    }
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

impl MoleculeAst {
    pub fn add_atom(&mut self, atom: AtomAst) -> AtomIdx {
        AtomIdx::from(self.graph.add_node(atom))
    }

    pub fn add_bond(
        &mut self,
        source: AtomIdx,
        target: AtomIdx,
        bond: BondAst,
    ) -> BondIdx {
        BondIdx::from(
            self.graph
                .add_edge(NodeId::from(source), NodeId::from(target), bond),
        )
    }

    pub fn add_dative_bond(
        &mut self,
        donor: AtomIdx,
        acceptor: AtomIdx,
        bond: BondAst,
    ) -> DativeBondIdx {
        let id = self.dative_bonds.add(
            vec![NodeId::from(donor), NodeId::from(acceptor)],
            bond,
        );
        DativeBondIdx::from(id)
    }

    pub fn add_noncovalent_bond(
        &mut self,
        a: AtomIdx,
        b: AtomIdx,
        bond: BondAst,
    ) -> NoncovalentBondIdx {
        let id = self.noncovalent_bonds.add(
            vec![NodeId::from(a), NodeId::from(b)],
            bond,
        );
        NoncovalentBondIdx::from(id)
    }

    pub fn add_aromatic_system(
        &mut self,
        atoms: Vec<AtomIdx>,
        data: AromaticSystemAst,
    ) -> AromaticSystemIdx {
        let participants: Vec<NodeId> = atoms.into_iter().map(NodeId::from).collect();
        AromaticSystemIdx::from(self.aromatic_systems.add(participants, data))
    }

    pub fn add_multicenter_bond(
        &mut self,
        atoms: Vec<AtomIdx>,
        data: MulticenterBondAst,
    ) -> MulticenterBondIdx {
        let participants: Vec<NodeId> = atoms.into_iter().map(NodeId::from).collect();
        MulticenterBondIdx::from(self.multicenter_bonds.add(participants, data))
    }
}

// ---------------------------------------------------------------------------
// Derived queries
// ---------------------------------------------------------------------------

impl MoleculeAst {
    pub fn is_ground(&self) -> bool {
        self.graph.node_ids().all(|id| self.graph[id].is_ground())
            && self.graph.edge_ids().all(|id| self.graph[id].is_ground())
            && self
                .dative_bonds
                .relation_ids()
                .all(|id| self.dative_bonds.data(id).unwrap().is_ground())
            && self
                .noncovalent_bonds
                .relation_ids()
                .all(|id| self.noncovalent_bonds.data(id).unwrap().is_ground())
            && self.constraints.iter().all(|c| c.is_ground_assertion())
    }

    pub fn bond_order_sum(&self, atom: AtomIdx) -> Option<u8> {
        let node = NodeId::from(atom);
        let mut sum: u8 = 0;
        for neighbor in self.graph.neighbors(node) {
            match self.graph[neighbor.edge].order {
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
            let data = self.dative_bonds.data(rel_id).unwrap();
            let order = match data.order {
                ValueAst::Lit(n) => n as u8,
                _ => continue,
            };
            let participants = self.dative_bonds.participants(rel_id).unwrap();
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
        let mut ast = MoleculeAst::default();
        ast.add_atom(ground_atom());
        ast
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
        let mut ast = MoleculeAst::default();
        ast.add_atom(AtomAst::new(ElementAst::Undetermined));
        assert!(!ast.is_ground());
    }

    #[test]
    fn test_molecule_ast_is_ground_wildcard_bond() {
        let mut ast = MoleculeAst::default();
        let a = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::C));
        let b = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::O));
        ast.add_bond(a, b, BondAst::new(ValueAst::Undetermined));
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
        let mut ast = MoleculeAst::default();
        ast.add_atom(AtomAst::new(ElementAst::Undetermined));
        assert_eq!(GroundMolecule::new(ast), Err(GroundError));
    }

    #[test]
    fn test_molecule_ast_bond_order_sum() {
        let mut ast = MoleculeAst::default();
        let a = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::C));
        let b = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::O));
        ast.add_bond(a, b, BondAst::from_order(2));
        assert_eq!(ast.bond_order_sum(a), Some(2));
        assert_eq!(ast.bond_order_sum(b), Some(2));
    }

    #[test]
    fn test_molecule_ast_bond_order_sum_no_bonds() {
        let mut ast = MoleculeAst::default();
        let a = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::C));
        assert_eq!(ast.bond_order_sum(a), Some(0));
    }

    #[test]
    fn test_molecule_ast_bond_order_sum_wildcard() {
        let mut ast = MoleculeAst::default();
        let a = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::C));
        let b = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::O));
        ast.add_bond(a, b, BondAst::new(ValueAst::Undetermined));
        assert_eq!(ast.bond_order_sum(a), None);
    }

    #[test]
    fn test_molecule_ast_is_in_aromatic_system() {
        let mut ast = MoleculeAst::default();
        let a = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::C));
        let b = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::C));
        let c = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::C));
        ast.add_aromatic_system(vec![a, b], AromaticSystemAst {});
        assert!(ast.is_in_aromatic_system(a));
        assert!(ast.is_in_aromatic_system(b));
        assert!(!ast.is_in_aromatic_system(c));
    }

    #[test]
    fn test_molecule_ast_dative_bond_order_sums() {
        let mut ast = MoleculeAst::default();
        let n = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::N));
        let b = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::B));
        ast.add_dative_bond(n, b, BondAst::from_order(1));
        assert_eq!(ast.dative_bond_order_sums(n), (1, 0));
        assert_eq!(ast.dative_bond_order_sums(b), (0, 1));
    }

    #[test]
    fn test_molecule_ast_neighbors() {
        let mut ast = MoleculeAst::default();
        let a = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::C));
        let b = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::O));
        let c = ast.add_atom(AtomAst::from_element(umol_shared::element::Element::N));
        ast.add_bond(a, b, BondAst::from_order(1));
        ast.add_bond(a, c, BondAst::from_order(2));
        assert_eq!(ast.neighbors(a).len(), 2);
        assert_eq!(ast.neighbors(b).len(), 1);
        assert_eq!(ast.neighbors(c).len(), 1);
    }
}
