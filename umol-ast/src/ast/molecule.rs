//! Molecule structural AST.

mod builder;

use std::collections::HashSet;
use std::ops::Index;
use std::sync::Arc;

pub use builder::MoleculeBuilder;
use umol_graph_core::relation::RelationId;
use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm,
    CycleEnumerationAlgorithm, EdgeId, FixedRelationSet, Graph, MaxIndependentSetAlgorithm,
    MaxMatchingAlgorithm, MatchingEnumerationAlgorithm, NodeId, ShortestCycleAlgorithm,
    SubgraphIsomorphismAlgorithm, VarRelationSet,
};

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::automorphism::AtomAutomorphism;
use super::bond::BondAst;
use super::constraint::Constraints;
use super::dative::DativeBondAst;
use super::error::RewriteError;
use super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use super::matching::BondMatching;
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;
use super::rings::{self, RingFamily, RingSet};
use super::subgraph::MoleculeSubgraph;
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

    // -- Ring enumeration -----------------------------------------------------

    pub fn rings(
        &self,
        family: RingFamily,
        max_ring_size: usize,
        atom_filter: impl Fn(AtomIdx) -> bool,
    ) -> RingSet {
        rings::enumerate_rings(&self.graph, family, max_ring_size, atom_filter)
    }

    // -- Read: induced subsets ----------------------------------------------

    pub fn degree(&self, atom: AtomIdx) -> usize {
        self.graph.degree(NodeId::from(atom))
    }

    pub fn connected_components(
        &self,
        alg: ConnectedComponentsAlgorithm,
    ) -> Vec<Vec<AtomIdx>> {
        self.graph
            .connected_components(alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomIdx::from).collect())
            .collect()
    }

    pub fn biconnected_components(
        &self,
        alg: BiconnectedComponentsAlgorithm,
    ) -> Vec<Vec<AtomIdx>> {
        self.graph
            .biconnected_components(alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomIdx::from).collect())
            .collect()
    }

    pub fn shortest_cycle_through_bond(
        &self,
        bond: BondIdx,
        alg: ShortestCycleAlgorithm,
    ) -> Option<usize> {
        self.graph
            .shortest_cycle_through_edge(EdgeId::from(bond), alg)
    }

    pub fn shortest_cycle_through_atom(
        &self,
        atom: AtomIdx,
        alg: ShortestCycleAlgorithm,
    ) -> Option<usize> {
        self.graph
            .shortest_cycle_through_node(NodeId::from(atom), alg)
    }

    pub fn enumerate_cycles(
        &self,
        max_size: usize,
        alg: CycleEnumerationAlgorithm,
    ) -> Vec<Vec<AtomIdx>> {
        self.graph
            .enumerate_cycles(max_size, alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomIdx::from).collect())
            .collect()
    }

    pub fn maximum_independent_set(&self, alg: MaxIndependentSetAlgorithm) -> Vec<AtomIdx> {
        self.graph
            .maximum_independent_set(alg)
            .into_iter()
            .map(AtomIdx::from)
            .collect()
    }

    pub fn maximum_matching(&self, alg: MaxMatchingAlgorithm) -> BondMatching {
        BondMatching(self.graph.maximum_matching(alg))
    }

    pub fn enumerate_perfect_matchings(
        &self,
        alg: MatchingEnumerationAlgorithm,
    ) -> Vec<BondMatching> {
        self.graph
            .enumerate_perfect_matchings(alg)
            .into_iter()
            .map(BondMatching)
            .collect()
    }

    pub fn enumerate_maximum_matchings(
        &self,
        alg: MatchingEnumerationAlgorithm,
    ) -> Vec<BondMatching> {
        self.graph
            .enumerate_maximum_matchings(alg)
            .into_iter()
            .map(BondMatching)
            .collect()
    }

    pub fn automorphisms<C: Ord + Copy>(
        &self,
        atom_color: impl Fn(AtomIdx) -> C,
        alg: AutomorphismAlgorithm,
    ) -> AtomAutomorphism {
        AtomAutomorphism(self.graph.automorphisms(|n| atom_color(AtomIdx::from(n)), alg))
    }

    pub fn subgraph_isomorphisms(
        &self,
        query: &MoleculeAst,
        atom_match: &mut impl FnMut(AtomIdx, AtomIdx) -> bool,
        bond_match: &mut impl FnMut(BondIdx, BondIdx) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<AtomIdx>> {
        self.graph
            .subgraph_isomorphisms(
                &query.graph,
                &mut |tn, qn| atom_match(AtomIdx::from(tn), AtomIdx::from(qn)),
                &mut |te, qe| bond_match(BondIdx::from(te), BondIdx::from(qe)),
                alg,
            )
            .into_iter()
            .map(|m| m.into_iter().map(AtomIdx::from).collect())
            .collect()
    }

    pub fn subgraph_isomorphisms_at(
        &self,
        query: &MoleculeAst,
        anchor: (AtomIdx, AtomIdx),
        atom_match: &mut impl FnMut(AtomIdx, AtomIdx) -> bool,
        bond_match: &mut impl FnMut(BondIdx, BondIdx) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<AtomIdx>> {
        self.graph
            .subgraph_isomorphisms_at(
                &query.graph,
                (NodeId::from(anchor.0), NodeId::from(anchor.1)),
                &mut |tn, qn| atom_match(AtomIdx::from(tn), AtomIdx::from(qn)),
                &mut |te, qe| bond_match(BondIdx::from(te), BondIdx::from(qe)),
                alg,
            )
            .into_iter()
            .map(|m| m.into_iter().map(AtomIdx::from).collect())
            .collect()
    }

    // -- Read: induced subsets (bond/relation) --------------------------------

    pub fn induced_subgraph(&self, atoms: &[AtomIdx]) -> MoleculeSubgraph {
        let keep: HashSet<AtomIdx> = atoms.iter().copied().collect();
        let remove_atoms: Vec<AtomIdx> = (0..self.atom_count())
            .map(AtomIdx::from)
            .filter(|a| !keep.contains(a))
            .collect();
        let remove_bonds: Vec<BondIdx> = self
            .bonds()
            .iter()
            .filter(|b| !keep.contains(&b.src) || !keep.contains(&b.tgt))
            .map(|b| b.idx)
            .collect();
        let mut builder = self.edit();
        let remap = builder.remove(&remove_atoms, &remove_bonds);
        let ast = builder.build();

        let atom_map: Vec<AtomIdx> = (0..self.atom_count())
            .map(AtomIdx::from)
            .filter(|&a| remap.atom(a).is_some())
            .collect();
        let bond_map: Vec<BondIdx> = (0..self.bond_count())
            .map(BondIdx::from)
            .filter(|&b| remap.bond(b).is_some())
            .collect();
        let dative_bond_map: Vec<DativeBondIdx> = (0..self.dative_bond_count())
            .map(DativeBondIdx::from)
            .filter(|&d| remap.dative_bond(d).is_some())
            .collect();
        let aromatic_system_map: Vec<AromaticSystemIdx> = (0..self.aromatic_system_count())
            .map(AromaticSystemIdx::from)
            .filter(|&a| remap.aromatic_system(a).is_some())
            .collect();
        let multicenter_bond_map: Vec<MulticenterBondIdx> = (0..self.multicenter_bond_count())
            .map(MulticenterBondIdx::from)
            .filter(|&m| remap.multicenter_bond(m).is_some())
            .collect();
        let noncovalent_bond_map: Vec<NoncovalentBondIdx> = (0..self.noncovalent_bond_count())
            .map(NoncovalentBondIdx::from)
            .filter(|&n| remap.noncovalent_bond(n).is_some())
            .collect();

        MoleculeSubgraph {
            ast,
            atom_map,
            bond_map,
            dative_bond_map,
            aromatic_system_map,
            multicenter_bond_map,
            noncovalent_bond_map,
        }
    }

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

// -- DPO rewriting -----------------------------------------------------------

impl MoleculeAst {
    /// Apply a DPO reaction rule to this molecule given a match assignment.
    ///
    /// Phases: (1) add R\K entities, (2) modify K attributes, (3) remove
    /// L\K overlay relations, (4) remove L\K atoms/bonds. See discussion
    /// doc 90 for the full case analysis.
    pub fn apply_rule(
        &self,
        _rule: &super::reaction::ReactionRuleAst,
        _assignment: &super::reaction::Assignment,
    ) -> Result<MoleculeAst, RewriteError> {
        todo!("DPO rewrite implementation — see discussion/90-reactions-relation-mutation-2026-04-21.md")
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
    use rstest::*;
    use umol_graph_core::{
        BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm,
        CycleEnumerationAlgorithm, MaxIndependentSetAlgorithm, MaxMatchingAlgorithm,
        MatchingEnumerationAlgorithm, ShortestCycleAlgorithm,
    };
    use umol_shared::element::Element;

    use super::super::atom::ElementAst;
    use super::super::constraint::{Constraint, MoleculeConstraint};
    use super::super::dative::DativeBondAst;
    use super::super::multicenter::MulticenterBondAst;
    use super::super::noncovalent::{NoncovalentBondAst, NoncovalentKind};
    use super::super::rings::RingFamily;
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

    fn chain(n: usize) -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); n];
        let bonds: Vec<_> = (0..n.saturating_sub(1))
            .map(|i| {
                (
                    AtomIdx(i as u32),
                    AtomIdx((i + 1) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], Constraints::default())
    }

    fn ring(n: usize) -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); n];
        let bonds: Vec<_> = (0..n)
            .map(|i| {
                (
                    AtomIdx(i as u32),
                    AtomIdx(((i + 1) % n) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], Constraints::default())
    }

    fn two_components() -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); 4];
        let bonds = vec![
            (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
            (AtomIdx(2), AtomIdx(3), BondAst::from_order(1)),
        ];
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], Constraints::default())
    }

    #[rstest]
    #[case::isolated(chain(1), AtomIdx(0), 0)]
    #[case::chain_end(chain(3), AtomIdx(0), 1)]
    #[case::chain_mid(chain(3), AtomIdx(1), 2)]
    #[case::ring_vertex(ring(6), AtomIdx(0), 2)]
    fn test_molecule_ast_degree(
        #[case] ast: MoleculeAst,
        #[case] atom: AtomIdx,
        #[case] expected: usize,
    ) {
        assert_eq!(ast.degree(atom), expected);
    }

    #[rstest]
    #[case::single(chain(3), 1)]
    #[case::two(two_components(), 2)]
    #[case::empty(MoleculeAst::default(), 0)]
    fn test_molecule_ast_connected_components(#[case] ast: MoleculeAst, #[case] expected: usize) {
        let cc = ast.connected_components(ConnectedComponentsAlgorithm::Bfs);
        assert_eq!(cc.len(), expected);
    }

    #[rstest]
    #[case::ring_6(ring(6), 1)]
    #[case::chain(chain(5), 0)]
    fn test_molecule_ast_biconnected_components(
        #[case] ast: MoleculeAst,
        #[case] expected: usize,
    ) {
        let bcc = ast.biconnected_components(BiconnectedComponentsAlgorithm::Tarjan);
        assert_eq!(bcc.len(), expected);
    }

    #[rstest]
    #[case::ring_bond(ring(6), BondIdx(0), Some(6))]
    #[case::chain_bond(chain(3), BondIdx(0), None)]
    fn test_molecule_ast_shortest_cycle_through_bond(
        #[case] ast: MoleculeAst,
        #[case] bond: BondIdx,
        #[case] expected: Option<usize>,
    ) {
        assert_eq!(
            ast.shortest_cycle_through_bond(bond, ShortestCycleAlgorithm::Bfs),
            expected
        );
    }

    #[rstest]
    #[case::ring_atom(ring(6), AtomIdx(0), Some(6))]
    #[case::chain_atom(chain(3), AtomIdx(1), None)]
    fn test_molecule_ast_shortest_cycle_through_atom(
        #[case] ast: MoleculeAst,
        #[case] atom: AtomIdx,
        #[case] expected: Option<usize>,
    ) {
        assert_eq!(
            ast.shortest_cycle_through_atom(atom, ShortestCycleAlgorithm::Bfs),
            expected
        );
    }

    #[rstest]
    #[case::hexagon(ring(6), 6, 1)]
    #[case::hexagon_cutoff(ring(6), 5, 0)]
    #[case::chain(chain(5), 10, 0)]
    #[case::empty(MoleculeAst::default(), 10, 0)]
    fn test_molecule_ast_enumerate_cycles(
        #[case] ast: MoleculeAst,
        #[case] max_size: usize,
        #[case] expected: usize,
    ) {
        let cycles = ast.enumerate_cycles(max_size, CycleEnumerationAlgorithm::Vismara);
        assert_eq!(cycles.len(), expected);
    }

    #[rstest]
    #[case::triangle(ring(3), 1)]
    #[case::chain_3(chain(3), 2)]
    fn test_molecule_ast_maximum_independent_set(
        #[case] ast: MoleculeAst,
        #[case] expected: usize,
    ) {
        let mis = ast.maximum_independent_set(MaxIndependentSetAlgorithm::BranchAndBound);
        assert_eq!(mis.len(), expected);
    }

    #[rstest]
    #[case::chain_4(chain(4), 2)]
    #[case::ring_6(ring(6), 3)]
    #[case::single(chain(1), 0)]
    fn test_molecule_ast_maximum_matching(
        #[case] ast: MoleculeAst,
        #[case] expected_size: usize,
    ) {
        let m = ast.maximum_matching(MaxMatchingAlgorithm::Edmonds);
        assert_eq!(m.size(), expected_size);
    }

    #[test]
    fn test_bond_matching_mate() {
        let ast = chain(4);
        let m = ast.maximum_matching(MaxMatchingAlgorithm::Edmonds);
        assert!(m.is_matched(AtomIdx(0)));
        let mate = m.mate(AtomIdx(0));
        assert!(mate.is_some());
    }

    #[rstest]
    #[case::ring_6(ring(6), 2)]
    fn test_molecule_ast_enumerate_perfect_matchings(
        #[case] ast: MoleculeAst,
        #[case] expected: usize,
    ) {
        let ms = ast.enumerate_perfect_matchings(MatchingEnumerationAlgorithm::BranchAndBound);
        assert_eq!(ms.len(), expected);
        for m in &ms {
            assert!(m.is_perfect(ast.atom_count()));
        }
    }

    #[rstest]
    #[case::ring_6(ring(6), 1)]
    #[case::chain_3(chain(3), 2)]
    fn test_molecule_ast_automorphisms(
        #[case] ast: MoleculeAst,
        #[case] expected_orbits: usize,
    ) {
        let auto = ast.automorphisms(
            |_| 0u8,
            umol_graph_core::AutomorphismAlgorithm::Nauty,
        );
        assert_eq!(auto.num_orbits(), expected_orbits);
        assert_eq!(auto.atom_count(), ast.atom_count());
    }

    #[test]
    fn test_atom_automorphism_same_orbit() {
        let ast = ring(6);
        let auto = ast.automorphisms(
            |_| 0u8,
            umol_graph_core::AutomorphismAlgorithm::Nauty,
        );
        assert!(auto.same_orbit(AtomIdx(0), AtomIdx(3)));
    }

    #[test]
    fn test_molecule_ast_subgraph_isomorphisms() {
        let target = ring(6);
        let query = chain(2);
        let matches = target.subgraph_isomorphisms(
            &query,
            &mut |_, _| true,
            &mut |_, _| true,
            umol_graph_core::SubgraphIsomorphismAlgorithm::Vf2,
        );
        assert_eq!(matches.len(), 12);
    }

    #[test]
    fn test_molecule_ast_subgraph_isomorphisms_at() {
        let target = ring(6);
        let query = chain(2);
        let matches = target.subgraph_isomorphisms_at(
            &query,
            (AtomIdx(0), AtomIdx(0)),
            &mut |_, _| true,
            &mut |_, _| true,
            umol_graph_core::SubgraphIsomorphismAlgorithm::Vf2,
        );
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn test_molecule_ast_induced_subgraph() {
        let ast = rich_molecule();
        let sub = ast.induced_subgraph(&[AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
        assert_eq!(sub.ast.atom_count(), 3);
        assert_eq!(sub.ast.bond_count(), 2);
        assert_eq!(sub.atom_map, vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
        assert_eq!(sub.bond_map, vec![BondIdx(0), BondIdx(1)]);
        assert_eq!(sub.ast.aromatic_system_count(), 1);
        assert_eq!(sub.aromatic_system_map, vec![AromaticSystemIdx(0)]);
        assert_eq!(sub.ast.multicenter_bond_count(), 1);
        assert_eq!(sub.multicenter_bond_map, vec![MulticenterBondIdx(0)]);
        assert_eq!(sub.ast.dative_bond_count(), 0);
        assert!(sub.dative_bond_map.is_empty());
        assert_eq!(sub.ast.noncovalent_bond_count(), 0);
        assert!(sub.noncovalent_bond_map.is_empty());
    }

    #[test]
    fn test_molecule_ast_induced_subgraph_preserves_dative() {
        let ast = rich_molecule();
        let sub = ast.induced_subgraph(&[AtomIdx(2), AtomIdx(3)]);
        assert_eq!(sub.ast.atom_count(), 2);
        assert_eq!(sub.ast.dative_bond_count(), 1);
        assert_eq!(sub.dative_bond_map, vec![DativeBondIdx(0)]);
    }

    #[test]
    fn test_builder_remove_aromatic_systems() {
        let ast = rich_molecule();
        let mut b = ast.edit();
        b.remove_aromatic_systems(&[AromaticSystemIdx(0)]);
        let result = b.build();
        assert_eq!(result.aromatic_system_count(), 0);
        assert_eq!(result.atom_count(), 4);
        assert_eq!(result.bond_count(), 3);
    }

    #[test]
    fn test_builder_remove_dative_bonds() {
        let ast = rich_molecule();
        let mut b = ast.edit();
        b.remove_dative_bonds(&[DativeBondIdx(0)]);
        let result = b.build();
        assert_eq!(result.dative_bond_count(), 0);
        assert_eq!(result.atom_count(), 4);
    }

    #[test]
    fn test_builder_remove_multicenter_bonds() {
        let ast = rich_molecule();
        let mut b = ast.edit();
        b.remove_multicenter_bonds(&[MulticenterBondIdx(0)]);
        let result = b.build();
        assert_eq!(result.multicenter_bond_count(), 0);
    }

    #[test]
    fn test_builder_remove_noncovalent_bonds() {
        let ast = rich_molecule();
        let mut b = ast.edit();
        b.remove_noncovalent_bonds(&[NoncovalentBondIdx(0)]);
        let result = b.build();
        assert_eq!(result.noncovalent_bond_count(), 0);
    }

    #[test]
    fn test_builder_atom_mut() {
        let ast = rich_molecule();
        let mut b = ast.edit();
        b.atom_mut(AtomIdx(0)).element = ElementAst::Lit(Element::N);
        let result = b.build();
        assert_eq!(result[AtomIdx(0)].element, ElementAst::Lit(Element::N));
        assert_eq!(ast[AtomIdx(0)].element, ElementAst::Lit(Element::C));
    }

    #[test]
    fn test_builder_bond_mut() {
        let ast = rich_molecule();
        let mut b = ast.edit();
        b.bond_mut(BondIdx(0)).order = ValueAst::Lit(3);
        let result = b.build();
        assert_eq!(result[BondIdx(0)].order, ValueAst::Lit(3));
        assert_eq!(ast[BondIdx(0)].order, ValueAst::Lit(1));
    }

    #[test]
    fn test_builder_constraints_mut() {
        let ast = rich_molecule();
        let mut b = ast.edit();
        b.constraints_mut().push_atom(AtomIdx(0), super::super::constraint::AtomConstraint::Degree(ValueAst::Lit(2)));
        let result = b.build();
        assert_eq!(result.constraints().atom(AtomIdx(0)).len(), 1);
        assert!(ast.constraints().atom(AtomIdx(0)).is_empty());
    }

    #[rstest]
    #[case::hexagon(ring(6), 6, 1)]
    #[case::hexagon_cutoff(ring(6), 5, 0)]
    #[case::chain(chain(5), 10, 0)]
    #[case::empty(MoleculeAst::default(), 10, 0)]
    fn test_molecule_ast_rings(
        #[case] ast: MoleculeAst,
        #[case] max_ring_size: usize,
        #[case] expected: usize,
    ) {
        let rs = ast.rings(RingFamily::Simple, max_ring_size, |_| true);
        assert_eq!(rs.count(), expected);
    }

    #[test]
    fn test_molecule_ast_rings_atom_filter() {
        let ast = ring(6);
        let rs = ast.rings(RingFamily::Simple, 10, |a| a.0 < 3);
        assert_eq!(rs.count(), 0);
    }

    #[test]
    fn test_molecule_ast_rings_induced() {
        // K4 = complete graph on 4 nodes (6 edges, all pairs connected)
        let atoms = vec![AtomAst::from_element(Element::C); 4];
        let bonds = vec![
            (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
            (AtomIdx(0), AtomIdx(2), BondAst::from_order(1)),
            (AtomIdx(0), AtomIdx(3), BondAst::from_order(1)),
            (AtomIdx(1), AtomIdx(2), BondAst::from_order(1)),
            (AtomIdx(1), AtomIdx(3), BondAst::from_order(1)),
            (AtomIdx(2), AtomIdx(3), BondAst::from_order(1)),
        ];
        let ast = MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], Constraints::default());
        let simple = ast.rings(RingFamily::Simple, 4, |_| true);
        let induced = ast.rings(RingFamily::Induced, 4, |_| true);
        // K4 has 4 relevant triangles; all are induced (no chords in a triangle)
        assert_eq!(simple.count(), 4);
        assert_eq!(induced.count(), 4);
    }

    #[test]
    fn test_molecule_ast_rings_induced_naphthalene() {
        let atoms = vec![AtomAst::from_element(Element::C); 10];
        #[rustfmt::skip]
        let bonds = vec![
            (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
            (AtomIdx(1), AtomIdx(2), BondAst::from_order(1)),
            (AtomIdx(2), AtomIdx(3), BondAst::from_order(1)),
            (AtomIdx(3), AtomIdx(4), BondAst::from_order(1)),
            (AtomIdx(4), AtomIdx(5), BondAst::from_order(1)),
            (AtomIdx(5), AtomIdx(0), BondAst::from_order(1)),
            (AtomIdx(3), AtomIdx(6), BondAst::from_order(1)),
            (AtomIdx(6), AtomIdx(7), BondAst::from_order(1)),
            (AtomIdx(7), AtomIdx(8), BondAst::from_order(1)),
            (AtomIdx(8), AtomIdx(9), BondAst::from_order(1)),
            (AtomIdx(9), AtomIdx(4), BondAst::from_order(1)),
        ];
        let ast = MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], Constraints::default());
        let simple = ast.rings(RingFamily::Simple, 10, |_| true);
        assert_eq!(simple.count(), 2);
        let induced = ast.rings(RingFamily::Induced, 10, |_| true);
        assert_eq!(induced.count(), 2);
    }

    #[test]
    fn test_rings_membership() {
        let ast = ring(6);
        let rs = ast.rings(RingFamily::Simple, 6, |_| true);
        assert!(rs.contains_atom(AtomIdx(0)));
        assert!(rs.contains_bond(BondIdx(0)));
        assert_eq!(rs.atom_smallest_ring_size(AtomIdx(0)), Some(6));
    }

    #[test]
    fn test_dpo_add_then_remove() {
        let ast = rich_molecule();
        let mut b = ast.edit();
        let new_a = b.add_atom(AtomAst::from_element(Element::Br));
        b.add_bond(AtomIdx(0), new_a, BondAst::from_order(1));
        b.remove_aromatic_systems(&[AromaticSystemIdx(0)]);
        let _remap = b.remove(&[AtomIdx(3)], &[BondIdx(2)]);
        let result = b.build();
        assert_eq!(result.atom_count(), 4);
        assert_eq!(result.bond_count(), 3);
        assert_eq!(result.aromatic_system_count(), 0);
        assert_eq!(result.dative_bond_count(), 0);
        assert_eq!(result.noncovalent_bond_count(), 0);
    }
}
