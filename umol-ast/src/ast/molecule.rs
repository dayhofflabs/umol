//! Molecule structural AST.

use std::collections::HashSet;
use std::ops::Index;
use std::sync::Arc;

pub use builder::MoleculeBuilder;
use umol_graph_core::relation::RelationId;
use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm,
    CycleEnumerationAlgorithm, EdgeId, FixedRelationSet, Graph, MatchingEnumerationAlgorithm,
    MaxIndependentSetAlgorithm, MaxMatchingAlgorithm, NodeId, ShortestCycleAlgorithm,
    SubgraphIsomorphismAlgorithm, VarRelationSet,
};

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::automorphism::AtomAutomorphism;
use super::bond::BondAst;
use super::constraint::{Constraint, Constraints};
use super::dative::{DativeBondAst, DativeBondDirection};
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

mod builder;
mod rewrite;

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
                .map(|(donor, acceptor, mut d)| {
                    d.direction = if donor.0 <= acceptor.0 {
                        DativeBondDirection::Forward
                    } else {
                        DativeBondDirection::Reverse
                    };
                    ([NodeId::from(donor), NodeId::from(acceptor)], d)
                })
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

    // region: Read: topology

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

    // endregion: Read: topology

    // region: Read: atoms

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn atoms(&self) -> AtomViews<'_> {
        AtomViews::new(&self.atoms)
    }

    pub fn atom(&self, idx: AtomIdx) -> AtomView<'_> {
        self.atoms().get(idx)
    }

    // endregion: Read: atoms

    // region: Read: bonds

    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    pub fn bonds(&self) -> BondViews<'_> {
        BondViews::new(&self.bonds, &self.graph)
    }

    pub fn bond(&self, idx: BondIdx) -> BondView<'_> {
        self.bonds().get(idx)
    }

    // endregion: Read: bonds

    // region: Read: dative bonds

    pub fn dative_bond_count(&self) -> usize {
        self.dative_bonds.relation_count()
    }

    pub fn dative_bonds(&self) -> DativeBondViews<'_> {
        DativeBondViews::new(&self.dative_bonds)
    }

    pub fn dative_bond(&self, idx: DativeBondIdx) -> DativeBondView<'_> {
        self.dative_bonds().get(idx)
    }

    // endregion: Read: dative bonds

    // region: Read: aromatic systems

    pub fn aromatic_system_count(&self) -> usize {
        self.aromatic_systems.relation_count()
    }

    pub fn aromatic_systems(&self) -> AromaticSystemViews<'_> {
        AromaticSystemViews::new(&self.aromatic_systems, &self.graph)
    }

    pub fn aromatic_system(&self, idx: AromaticSystemIdx) -> AromaticSystemView<'_> {
        self.aromatic_systems().get(idx)
    }

    // endregion: Read: aromatic systems

    // region: Read: multicenter bonds

    pub fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.relation_count()
    }

    pub fn multicenter_bonds(&self) -> MulticenterBondViews<'_> {
        MulticenterBondViews::new(&self.multicenter_bonds)
    }

    pub fn multicenter_bond(&self, idx: MulticenterBondIdx) -> MulticenterBondView<'_> {
        self.multicenter_bonds().get(idx)
    }

    // endregion: Read: multicenter bonds

    // region: Read: noncovalent bonds

    pub fn noncovalent_bond_count(&self) -> usize {
        self.noncovalent_bonds.relation_count()
    }

    pub fn noncovalent_bonds(&self) -> NoncovalentBondViews<'_> {
        NoncovalentBondViews::new(&self.noncovalent_bonds)
    }

    pub fn noncovalent_bond(&self, idx: NoncovalentBondIdx) -> NoncovalentBondView<'_> {
        self.noncovalent_bonds().get(idx)
    }

    // endregion: Read: noncovalent bonds

    // region: Read: connecting relations

    pub fn connecting_bond(&self, a: AtomIdx, b: AtomIdx) -> Option<BondIdx> {
        self.graph
            .find_edge(NodeId::from(a), NodeId::from(b))
            .map(BondIdx::from)
    }

    pub fn connecting_dative_bond(
        &self,
        donor: AtomIdx,
        acceptor: AtomIdx,
    ) -> Option<DativeBondIdx> {
        self.dative_bonds_incident(donor).find(|&idx| {
            let v = self.dative_bond(idx);
            v.donor == donor && v.acceptor == acceptor
        })
    }

    pub fn connecting_noncovalent_bond(
        &self,
        a: AtomIdx,
        b: AtomIdx,
    ) -> Option<NoncovalentBondIdx> {
        self.noncovalent_bonds_incident(a).find(|&idx| {
            let v = self.noncovalent_bond(idx);
            (v.atoms[0] == a && v.atoms[1] == b) || (v.atoms[0] == b && v.atoms[1] == a)
        })
    }

    pub fn connecting_aromatic_system(
        &self,
        atoms: &HashSet<AtomIdx>,
    ) -> Option<AromaticSystemIdx> {
        let &first = atoms.iter().next()?;
        self.aromatic_systems_incident(first).find(|&idx| {
            let v = self.aromatic_system(idx);
            let parts: HashSet<AtomIdx> = v.atoms().collect();
            parts == *atoms
        })
    }

    pub fn connecting_multicenter_bond(
        &self,
        atoms: &HashSet<AtomIdx>,
    ) -> Option<MulticenterBondIdx> {
        let &first = atoms.iter().next()?;
        self.multicenter_bonds_incident(first).find(|&idx| {
            let v = self.multicenter_bond(idx);
            let parts: HashSet<AtomIdx> = v.atoms().collect();
            parts == *atoms
        })
    }

    // endregion: Read: connecting relations

    // region: Read: incidence

    pub fn dative_bonds_incident(&self, atom: AtomIdx) -> impl Iterator<Item = DativeBondIdx> + '_ {
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

    // endregion: Read: incidence

    // region: Read: induced subsets (bond/relation)

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
        self.graph
            .induced_edges(&nodes)
            .map(BondIdx::from)
            .collect()
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

    // endregion: Read: induced subsets (bond/relation)

    // region: Ring enumeration

    pub fn rings(
        &self,
        family: RingFamily,
        max_ring_size: usize,
        atom_filter: impl Fn(AtomIdx) -> bool,
    ) -> RingSet {
        rings::enumerate_rings(&self.graph, family, max_ring_size, atom_filter)
    }

    // endregion: Ring enumeration

    // region: Read: graph algorithms

    pub fn degree(&self, atom: AtomIdx) -> usize {
        self.graph.degree(NodeId::from(atom))
    }

    pub fn connected_components(&self, alg: ConnectedComponentsAlgorithm) -> Vec<Vec<AtomIdx>> {
        self.graph
            .connected_components(alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomIdx::from).collect())
            .collect()
    }

    pub fn biconnected_components(&self, alg: BiconnectedComponentsAlgorithm) -> Vec<Vec<AtomIdx>> {
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
        AtomAutomorphism(
            self.graph
                .automorphisms(|n| atom_color(AtomIdx::from(n)), alg),
        )
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

    // endregion: Read: graph algorithms

    // region: Entity mutation: atoms

    pub fn atom_mut(&mut self, idx: AtomIdx) -> AtomViewMut<'_> {
        let data = &mut Arc::make_mut(&mut self.atoms)[idx.index()];
        AtomViewMut { idx, data }
    }

    pub fn atoms_mut(&mut self) -> impl Iterator<Item = &mut AtomAst> {
        Arc::make_mut(&mut self.atoms).iter_mut()
    }

    // endregion: Entity mutation: atoms

    // region: Entity mutation: bonds

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

    // endregion: Entity mutation: bonds

    // region: Entity mutation: dative bonds

    pub fn dative_bond_mut(&mut self, idx: DativeBondIdx) -> &mut DativeBondAst {
        Arc::make_mut(&mut self.dative_bonds).data_mut(RelationId::from(idx))
    }

    pub fn dative_bonds_mut(&mut self) -> impl Iterator<Item = &mut DativeBondAst> {
        Arc::make_mut(&mut self.dative_bonds).data_iter_mut()
    }

    // endregion: Entity mutation: dative bonds

    // region: Entity mutation: aromatic systems

    pub fn aromatic_system_mut(&mut self, idx: AromaticSystemIdx) -> &mut AromaticSystemAst {
        Arc::make_mut(&mut self.aromatic_systems).data_mut(RelationId::from(idx))
    }

    pub fn aromatic_systems_mut(&mut self) -> impl Iterator<Item = &mut AromaticSystemAst> {
        Arc::make_mut(&mut self.aromatic_systems).data_iter_mut()
    }

    // endregion: Entity mutation: aromatic systems

    // region: Entity mutation: multicenter bonds

    pub fn multicenter_bond_mut(&mut self, idx: MulticenterBondIdx) -> &mut MulticenterBondAst {
        Arc::make_mut(&mut self.multicenter_bonds).data_mut(RelationId::from(idx))
    }

    pub fn multicenter_bonds_mut(&mut self) -> impl Iterator<Item = &mut MulticenterBondAst> {
        Arc::make_mut(&mut self.multicenter_bonds).data_iter_mut()
    }

    // endregion: Entity mutation: multicenter bonds

    // region: Entity mutation: noncovalent bonds

    pub fn noncovalent_bond_mut(&mut self, idx: NoncovalentBondIdx) -> &mut NoncovalentBondAst {
        Arc::make_mut(&mut self.noncovalent_bonds).data_mut(RelationId::from(idx))
    }

    pub fn noncovalent_bonds_mut(&mut self) -> impl Iterator<Item = &mut NoncovalentBondAst> {
        Arc::make_mut(&mut self.noncovalent_bonds).data_iter_mut()
    }

    // endregion: Entity mutation: noncovalent bonds

    // region: Constraints

    pub fn constraints(&self) -> &Constraints {
        &self.constraints
    }

    pub fn constraints_mut(&mut self) -> &mut Constraints {
        &mut self.constraints
    }

    /// Recursively reduce every contained `ValueAst` to canonical form
    /// via [`ValueAst::simplify`]. Walks every entity (atoms, bonds,
    /// dative/aromatic/multicenter/noncovalent), each entity's inline
    /// constraint store, and the molecule-scope `Constraints` tree —
    /// including `SubPattern` patterns recursively. Entity counts and
    /// topology are unchanged.
    pub fn simplify_values(&mut self) {
        for atom in self.atoms_mut() {
            atom.simplify_values();
        }
        for bond in self.bonds_mut() {
            bond.simplify_values();
        }
        for db in self.dative_bonds_mut() {
            db.simplify_values();
        }
        for ar in self.aromatic_systems_mut() {
            ar.simplify_values();
        }
        for mc in self.multicenter_bonds_mut() {
            mc.simplify_values();
        }
        for nc in self.noncovalent_bonds_mut() {
            nc.simplify_values();
        }
        self.constraints.simplify_each();
    }

    /// Drain every entity's inline `constraints` store into `self.constraints`
    /// as `Constraint::Atom` / `Bond` / `DativeBond` / `AromaticSystem` /
    /// `MulticenterBond` / `NoncovalentBond` entries. Each entity kind is
    /// walked unconditionally; for entities whose narrow constraint enum is
    /// currently uninhabited (aromatic system, multicenter, noncovalent),
    /// the iteration is empty until a variant is added — at which point
    /// new variants are lifted automatically with no code change here.
    /// The order of inserted entries in `self.constraints` is unspecified.
    pub fn lift_constraints(&mut self) {
        let atom_count = self.atom_count();
        let bond_count = self.bond_count();
        let dative_count = self.dative_bond_count();
        let aromatic_count = self.aromatic_system_count();
        let multicenter_count = self.multicenter_bond_count();
        let noncovalent_count = self.noncovalent_bond_count();

        let mut additions: Vec<Constraint> = Vec::new();
        for i in 0..atom_count {
            let idx = AtomIdx::from(i);
            for c in self.atom_mut(idx).data.constraints.take() {
                additions.push(Constraint::Atom(idx, c));
            }
        }
        for i in 0..bond_count {
            let idx = BondIdx::from(i);
            for c in self.bond_mut(idx).data.constraints.take() {
                additions.push(Constraint::Bond(idx, c));
            }
        }
        for i in 0..dative_count {
            let idx = DativeBondIdx::from(i);
            for c in self.dative_bond_mut(idx).constraints.take() {
                additions.push(Constraint::DativeBond(idx, c));
            }
        }
        for i in 0..aromatic_count {
            let idx = AromaticSystemIdx::from(i);
            for c in self.aromatic_system_mut(idx).constraints.take() {
                additions.push(Constraint::AromaticSystem(idx, c));
            }
        }
        for i in 0..multicenter_count {
            let idx = MulticenterBondIdx::from(i);
            for c in self.multicenter_bond_mut(idx).constraints.take() {
                additions.push(Constraint::MulticenterBond(idx, c));
            }
        }
        for i in 0..noncovalent_count {
            let idx = NoncovalentBondIdx::from(i);
            for c in self.noncovalent_bond_mut(idx).constraints.take() {
                additions.push(Constraint::NoncovalentBond(idx, c));
            }
        }
        for c in additions {
            self.constraints.push(c);
        }
    }

    /// Push every entity inline constraints from `self.constraints`
    /// into the targeted entity's inline `constraints` store via `add`
    /// (last-wins per kind), removing it from the molecule list.
    /// Combinator subtrees, `Relational`, and `Molecule` entries are left
    /// in place.
    ///
    /// The `Constraint` arm is exhaustively matched: adding a new variant
    /// or making any uninhabited entity-leaf inner enum (aromatic,
    /// multicenter, noncovalent) inhabited is a compile-time forcing
    /// function on this method.
    pub fn inline_constraints(&mut self) {
        let entries = self.constraints.take();
        let mut leftover: Vec<Constraint> = Vec::new();
        for c in entries {
            match c {
                Constraint::Atom(idx, inner) => {
                    self.atom_mut(idx).data.constraints.add(inner);
                }
                Constraint::Bond(idx, inner) => {
                    self.bond_mut(idx).data.constraints.add(inner);
                }
                Constraint::DativeBond(idx, inner) => {
                    self.dative_bond_mut(idx).constraints.add(inner);
                }
                Constraint::AromaticSystem(_, inner) => match inner {},
                Constraint::MulticenterBond(_, inner) => match inner {},
                Constraint::NoncovalentBond(_, inner) => match inner {},
                c @ (Constraint::Relational(_)
                | Constraint::Molecule(_)
                | Constraint::And(_)
                | Constraint::Or(_)
                | Constraint::Not(_)) => leftover.push(c),
            }
        }
        for c in leftover {
            self.constraints.push(c);
        }
    }

    // endregion: Constraints

    // region: Topological mutation (via builder)

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
    // endregion: Topological mutation (via builder)
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
mod tests;
