//! Read-only views over `MoleculeAst` topology and relations.
//!
//! View records bundle an index with the underlying data so consumers
//! never assemble (idx, data, participants) tuples by hand. Namespace
//! types group per-relation accessors (`count`, `ids`, `iter`, `get`,
//! and `Index`) without burying them on `MoleculeAst` itself.

use std::ops::Index;

use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm,
    CycleEnumerationAlgorithm, EdgeId, FixedRelationSet, Graph, MatchingEnumerationAlgorithm,
    MaxIndependentSetAlgorithm, MaxMatchingAlgorithm, NodeId, RelationId, ShortestCycleAlgorithm,
    SubgraphIsomorphismAlgorithm, VarRelationSet,
};

use super::automorphism::AtomAutomorphism;
use super::matching::BondMatching;

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::constraint::{
    AromaticValenceAst, AtomConstraint, AtomConstraintKind, MulticenterValenceAst,
};
use super::dative::DativeBondAst;
use super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use super::molecule::MoleculeAst;
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;
use super::value::ValueAst;

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
        self.atoms
            .iter()
            .enumerate()
            .map(move |(i, data)| AtomView {
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
    /// Iterator over incident bonds and their neighbor atoms. Equivalent to
    /// `self.ast.neighbors(self.idx)` but exposed on the view so closures
    /// that take `&AtomView` (e.g. perception electron-counting) can inspect
    /// bonds without reaching back to the molecule.
    pub fn neighbors(&self) -> impl Iterator<Item = NeighborView<'a>> {
        self.ast.neighbors(self.idx)
    }

    /// Localized valence summed from incident bond orders. `None` if any incident
    /// bond's order is not a non-negative `Lit`.
    pub fn bond_order_sum(&self) -> Option<u32> {
        let mut sum: u32 = 0;
        for n in self.ast.neighbors(self.idx) {
            match n.data.order {
                ValueAst::Lit(v) if v >= 0 => sum += v as u32,
                _ => return None,
            }
        }
        Some(sum)
    }

    /// Localized valence constraint, if asserted.
    pub fn valence_constraint(&self) -> Option<&'a ValueAst> {
        atom_constraint_value(self.data, AtomConstraintKind::Valence, |c| match c {
            AtomConstraint::Valence(v) => Some(v),
            _ => None,
        })
    }

    /// Sum of `order` over incident dative bonds where this atom is the sole
    /// donor (i.e. the dative is single-donor). Multi-donor datives contribute
    /// nothing per individual donor atom — the donated pair is collective and
    /// has no well-defined per-atom share. `None` if any contributing dative's
    /// `order` is not a non-negative `Lit`.
    pub fn donated_pairs(&self) -> Option<u32> {
        let mut sum: u32 = 0;
        for id in self.ast.dative_bonds_incident(self.idx) {
            let view = self.ast.dative_bond(id);
            let donors: Vec<_> = view.donors().collect();
            if donors.len() != 1 || donors[0] != self.idx {
                continue;
            }
            match view.data.order {
                ValueAst::Lit(v) if v >= 0 => sum += v as u32,
                _ => return None,
            }
        }
        Some(sum)
    }

    pub fn donated_pairs_constraint(&self) -> Option<&'a ValueAst> {
        atom_constraint_value(self.data, AtomConstraintKind::DonatedPairs, |c| match c {
            AtomConstraint::DonatedPairs(v) => Some(v),
            _ => None,
        })
    }

    /// Sum of `order` over incident dative bonds where this atom is the
    /// acceptor. `None` if any contributing dative's `order` is not a
    /// non-negative `Lit`.
    pub fn accepted_pairs(&self) -> Option<u32> {
        let mut sum: u32 = 0;
        for id in self.ast.dative_bonds_incident(self.idx) {
            let view = self.ast.dative_bond(id);
            if view.acceptor != self.idx {
                continue;
            }
            match view.data.order {
                ValueAst::Lit(v) if v >= 0 => sum += v as u32,
                _ => return None,
            }
        }
        Some(sum)
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
        self.ast
            .aromatic_systems_incident(self.idx)
            .next()
            .is_some()
    }

    pub fn aromatic_valence_constraint(&self) -> Option<&'a AromaticValenceAst> {
        atom_constraint_value(
            self.data,
            AtomConstraintKind::AromaticValence,
            |c| match c {
                AtomConstraint::AromaticValence(v) => Some(v),
                _ => None,
            },
        )
    }

    /// Sum of per-atom contributions across incident multicenter bonds.
    /// `None` if any contribution is not a non-negative `Lit`.
    pub fn multicenter_contribution(&self) -> Option<u32> {
        let mut sum: u32 = 0;
        for mc_id in self.ast.multicenter_bonds_incident(self.idx) {
            let view = self.ast.multicenter_bond(mc_id);
            let pos = view.atoms().position(|a| a == self.idx)?;
            match view.data.electrons.get(pos)? {
                ValueAst::Lit(v) if *v >= 0 => sum += *v as u32,
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

/// Namespace accessor for bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct BondViews<'a> {
    ast: &'a MoleculeAst,
    bonds: &'a [BondAst],
}

impl<'a> BondViews<'a> {
    pub(super) fn new(ast: &'a MoleculeAst, bonds: &'a [BondAst]) -> Self {
        Self { ast, bonds }
    }

    pub fn count(&self) -> usize {
        self.bonds.len()
    }

    pub fn ids(&self) -> impl Iterator<Item = BondIdx> {
        (0..self.bonds.len() as u32).map(BondIdx)
    }

    pub fn iter(&self) -> impl Iterator<Item = BondView<'a>> {
        let ast = self.ast;
        let bonds = self.bonds;
        let graph = ast.raw_graph();
        graph.edge_ids().map(move |id| {
            let [s, t] = graph.edge_endpoints(id);
            BondView {
                idx: BondIdx::from(id),
                atoms: [AtomIdx::from(s), AtomIdx::from(t)],
                data: &bonds[id.index()],
                ast,
            }
        })
    }

    pub fn get(&self, idx: BondIdx) -> BondView<'a> {
        let [s, t] = self.ast.raw_graph().edge_endpoints(EdgeId::from(idx));
        BondView {
            idx,
            atoms: [AtomIdx::from(s), AtomIdx::from(t)],
            data: &self.bonds[idx.index()],
            ast: self.ast,
        }
    }
}

impl<'a> Index<BondIdx> for BondViews<'a> {
    type Output = BondAst;
    fn index(&self, idx: BondIdx) -> &BondAst {
        &self.bonds[idx.index()]
    }
}

/// Borrowed view of a bond: its index, the two participating atoms, and data.
#[derive(Clone, Copy, Debug)]
pub struct BondView<'a> {
    pub idx: BondIdx,
    atoms: [AtomIdx; 2],
    pub data: &'a BondAst,
    #[allow(dead_code)]
    ast: &'a MoleculeAst,
}

impl<'a> BondView<'a> {
    /// The two atoms incident to this bond.
    pub fn atoms(&self) -> [AtomIdx; 2] {
        self.atoms
    }
}

/// Mutable borrowed view of a bond.
#[derive(Debug)]
pub struct BondViewMut<'a> {
    pub idx: BondIdx,
    atoms: [AtomIdx; 2],
    pub data: &'a mut BondAst,
}

impl<'a> BondViewMut<'a> {
    pub(super) fn new(idx: BondIdx, atoms: [AtomIdx; 2], data: &'a mut BondAst) -> Self {
        Self { idx, atoms, data }
    }

    /// The two atoms incident to this bond.
    pub fn atoms(&self) -> [AtomIdx; 2] {
        self.atoms
    }
}

/// Namespace accessor for dative-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct DativeBondViews<'a> {
    ast: &'a MoleculeAst,
    set: &'a VarRelationSet<DativeBondAst>,
}

impl<'a> DativeBondViews<'a> {
    pub(super) fn new(ast: &'a MoleculeAst, set: &'a VarRelationSet<DativeBondAst>) -> Self {
        Self { ast, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = DativeBondIdx> {
        self.set.relation_ids().map(DativeBondIdx::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = DativeBondView<'a>> {
        let ast = self.ast;
        let set = self.set;
        set.relation_ids().map(move |rid| {
            let atoms = set.participants(rid);
            let data = set.data(rid);
            let acceptor = AtomIdx::from(atoms[data.acceptor_slot as usize]);
            DativeBondView {
                idx: DativeBondIdx::from(rid),
                data,
                acceptor,
                atoms,
                ast,
            }
        })
    }

    pub fn get(&self, idx: DativeBondIdx) -> DativeBondView<'a> {
        let rid = RelationId::from(idx);
        let atoms = self.set.participants(rid);
        let data = self.set.data(rid);
        let acceptor = AtomIdx::from(atoms[data.acceptor_slot as usize]);
        DativeBondView {
            idx,
            data,
            acceptor,
            atoms,
            ast: self.ast,
        }
    }
}

impl<'a> Index<DativeBondIdx> for DativeBondViews<'a> {
    type Output = DativeBondAst;
    fn index(&self, idx: DativeBondIdx) -> &DativeBondAst {
        self.set.data(RelationId::from(idx))
    }
}

/// Borrowed view of a dative bond: index, the designated acceptor atom,
/// and underlying `DativeBondAst`. Donor atoms and the full participant
/// set are reachable through `donors()` and `atoms()`.
#[derive(Clone, Copy, Debug)]
pub struct DativeBondView<'a> {
    pub idx: DativeBondIdx,
    pub acceptor: AtomIdx,
    atoms: &'a [NodeId],
    pub data: &'a DativeBondAst,
    #[allow(dead_code)]
    ast: &'a MoleculeAst,
}

impl<'a> DativeBondView<'a> {
    /// All atoms in this dative bond (donors + acceptor), sorted by `AtomIdx`.
    pub fn atoms(&self) -> impl Iterator<Item = AtomIdx> + '_ {
        self.atoms.iter().map(|&n| AtomIdx::from(n))
    }

    /// Donor atoms (participants minus the acceptor slot).
    pub fn donors(&self) -> impl Iterator<Item = AtomIdx> + '_ {
        let acceptor_slot = self.data.acceptor_slot as usize;
        self.atoms
            .iter()
            .enumerate()
            .filter(move |(i, _)| *i != acceptor_slot)
            .map(|(_, &n)| AtomIdx::from(n))
    }
}

/// Namespace accessor for aromatic-system views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct AromaticSystemViews<'a> {
    ast: &'a MoleculeAst,
    set: &'a VarRelationSet<AromaticSystemAst>,
}

impl<'a> AromaticSystemViews<'a> {
    pub(super) fn new(ast: &'a MoleculeAst, set: &'a VarRelationSet<AromaticSystemAst>) -> Self {
        Self { ast, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = AromaticSystemIdx> {
        self.set.relation_ids().map(AromaticSystemIdx::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = AromaticSystemView<'a>> {
        let ast = self.ast;
        let set = self.set;
        set.relation_ids().map(move |rid| AromaticSystemView {
            idx: AromaticSystemIdx::from(rid),
            data: set.data(rid),
            atoms: set.participants(rid),
            ast,
        })
    }

    pub fn get(&self, idx: AromaticSystemIdx) -> AromaticSystemView<'a> {
        let rid = RelationId::from(idx);
        AromaticSystemView {
            idx,
            data: self.set.data(rid),
            atoms: self.set.participants(rid),
            ast: self.ast,
        }
    }
}

impl<'a> Index<AromaticSystemIdx> for AromaticSystemViews<'a> {
    type Output = AromaticSystemAst;
    fn index(&self, idx: AromaticSystemIdx) -> &AromaticSystemAst {
        self.set.data(RelationId::from(idx))
    }
}

/// Borrowed view of an aromatic system: its index, the `AromaticSystemAst`,
/// and accessors for member atoms and induced ring bonds via `atoms()` and
/// `bonds()`.
#[derive(Clone, Copy, Debug)]
pub struct AromaticSystemView<'a> {
    pub idx: AromaticSystemIdx,
    atoms: &'a [NodeId],
    pub data: &'a AromaticSystemAst,
    ast: &'a MoleculeAst,
}

impl<'a> AromaticSystemView<'a> {
    pub fn atoms(&self) -> impl Iterator<Item = AtomIdx> + '_ {
        self.atoms.iter().map(|&n| AtomIdx::from(n))
    }

    pub fn bonds(&self) -> impl Iterator<Item = BondIdx> + '_ {
        self.ast
            .raw_graph()
            .induced_edges(self.atoms)
            .map(BondIdx::from)
    }
}

/// Namespace accessor for multicenter-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct MulticenterBondViews<'a> {
    ast: &'a MoleculeAst,
    set: &'a VarRelationSet<MulticenterBondAst>,
}

impl<'a> MulticenterBondViews<'a> {
    pub(super) fn new(ast: &'a MoleculeAst, set: &'a VarRelationSet<MulticenterBondAst>) -> Self {
        Self { ast, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = MulticenterBondIdx> {
        self.set.relation_ids().map(MulticenterBondIdx::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = MulticenterBondView<'a>> {
        let ast = self.ast;
        let set = self.set;
        set.relation_ids().map(move |rid| MulticenterBondView {
            idx: MulticenterBondIdx::from(rid),
            data: set.data(rid),
            atoms: set.participants(rid),
            ast,
        })
    }

    pub fn get(&self, idx: MulticenterBondIdx) -> MulticenterBondView<'a> {
        let rid = RelationId::from(idx);
        MulticenterBondView {
            idx,
            data: self.set.data(rid),
            atoms: self.set.participants(rid),
            ast: self.ast,
        }
    }
}

impl<'a> Index<MulticenterBondIdx> for MulticenterBondViews<'a> {
    type Output = MulticenterBondAst;
    fn index(&self, idx: MulticenterBondIdx) -> &MulticenterBondAst {
        self.set.data(RelationId::from(idx))
    }
}

/// Borrowed view of a multicenter bond: its index, member atoms via
/// `atoms()`, and underlying `MulticenterBondAst`.
#[derive(Clone, Copy, Debug)]
pub struct MulticenterBondView<'a> {
    pub idx: MulticenterBondIdx,
    atoms: &'a [NodeId],
    pub data: &'a MulticenterBondAst,
    #[allow(dead_code)]
    ast: &'a MoleculeAst,
}

impl<'a> MulticenterBondView<'a> {
    pub fn atoms(&self) -> impl Iterator<Item = AtomIdx> + '_ {
        self.atoms.iter().map(|&n| AtomIdx::from(n))
    }
}

/// Namespace accessor for noncovalent-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct NoncovalentBondViews<'a> {
    ast: &'a MoleculeAst,
    set: &'a FixedRelationSet<NoncovalentBondAst, 2>,
}

impl<'a> NoncovalentBondViews<'a> {
    pub(super) fn new(
        ast: &'a MoleculeAst,
        set: &'a FixedRelationSet<NoncovalentBondAst, 2>,
    ) -> Self {
        Self { ast, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = NoncovalentBondIdx> {
        self.set.relation_ids().map(NoncovalentBondIdx::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = NoncovalentBondView<'a>> {
        let ast = self.ast;
        let set = self.set;
        set.relation_ids().map(move |rid| NoncovalentBondView {
            idx: NoncovalentBondIdx::from(rid),
            data: set.data(rid),
            atoms: {
                let parts = set.participants(rid);
                [AtomIdx::from(parts[0]), AtomIdx::from(parts[1])]
            },
            ast,
        })
    }

    pub fn get(&self, idx: NoncovalentBondIdx) -> NoncovalentBondView<'a> {
        let rid = RelationId::from(idx);
        let parts = self.set.participants(rid);
        NoncovalentBondView {
            idx,
            data: self.set.data(rid),
            atoms: [AtomIdx::from(parts[0]), AtomIdx::from(parts[1])],
            ast: self.ast,
        }
    }
}

impl<'a> Index<NoncovalentBondIdx> for NoncovalentBondViews<'a> {
    type Output = NoncovalentBondAst;
    fn index(&self, idx: NoncovalentBondIdx) -> &NoncovalentBondAst {
        self.set.data(RelationId::from(idx))
    }
}

/// Borrowed view of a noncovalent bond: the two participating atoms plus data.
#[derive(Clone, Copy, Debug)]
pub struct NoncovalentBondView<'a> {
    pub idx: NoncovalentBondIdx,
    atoms: [AtomIdx; 2],
    pub data: &'a NoncovalentBondAst,
    #[allow(dead_code)]
    ast: &'a MoleculeAst,
}

impl<'a> NoncovalentBondView<'a> {
    /// The two atoms in this noncovalent interaction.
    pub fn atoms(&self) -> [AtomIdx; 2] {
        self.atoms
    }
}

/// Neighbor-side view of a bond: the atom on the other end (`atom`), the
/// bond index, the bond data, and the parent `MoleculeAst` for navigation
/// to the neighbor's full atom view. Yielded by `MoleculeAst::neighbors`.
#[derive(Clone, Copy, Debug)]
pub struct NeighborView<'a> {
    pub bond: BondIdx,
    pub atom: AtomIdx,
    pub data: &'a BondAst,
    #[allow(dead_code)]
    ast: &'a MoleculeAst,
}

impl<'a> NeighborView<'a> {
    pub(super) fn new(
        bond: BondIdx,
        atom: AtomIdx,
        data: &'a BondAst,
        ast: &'a MoleculeAst,
    ) -> Self {
        Self {
            bond,
            atom,
            data,
            ast,
        }
    }
}

/// AtomIdx/BondIdx-typed adapter over the underlying `Graph`. Holds the
/// pure-graph algorithms (connectivity, cycles, matchings, isomorphisms)
/// without exposing graph-core's `NodeId` / `EdgeId` types in the public
/// API. Construct via `MoleculeAst::graph()`.
#[derive(Clone, Copy)]
pub struct GraphView<'a> {
    graph: &'a Graph,
}

impl<'a> GraphView<'a> {
    pub(super) fn new(graph: &'a Graph) -> Self {
        Self { graph }
    }

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
        AtomAutomorphism(
            self.graph
                .automorphisms(|n| atom_color(AtomIdx::from(n)), alg),
        )
    }

    pub fn subgraph_isomorphisms(
        &self,
        query: &GraphView<'_>,
        atom_match: &mut impl FnMut(AtomIdx, AtomIdx) -> bool,
        bond_match: &mut impl FnMut(BondIdx, BondIdx) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<AtomIdx>> {
        self.graph
            .subgraph_isomorphisms(
                query.graph,
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
        query: &GraphView<'_>,
        anchor: (AtomIdx, AtomIdx),
        atom_match: &mut impl FnMut(AtomIdx, AtomIdx) -> bool,
        bond_match: &mut impl FnMut(BondIdx, BondIdx) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<AtomIdx>> {
        self.graph
            .subgraph_isomorphisms_at(
                query.graph,
                (NodeId::from(anchor.0), NodeId::from(anchor.1)),
                &mut |tn, qn| atom_match(AtomIdx::from(tn), AtomIdx::from(qn)),
                &mut |te, qe| bond_match(BondIdx::from(te), BondIdx::from(qe)),
                alg,
            )
            .into_iter()
            .map(|m| m.into_iter().map(AtomIdx::from).collect())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::{AtomConstraint, Constraints};
    use crate::ast::dative::DativeBondAst;
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use crate::ast::value::ValueAst;
    use crate::mol;

    /// 4-atom molecule with one of every relation kind:
    /// atoms C C N O; bonds 0-1 single, 1-2 double, 2-3 single;
    /// dative donor=2 → acceptor=3; aromatic system [0,1,2];
    /// multicenter bond [0,1,2]; noncovalent H-bond 0-3.
    #[fixture]
    fn molecule() -> MoleculeAst {
        MoleculeAst::from_parts(
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
            vec![(vec![AtomIdx(2)], AtomIdx(3), DativeBondAst::from_order(1))],
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

    // --- AtomViews ---

    #[rstest]
    fn test_atom_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.atoms().count(), 4);
    }

    #[rstest]
    fn test_atom_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.atoms().ids().collect::<Vec<_>>(),
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)],
        );
    }

    #[rstest]
    fn test_atom_views_iter(molecule: MoleculeAst) {
        let views = molecule.atoms();
        let collected: Vec<(AtomIdx, AtomAst)> =
            views.iter().map(|v| (v.idx, v.data.clone())).collect();
        assert_eq!(
            collected,
            vec![
                (AtomIdx(0), AtomAst::from_element(Element::C)),
                (AtomIdx(1), AtomAst::from_element(Element::C)),
                (AtomIdx(2), AtomAst::from_element(Element::N)),
                (AtomIdx(3), AtomAst::from_element(Element::O)),
            ],
        );
    }

    #[rstest]
    fn test_atom_views_get(molecule: MoleculeAst) {
        let view = molecule.atoms().get(AtomIdx(2));
        assert_eq!(view.idx, AtomIdx(2));
        assert_eq!(*view.data, AtomAst::from_element(Element::N));
    }

    #[rstest]
    fn test_atom_views_index(molecule: MoleculeAst) {
        let atom: &AtomAst = &molecule.atoms()[AtomIdx(2)];
        assert_eq!(*atom, AtomAst::from_element(Element::N));
    }

    // --- AtomView ---

    #[rstest]
    fn test_atom_view_neighbors(molecule: MoleculeAst) {
        let view = molecule.atom(AtomIdx(1));
        let collected: Vec<(BondIdx, AtomIdx, BondAst)> = view
            .neighbors()
            .map(|n| (n.bond, n.atom, n.data.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondIdx(0), AtomIdx(0), BondAst::from_order(1)),
                (BondIdx(1), AtomIdx(2), BondAst::from_order(2)),
            ],
        );
    }

    #[rstest]
    #[case::no_incident(AtomIdx(3), Some(0))]
    #[case::single(AtomIdx(0), Some(1))]
    #[case::three_around_center(AtomIdx(1), Some(3))]
    #[case::double(AtomIdx(2), Some(2))]
    fn test_atom_view_bond_order_sum(#[case] center: AtomIdx, #[case] expected: Option<u32>) {
        let ast = mol!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#);
        assert_eq!(ast.atom(center).bond_order_sum(), expected);
    }

    #[rstest]
    fn test_atom_view_bond_order_sum_undetermined() {
        let ast = mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "*"]]}"#);
        assert_eq!(ast.atom(AtomIdx(0)).bond_order_sum(), None);
    }

    #[rstest]
    fn test_atom_view_valence_constraint() {
        let mut atom = AtomAst::from_element(Element::C);
        atom.constraints.add(AtomConstraint::valence(4));
        let ast = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            ast.atom(AtomIdx(0)).valence_constraint(),
            Some(&ValueAst::Lit(4)),
        );
    }

    #[rstest]
    fn test_atom_view_valence_constraint_absent(molecule: MoleculeAst) {
        assert!(molecule.atom(AtomIdx(0)).valence_constraint().is_none());
    }

    #[rstest]
    #[case::donor(AtomIdx(0), Some(1))]
    #[case::acceptor(AtomIdx(1), Some(0))]
    fn test_atom_view_donated_pairs(#[case] atom: AtomIdx, #[case] expected: Option<u32>) {
        let ast = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![(vec![AtomIdx(0)], AtomIdx(1), DativeBondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(ast.atom(atom).donated_pairs(), expected);
    }

    #[rstest]
    fn test_atom_view_donated_pairs_constraint() {
        let mut atom = AtomAst::from_element(Element::N);
        atom.constraints.add(AtomConstraint::donated_pairs(1));
        let ast = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            ast.atom(AtomIdx(0)).donated_pairs_constraint(),
            Some(&ValueAst::Lit(1)),
        );
    }

    #[rstest]
    #[case::donor(AtomIdx(0), Some(0))]
    #[case::acceptor(AtomIdx(1), Some(1))]
    fn test_atom_view_accepted_pairs(#[case] atom: AtomIdx, #[case] expected: Option<u32>) {
        let ast = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![(vec![AtomIdx(0)], AtomIdx(1), DativeBondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(ast.atom(atom).accepted_pairs(), expected);
    }

    #[rstest]
    fn test_atom_view_accepted_pairs_constraint() {
        let mut atom = AtomAst::from_element(Element::C);
        atom.constraints.add(AtomConstraint::accepted_pairs(2));
        let ast = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            ast.atom(AtomIdx(0)).accepted_pairs_constraint(),
            Some(&ValueAst::Lit(2)),
        );
    }

    #[rstest]
    #[case::lit(ValueAst::Lit(2), Some(2))]
    #[case::undetermined(ValueAst::Undetermined, None)]
    fn test_atom_view_aromatic_contribution(
        #[case] entry: ValueAst,
        #[case] expected: Option<u32>,
    ) {
        let ast = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![(
                vec![AtomIdx(0), AtomIdx(1)],
                AromaticSystemAst::new(vec![entry, ValueAst::Lit(1)]),
            )],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(ast.atom(AtomIdx(0)).aromatic_contribution(), expected);
    }

    #[rstest]
    fn test_atom_view_aromatic_contribution_not_in_system() {
        let ast = mol!(r#"{:atoms ["C"] :bonds []}"#);
        assert_eq!(ast.atom(AtomIdx(0)).aromatic_contribution(), Some(0));
    }

    #[rstest]
    #[case::in_system(AtomIdx(0), true)]
    #[case::not_in_system(AtomIdx(3), false)]
    fn test_atom_view_is_in_aromatic_system(
        molecule: MoleculeAst,
        #[case] atom: AtomIdx,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.atom(atom).is_in_aromatic_system(), expected);
    }

    #[rstest]
    fn test_atom_view_aromatic_valence_constraint() {
        let mut atom = AtomAst::from_element(Element::C);
        atom.constraints.add(AtomConstraint::aromatic_valence(
            AromaticValenceAst::Aromatic(ValueAst::Lit(1)),
        ));
        let ast = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            ast.atom(AtomIdx(0)).aromatic_valence_constraint(),
            Some(&AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
        );
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
        let multicenter: Vec<_> = bonds
            .into_iter()
            .map(|(parts, electrons)| (parts, MulticenterBondAst::new(electrons)))
            .collect();
        let ast = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
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
    fn test_atom_view_multicenter_valence_constraint() {
        let mut atom = AtomAst::from_element(Element::C);
        atom.constraints.add(AtomConstraint::multicenter_valence(
            MulticenterValenceAst::Multicenter(ValueAst::Lit(2)),
        ));
        let ast = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            ast.atom(AtomIdx(0)).multicenter_valence_constraint(),
            Some(&MulticenterValenceAst::Multicenter(ValueAst::Lit(2))),
        );
    }

    // --- BondViews ---

    #[rstest]
    fn test_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.bonds().count(), 3);
    }

    #[rstest]
    fn test_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.bonds().ids().collect::<Vec<_>>(),
            vec![BondIdx(0), BondIdx(1), BondIdx(2)],
        );
    }

    #[rstest]
    fn test_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(BondIdx, [AtomIdx; 2], BondAst)> = molecule
            .bonds()
            .iter()
            .map(|v| (v.idx, v.atoms(), v.data.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondIdx(0), [AtomIdx(0), AtomIdx(1)], BondAst::from_order(1)),
                (BondIdx(1), [AtomIdx(1), AtomIdx(2)], BondAst::from_order(2)),
                (BondIdx(2), [AtomIdx(2), AtomIdx(3)], BondAst::from_order(1)),
            ],
        );
    }

    #[rstest]
    fn test_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.bonds().get(BondIdx(1));
        assert_eq!(view.idx, BondIdx(1));
        assert_eq!(view.atoms(), [AtomIdx(1), AtomIdx(2)]);
        assert_eq!(*view.data, BondAst::from_order(2));
    }

    #[rstest]
    fn test_bond_views_index(molecule: MoleculeAst) {
        let bond: &BondAst = &molecule.bonds()[BondIdx(1)];
        assert_eq!(*bond, BondAst::from_order(2));
    }

    // --- BondView ---

    #[rstest]
    fn test_bond_view_atoms(molecule: MoleculeAst) {
        assert_eq!(molecule.bond(BondIdx(1)).atoms(), [AtomIdx(1), AtomIdx(2)]);
    }

    // --- DativeBondViews ---

    #[rstest]
    fn test_dative_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.dative_bonds().count(), 1);
    }

    #[rstest]
    fn test_dative_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.dative_bonds().ids().collect::<Vec<_>>(),
            vec![DativeBondIdx(0)],
        );
    }

    #[rstest]
    fn test_dative_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(DativeBondIdx, AtomIdx, DativeBondAst)> = molecule
            .dative_bonds()
            .iter()
            .map(|v| (v.idx, v.acceptor, v.data.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                DativeBondIdx(0),
                AtomIdx(3),
                DativeBondAst::from_order(1).with_acceptor_slot(1),
            )],
        );
    }

    #[rstest]
    fn test_dative_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.dative_bonds().get(DativeBondIdx(0));
        assert_eq!(view.idx, DativeBondIdx(0));
        assert_eq!(view.acceptor, AtomIdx(3));
    }

    #[rstest]
    fn test_dative_bond_views_index(molecule: MoleculeAst) {
        let dative: &DativeBondAst = &molecule.dative_bonds()[DativeBondIdx(0)];
        assert_eq!(dative.order, ValueAst::Lit(1));
    }

    // --- DativeBondView ---

    #[rstest]
    fn test_dative_bond_view_atoms(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .dative_bond(DativeBondIdx(0))
                .atoms()
                .collect::<Vec<_>>(),
            vec![AtomIdx(2), AtomIdx(3)],
        );
    }

    #[rstest]
    fn test_dative_bond_view_donors(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .dative_bond(DativeBondIdx(0))
                .donors()
                .collect::<Vec<_>>(),
            vec![AtomIdx(2)],
        );
    }

    #[rstest]
    fn test_dative_bond_view_acceptor(molecule: MoleculeAst) {
        assert_eq!(molecule.dative_bond(DativeBondIdx(0)).acceptor, AtomIdx(3));
    }

    // --- AromaticSystemViews ---

    #[rstest]
    fn test_aromatic_system_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.aromatic_systems().count(), 1);
    }

    #[rstest]
    fn test_aromatic_system_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.aromatic_systems().ids().collect::<Vec<_>>(),
            vec![AromaticSystemIdx(0)],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(AromaticSystemIdx, Vec<AtomIdx>)> = molecule
            .aromatic_systems()
            .iter()
            .map(|v| (v.idx, v.atoms().collect()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                AromaticSystemIdx(0),
                vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]
            )],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_get(molecule: MoleculeAst) {
        let view = molecule.aromatic_systems().get(AromaticSystemIdx(0));
        assert_eq!(view.idx, AromaticSystemIdx(0));
        assert_eq!(
            view.atoms().collect::<Vec<_>>(),
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_index(molecule: MoleculeAst) {
        let _: &AromaticSystemAst = &molecule.aromatic_systems()[AromaticSystemIdx(0)];
    }

    // --- AromaticSystemView ---

    #[rstest]
    fn test_aromatic_system_view_atoms(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemIdx(0))
                .atoms()
                .collect::<Vec<_>>(),
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
        );
    }

    #[rstest]
    fn test_aromatic_system_view_bonds(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemIdx(0))
                .bonds()
                .collect::<Vec<_>>(),
            vec![BondIdx(0), BondIdx(1)],
        );
    }

    // --- MulticenterBondViews ---

    #[rstest]
    fn test_multicenter_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.multicenter_bonds().count(), 1);
    }

    #[rstest]
    fn test_multicenter_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.multicenter_bonds().ids().collect::<Vec<_>>(),
            vec![MulticenterBondIdx(0)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(MulticenterBondIdx, Vec<AtomIdx>)> = molecule
            .multicenter_bonds()
            .iter()
            .map(|v| (v.idx, v.atoms().collect()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                MulticenterBondIdx(0),
                vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
            )],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.multicenter_bonds().get(MulticenterBondIdx(0));
        assert_eq!(view.idx, MulticenterBondIdx(0));
        assert_eq!(
            view.atoms().collect::<Vec<_>>(),
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_index(molecule: MoleculeAst) {
        let _: &MulticenterBondAst = &molecule.multicenter_bonds()[MulticenterBondIdx(0)];
    }

    // --- MulticenterBondView ---

    #[rstest]
    fn test_multicenter_bond_view_atoms(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .multicenter_bond(MulticenterBondIdx(0))
                .atoms()
                .collect::<Vec<_>>(),
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
        );
    }

    // --- NoncovalentBondViews ---

    #[rstest]
    fn test_noncovalent_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.noncovalent_bonds().count(), 1);
    }

    #[rstest]
    fn test_noncovalent_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.noncovalent_bonds().ids().collect::<Vec<_>>(),
            vec![NoncovalentBondIdx(0)],
        );
    }

    #[rstest]
    fn test_noncovalent_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(NoncovalentBondIdx, [AtomIdx; 2], NoncovalentBondAst)> = molecule
            .noncovalent_bonds()
            .iter()
            .map(|v| (v.idx, v.atoms(), v.data.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                NoncovalentBondIdx(0),
                [AtomIdx(0), AtomIdx(3)],
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
        );
    }

    #[rstest]
    fn test_noncovalent_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.noncovalent_bonds().get(NoncovalentBondIdx(0));
        assert_eq!(view.idx, NoncovalentBondIdx(0));
        assert_eq!(view.atoms(), [AtomIdx(0), AtomIdx(3)]);
    }

    #[rstest]
    fn test_noncovalent_bond_views_index(molecule: MoleculeAst) {
        let _: &NoncovalentBondAst = &molecule.noncovalent_bonds()[NoncovalentBondIdx(0)];
    }

    // --- NoncovalentBondView ---

    #[rstest]
    fn test_noncovalent_bond_view_atoms(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .noncovalent_bond(NoncovalentBondIdx(0))
                .atoms(),
            [AtomIdx(0), AtomIdx(3)],
        );
    }

    // --- NeighborView ---

    #[rstest]
    fn test_neighbor_view_fields(molecule: MoleculeAst) {
        let collected: Vec<(BondIdx, AtomIdx, BondAst)> = molecule
            .neighbors(AtomIdx(2))
            .map(|n| (n.bond, n.atom, n.data.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondIdx(1), AtomIdx(1), BondAst::from_order(2)),
                (BondIdx(2), AtomIdx(3), BondAst::from_order(1)),
            ],
        );
    }
}
