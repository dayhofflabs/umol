//! Read-only views over `MoleculeAst` topology and relations.
//!
//! View records bundle an index with the underlying data so consumers
//! never assemble (id, data, participants) tuples by hand. Namespace
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
use super::atom::{AtomAst, ElementAst, ImplicitHydrogensAst, IsotopeAst};
use super::bond::BondAst;
use super::constraint::{
    AromaticSystemConstraints, AtomConstraints, BondConstraints, DativeBondConstraints,
    MulticenterBondConstraints, NoncovalentBondConstraints,
};
use super::dative::DativeBondAst;
use super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::molecule::MoleculeAst;
use super::multicenter::MulticenterBondAst;
use super::noncovalent::{NoncovalentBondAst, NoncovalentBondKindAst};
use super::spin::SpinStateAst;
use super::value::ValueAst;

/// Namespace accessor for atom views on a `MoleculeAst`. Provides `count`,
/// `ids`, `iter`, `get`, and `Index` without burying them on `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct AtomViews<'a> {
    molecule: &'a MoleculeAst,
    atoms: &'a [AtomAst],
}

impl<'a> AtomViews<'a> {
    pub(super) fn new(molecule: &'a MoleculeAst, atoms: &'a [AtomAst]) -> Self {
        Self { molecule, atoms }
    }

    pub fn count(&self) -> usize {
        self.atoms.len()
    }

    pub fn ids(&self) -> impl Iterator<Item = AtomId> {
        (0..self.atoms.len() as u32).map(AtomId)
    }

    pub fn iter(&self) -> impl Iterator<Item = AtomView<'a>> {
        let molecule = self.molecule;
        self.atoms
            .iter()
            .enumerate()
            .map(move |(i, ast)| AtomView {
                id: AtomId(i as u32),
                ast,
                molecule,
            })
    }

    pub fn get(&self, id: AtomId) -> AtomView<'a> {
        AtomView {
            id,
            ast: &self.atoms[id.index()],
            molecule: self.molecule,
        }
    }
}

impl<'a> Index<AtomId> for AtomViews<'a> {
    type Output = AtomAst;
    fn index(&self, id: AtomId) -> &AtomAst {
        &self.atoms[id.index()]
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
    pub id: AtomId,
    pub ast: &'a AtomAst,
    molecule: &'a MoleculeAst,
}

impl<'a> AtomView<'a> {
    #[inline]
    pub fn element(&self) -> &'a ElementAst {
        &self.ast.element
    }

    #[inline]
    pub fn isotope_mass(&self) -> &'a IsotopeAst {
        &self.ast.isotope_mass
    }

    #[inline]
    pub fn charge(&self) -> &'a ValueAst {
        &self.ast.charge
    }

    #[inline]
    pub fn implicit_hydrogens(&self) -> &'a ImplicitHydrogensAst {
        &self.ast.implicit_hydrogens
    }

    #[inline]
    pub fn lone_pairs(&self) -> &'a ValueAst {
        &self.ast.lone_pairs
    }

    #[inline]
    pub fn spin(&self) -> &'a SpinStateAst {
        &self.ast.spin
    }

    #[inline]
    pub fn constraints(&self) -> &'a AtomConstraints {
        &self.ast.constraints
    }

    /// Iterator over incident bonds and their neighbor atoms. Equivalent to
    /// `self.molecule.neighbors(self.id)` but exposed on the view so closures
    /// that take `&AtomView` (e.g. perception electron-counting) can inspect
    /// bonds without reaching back to the molecule.
    pub fn neighbors(&self) -> impl Iterator<Item = NeighborView<'a>> {
        self.molecule.neighbors(self.id)
    }

    /// Localized valence summed from incident bond orders. `None` if any incident
    /// bond's order is not a non-negative `Lit`.
    pub fn valence(&self) -> Option<u32> {
        let mut sum: u32 = 0;
        for n in self.molecule.neighbors(self.id) {
            match n.ast.order {
                ValueAst::Lit(v) if v >= 0 => sum += v as u32,
                _ => return None,
            }
        }
        Some(sum)
    }

    /// Sum of `order` over incident dative bonds where this atom is the sole
    /// donor (i.e. the dative is single-donor). Multi-donor datives contribute
    /// nothing per individual donor atom — the donated pair is collective and
    /// has no well-defined per-atom share. `None` if any contributing dative's
    /// `order` is not a non-negative `Lit`.
    pub fn donated_pairs(&self) -> Option<u32> {
        let mut sum: u32 = 0;
        for id in self.molecule.dative_bonds_incident(self.id) {
            let view = self.molecule.dative_bond(id);
            let donors: Vec<_> = view.donors().collect();
            if donors.len() != 1 || donors[0] != self.id {
                continue;
            }
            match view.ast.order {
                ValueAst::Lit(v) if v >= 0 => sum += v as u32,
                _ => return None,
            }
        }
        Some(sum)
    }

    /// Sum of `order` over incident dative bonds where this atom is the
    /// acceptor. `None` if any contributing dative's `order` is not a
    /// non-negative `Lit`.
    pub fn accepted_pairs(&self) -> Option<u32> {
        let mut sum: u32 = 0;
        for id in self.molecule.dative_bonds_incident(self.id) {
            let view = self.molecule.dative_bond(id);
            if view.acceptor != self.id {
                continue;
            }
            match view.ast.order {
                ValueAst::Lit(v) if v >= 0 => sum += v as u32,
                _ => return None,
            }
        }
        Some(sum)
    }

    /// π contribution from the aromatic system this atom belongs to.
    /// `Some(0)` if the atom is not in any aromatic system. `None` if the
    /// recorded contribution is not a non-negative `Lit`.
    ///
    /// An atom belongs to at most one aromatic system; the first incident
    /// system is consulted.
    pub fn aromatic_valence(&self) -> Option<u32> {
        let Some(sys_id) = self.molecule.aromatic_systems_incident(self.id).next() else {
            return Some(0);
        };
        let view = self.molecule.aromatic_system(sys_id);
        let pos = view.atoms().position(|a| a == self.id)?;
        match view.ast.electrons.get(pos)? {
            ValueAst::Lit(v) if *v >= 0 => Some(*v as u32),
            _ => None,
        }
    }

    pub fn is_in_aromatic_system(&self) -> bool {
        self.molecule
            .aromatic_systems_incident(self.id)
            .next()
            .is_some()
    }

    /// Sum of per-atom contributions across incident multicenter bonds.
    /// `None` if any contribution is not a non-negative `Lit`.
    pub fn multicenter_valence(&self) -> Option<u32> {
        let mut sum: u32 = 0;
        for mc_id in self.molecule.multicenter_bonds_incident(self.id) {
            let view = self.molecule.multicenter_bond(mc_id);
            let pos = view.atoms().position(|a| a == self.id)?;
            match view.ast.electrons.get(pos)? {
                ValueAst::Lit(v) if *v >= 0 => sum += *v as u32,
                _ => return None,
            }
        }
        Some(sum)
    }

}

/// Mutable borrowed view of an atom.
#[derive(Debug)]
pub struct AtomViewMut<'a> {
    pub id: AtomId,
    pub ast: &'a mut AtomAst,
}

/// Namespace accessor for bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct BondViews<'a> {
    molecule: &'a MoleculeAst,
    bonds: &'a [BondAst],
}

impl<'a> BondViews<'a> {
    pub(super) fn new(molecule: &'a MoleculeAst, bonds: &'a [BondAst]) -> Self {
        Self { molecule, bonds }
    }

    pub fn count(&self) -> usize {
        self.bonds.len()
    }

    pub fn ids(&self) -> impl Iterator<Item = BondId> {
        (0..self.bonds.len() as u32).map(BondId)
    }

    pub fn iter(&self) -> impl Iterator<Item = BondView<'a>> {
        let molecule = self.molecule;
        let bonds = self.bonds;
        let graph = molecule.raw_graph();
        graph.edge_ids().map(move |id| {
            let [s, t] = graph.edge_endpoints(id);
            BondView {
                id: BondId::from(id),
                atoms: [AtomId::from(s), AtomId::from(t)],
                ast: &bonds[id.index()],
                molecule,
            }
        })
    }

    pub fn get(&self, id: BondId) -> BondView<'a> {
        let [s, t] = self.molecule.raw_graph().edge_endpoints(EdgeId::from(id));
        BondView {
            id,
            atoms: [AtomId::from(s), AtomId::from(t)],
            ast: &self.bonds[id.index()],
            molecule: self.molecule,
        }
    }
}

impl<'a> Index<BondId> for BondViews<'a> {
    type Output = BondAst;
    fn index(&self, id: BondId) -> &BondAst {
        &self.bonds[id.index()]
    }
}

/// Borrowed view of a bond: its index, the two participating atoms, and data.
#[derive(Clone, Copy, Debug)]
pub struct BondView<'a> {
    pub id: BondId,
    atoms: [AtomId; 2],
    pub ast: &'a BondAst,
    #[allow(dead_code)]
    molecule: &'a MoleculeAst,
}

impl<'a> BondView<'a> {
    #[inline]
    pub fn order(&self) -> &'a ValueAst {
        &self.ast.order
    }

    #[inline]
    pub fn charge(&self) -> &'a ValueAst {
        &self.ast.charge
    }

    #[inline]
    pub fn spin(&self) -> &'a SpinStateAst {
        &self.ast.spin
    }

    #[inline]
    pub fn constraints(&self) -> &'a BondConstraints {
        &self.ast.constraints
    }

    /// The two atoms incident to this bond.
    pub fn atoms(&self) -> [AtomId; 2] {
        self.atoms
    }
}

/// Mutable borrowed view of a bond.
#[derive(Debug)]
pub struct BondViewMut<'a> {
    pub id: BondId,
    atoms: [AtomId; 2],
    pub ast: &'a mut BondAst,
}

impl<'a> BondViewMut<'a> {
    pub(super) fn new(id: BondId, atoms: [AtomId; 2], ast: &'a mut BondAst) -> Self {
        Self { id, atoms, ast }
    }

    /// The two atoms incident to this bond.
    pub fn atoms(&self) -> [AtomId; 2] {
        self.atoms
    }
}

/// Namespace accessor for dative-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct DativeBondViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a VarRelationSet<DativeBondAst>,
}

impl<'a> DativeBondViews<'a> {
    pub(super) fn new(molecule: &'a MoleculeAst, set: &'a VarRelationSet<DativeBondAst>) -> Self {
        Self { molecule, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = DativeBondId> {
        self.set.relation_ids().map(DativeBondId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = DativeBondView<'a>> {
        let molecule = self.molecule;
        let set = self.set;
        set.relation_ids().map(move |rid| {
            let atoms = set.participants(rid);
            let ast = set.data(rid);
            let acceptor = AtomId::from(atoms[ast.acceptor_slot as usize]);
            DativeBondView {
                id: DativeBondId::from(rid),
                ast,
                acceptor,
                atoms,
                molecule,
            }
        })
    }

    pub fn get(&self, id: DativeBondId) -> DativeBondView<'a> {
        let rid = RelationId::from(id);
        let atoms = self.set.participants(rid);
        let ast = self.set.data(rid);
        let acceptor = AtomId::from(atoms[ast.acceptor_slot as usize]);
        DativeBondView {
            id,
            ast,
            acceptor,
            atoms,
            molecule: self.molecule,
        }
    }
}

impl<'a> Index<DativeBondId> for DativeBondViews<'a> {
    type Output = DativeBondAst;
    fn index(&self, id: DativeBondId) -> &DativeBondAst {
        self.set.data(RelationId::from(id))
    }
}

/// Borrowed view of a dative bond: index, the designated acceptor atom,
/// and underlying `DativeBondAst`. Donor atoms and the full participant
/// set are reachable through `donors()` and `atoms()`.
#[derive(Clone, Copy, Debug)]
pub struct DativeBondView<'a> {
    pub id: DativeBondId,
    pub acceptor: AtomId,
    atoms: &'a [NodeId],
    pub ast: &'a DativeBondAst,
    #[allow(dead_code)]
    molecule: &'a MoleculeAst,
}

impl<'a> DativeBondView<'a> {
    #[inline]
    pub fn acceptor_slot(&self) -> u8 {
        self.ast.acceptor_slot
    }

    #[inline]
    pub fn order(&self) -> &'a ValueAst {
        &self.ast.order
    }

    #[inline]
    pub fn constraints(&self) -> &'a DativeBondConstraints {
        &self.ast.constraints
    }

    /// All atoms in this dative bond (donors + acceptor), sorted by `AtomId`.
    pub fn atoms(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }

    /// Donor atoms (participants minus the acceptor slot).
    pub fn donors(&self) -> impl Iterator<Item = AtomId> + '_ {
        let acceptor_slot = self.ast.acceptor_slot as usize;
        self.atoms
            .iter()
            .enumerate()
            .filter(move |(i, _)| *i != acceptor_slot)
            .map(|(_, &n)| AtomId::from(n))
    }
}

/// Namespace accessor for aromatic-system views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct AromaticSystemViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a VarRelationSet<AromaticSystemAst>,
}

impl<'a> AromaticSystemViews<'a> {
    pub(super) fn new(molecule: &'a MoleculeAst, set: &'a VarRelationSet<AromaticSystemAst>) -> Self {
        Self { molecule, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = AromaticSystemId> {
        self.set.relation_ids().map(AromaticSystemId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = AromaticSystemView<'a>> {
        let molecule = self.molecule;
        let set = self.set;
        set.relation_ids().map(move |rid| AromaticSystemView {
            id: AromaticSystemId::from(rid),
            ast: set.data(rid),
            atoms: set.participants(rid),
            molecule,
        })
    }

    pub fn get(&self, id: AromaticSystemId) -> AromaticSystemView<'a> {
        let rid = RelationId::from(id);
        AromaticSystemView {
            id,
            ast: self.set.data(rid),
            atoms: self.set.participants(rid),
            molecule: self.molecule,
        }
    }
}

impl<'a> Index<AromaticSystemId> for AromaticSystemViews<'a> {
    type Output = AromaticSystemAst;
    fn index(&self, id: AromaticSystemId) -> &AromaticSystemAst {
        self.set.data(RelationId::from(id))
    }
}

/// Borrowed view of an aromatic system: its index, the `AromaticSystemAst`,
/// and accessors for member atoms and induced ring bonds via `atoms()` and
/// `bonds()`.
#[derive(Clone, Copy, Debug)]
pub struct AromaticSystemView<'a> {
    pub id: AromaticSystemId,
    atoms: &'a [NodeId],
    pub ast: &'a AromaticSystemAst,
    molecule: &'a MoleculeAst,
}

impl<'a> AromaticSystemView<'a> {
    #[inline]
    pub fn electrons(&self) -> &'a [ValueAst] {
        &self.ast.electrons
    }

    #[inline]
    pub fn charge(&self) -> &'a ValueAst {
        &self.ast.charge
    }

    #[inline]
    pub fn spin(&self) -> &'a SpinStateAst {
        &self.ast.spin
    }

    #[inline]
    pub fn constraints(&self) -> &'a AromaticSystemConstraints {
        &self.ast.constraints
    }

    pub fn atoms(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }

    pub fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
        self.molecule
            .raw_graph()
            .induced_edges(self.atoms)
            .map(BondId::from)
    }
}

/// Namespace accessor for multicenter-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct MulticenterBondViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a VarRelationSet<MulticenterBondAst>,
}

impl<'a> MulticenterBondViews<'a> {
    pub(super) fn new(molecule: &'a MoleculeAst, set: &'a VarRelationSet<MulticenterBondAst>) -> Self {
        Self { molecule, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = MulticenterBondId> {
        self.set.relation_ids().map(MulticenterBondId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = MulticenterBondView<'a>> {
        let molecule = self.molecule;
        let set = self.set;
        set.relation_ids().map(move |rid| MulticenterBondView {
            id: MulticenterBondId::from(rid),
            ast: set.data(rid),
            atoms: set.participants(rid),
            molecule,
        })
    }

    pub fn get(&self, id: MulticenterBondId) -> MulticenterBondView<'a> {
        let rid = RelationId::from(id);
        MulticenterBondView {
            id,
            ast: self.set.data(rid),
            atoms: self.set.participants(rid),
            molecule: self.molecule,
        }
    }
}

impl<'a> Index<MulticenterBondId> for MulticenterBondViews<'a> {
    type Output = MulticenterBondAst;
    fn index(&self, id: MulticenterBondId) -> &MulticenterBondAst {
        self.set.data(RelationId::from(id))
    }
}

/// Borrowed view of a multicenter bond: its index, member atoms via
/// `atoms()`, and underlying `MulticenterBondAst`.
#[derive(Clone, Copy, Debug)]
pub struct MulticenterBondView<'a> {
    pub id: MulticenterBondId,
    atoms: &'a [NodeId],
    pub ast: &'a MulticenterBondAst,
    #[allow(dead_code)]
    molecule: &'a MoleculeAst,
}

impl<'a> MulticenterBondView<'a> {
    #[inline]
    pub fn electrons(&self) -> &'a [ValueAst] {
        &self.ast.electrons
    }

    #[inline]
    pub fn charge(&self) -> &'a ValueAst {
        &self.ast.charge
    }

    #[inline]
    pub fn spin(&self) -> &'a SpinStateAst {
        &self.ast.spin
    }

    #[inline]
    pub fn constraints(&self) -> &'a MulticenterBondConstraints {
        &self.ast.constraints
    }

    pub fn atoms(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atoms.iter().map(|&n| AtomId::from(n))
    }
}

/// Namespace accessor for noncovalent-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct NoncovalentBondViews<'a> {
    molecule: &'a MoleculeAst,
    set: &'a FixedRelationSet<NoncovalentBondAst, 2>,
}

impl<'a> NoncovalentBondViews<'a> {
    pub(super) fn new(
        molecule: &'a MoleculeAst,
        set: &'a FixedRelationSet<NoncovalentBondAst, 2>,
    ) -> Self {
        Self { molecule, set }
    }

    pub fn count(&self) -> usize {
        self.set.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = NoncovalentBondId> {
        self.set.relation_ids().map(NoncovalentBondId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = NoncovalentBondView<'a>> {
        let molecule = self.molecule;
        let set = self.set;
        set.relation_ids().map(move |rid| NoncovalentBondView {
            id: NoncovalentBondId::from(rid),
            ast: set.data(rid),
            atoms: {
                let parts = set.participants(rid);
                [AtomId::from(parts[0]), AtomId::from(parts[1])]
            },
            molecule,
        })
    }

    pub fn get(&self, id: NoncovalentBondId) -> NoncovalentBondView<'a> {
        let rid = RelationId::from(id);
        let parts = self.set.participants(rid);
        NoncovalentBondView {
            id,
            ast: self.set.data(rid),
            atoms: [AtomId::from(parts[0]), AtomId::from(parts[1])],
            molecule: self.molecule,
        }
    }
}

impl<'a> Index<NoncovalentBondId> for NoncovalentBondViews<'a> {
    type Output = NoncovalentBondAst;
    fn index(&self, id: NoncovalentBondId) -> &NoncovalentBondAst {
        self.set.data(RelationId::from(id))
    }
}

/// Borrowed view of a noncovalent bond: the two participating atoms plus data.
#[derive(Clone, Copy, Debug)]
pub struct NoncovalentBondView<'a> {
    pub id: NoncovalentBondId,
    atoms: [AtomId; 2],
    pub ast: &'a NoncovalentBondAst,
    #[allow(dead_code)]
    molecule: &'a MoleculeAst,
}

impl<'a> NoncovalentBondView<'a> {
    #[inline]
    pub fn kind(&self) -> &'a NoncovalentBondKindAst {
        &self.ast.kind
    }

    #[inline]
    pub fn constraints(&self) -> &'a NoncovalentBondConstraints {
        &self.ast.constraints
    }

    /// The two atoms in this noncovalent interaction.
    pub fn atoms(&self) -> [AtomId; 2] {
        self.atoms
    }
}

/// Neighbor-side view of a bond: the atom on the other end (`atom`), the
/// bond index, the bond data, and the parent `MoleculeAst` for navigation
/// to the neighbor's full atom view. Yielded by `MoleculeAst::neighbors`.
#[derive(Clone, Copy, Debug)]
pub struct NeighborView<'a> {
    pub bond: BondId,
    pub atom: AtomId,
    pub ast: &'a BondAst,
    #[allow(dead_code)]
    molecule: &'a MoleculeAst,
}

impl<'a> NeighborView<'a> {
    pub(super) fn new(
        bond: BondId,
        atom: AtomId,
        ast: &'a BondAst,
        molecule: &'a MoleculeAst,
    ) -> Self {
        Self {
            bond,
            atom,
            ast,
            molecule,
        }
    }
}

/// AtomId/BondId-typed adapter over the underlying `Graph`. Holds the
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

    pub fn degree(&self, atom: AtomId) -> usize {
        self.graph.degree(NodeId::from(atom))
    }

    pub fn connected_components(&self, alg: ConnectedComponentsAlgorithm) -> Vec<Vec<AtomId>> {
        self.graph
            .connected_components(alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn biconnected_components(
        &self,
        alg: BiconnectedComponentsAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .biconnected_components(alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn shortest_cycle_through_bond(
        &self,
        bond: BondId,
        alg: ShortestCycleAlgorithm,
    ) -> Option<usize> {
        self.graph
            .shortest_cycle_through_edge(EdgeId::from(bond), alg)
    }

    pub fn shortest_cycle_through_atom(
        &self,
        atom: AtomId,
        alg: ShortestCycleAlgorithm,
    ) -> Option<usize> {
        self.graph
            .shortest_cycle_through_node(NodeId::from(atom), alg)
    }

    pub fn enumerate_cycles(
        &self,
        max_size: usize,
        alg: CycleEnumerationAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .enumerate_cycles(max_size, alg)
            .into_iter()
            .map(|c| c.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn maximum_independent_set(&self, alg: MaxIndependentSetAlgorithm) -> Vec<AtomId> {
        self.graph
            .maximum_independent_set(alg)
            .into_iter()
            .map(AtomId::from)
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
        atom_color: impl Fn(AtomId) -> C,
        alg: AutomorphismAlgorithm,
    ) -> AtomAutomorphism {
        AtomAutomorphism(
            self.graph
                .automorphisms(|n| atom_color(AtomId::from(n)), alg),
        )
    }

    pub fn subgraph_isomorphisms(
        &self,
        query: &GraphView<'_>,
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .subgraph_isomorphisms(
                query.graph,
                &mut |tn, qn| atom_match(AtomId::from(tn), AtomId::from(qn)),
                &mut |te, qe| bond_match(BondId::from(te), BondId::from(qe)),
                alg,
            )
            .into_iter()
            .map(|m| m.into_iter().map(AtomId::from).collect())
            .collect()
    }

    pub fn subgraph_isomorphisms_at(
        &self,
        query: &GraphView<'_>,
        anchor: (AtomId, AtomId),
        atom_match: &mut impl FnMut(AtomId, AtomId) -> bool,
        bond_match: &mut impl FnMut(BondId, BondId) -> bool,
        alg: SubgraphIsomorphismAlgorithm,
    ) -> Vec<Vec<AtomId>> {
        self.graph
            .subgraph_isomorphisms_at(
                query.graph,
                (NodeId::from(anchor.0), NodeId::from(anchor.1)),
                &mut |tn, qn| atom_match(AtomId::from(tn), AtomId::from(qn)),
                &mut |te, qe| bond_match(BondId::from(te), BondId::from(qe)),
                alg,
            )
            .into_iter()
            .map(|m| m.into_iter().map(AtomId::from).collect())
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
    use crate::ast::constraint::{
        AromaticValenceAst, AtomConstraint, Constraints, MulticenterValenceAst,
    };
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
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(2)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
            ],
            vec![(vec![AtomId(2)], AtomId(3), DativeBondAst::from_order(1))],
            vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemAst::default(),
            )],
            vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondAst::default(),
            )],
            vec![(
                AtomId(0),
                AtomId(3),
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
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)],
        );
    }

    #[rstest]
    fn test_atom_views_iter(molecule: MoleculeAst) {
        let views = molecule.atoms();
        let collected: Vec<(AtomId, AtomAst)> =
            views.iter().map(|v| (v.id, v.ast.clone())).collect();
        assert_eq!(
            collected,
            vec![
                (AtomId(0), AtomAst::from_element(Element::C)),
                (AtomId(1), AtomAst::from_element(Element::C)),
                (AtomId(2), AtomAst::from_element(Element::N)),
                (AtomId(3), AtomAst::from_element(Element::O)),
            ],
        );
    }

    #[rstest]
    fn test_atom_views_get(molecule: MoleculeAst) {
        let view = molecule.atoms().get(AtomId(2));
        assert_eq!(view.id, AtomId(2));
        assert_eq!(*view.ast, AtomAst::from_element(Element::N));
    }

    #[rstest]
    fn test_atom_views_index(molecule: MoleculeAst) {
        let atom: &AtomAst = &molecule.atoms()[AtomId(2)];
        assert_eq!(*atom, AtomAst::from_element(Element::N));
    }

    // --- AtomView ---

    #[rstest]
    fn test_atom_view_neighbors(molecule: MoleculeAst) {
        let view = molecule.atom(AtomId(1));
        let collected: Vec<(BondId, AtomId, BondAst)> = view
            .neighbors()
            .map(|n| (n.bond, n.atom, n.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondId(0), AtomId(0), BondAst::from_order(1)),
                (BondId(1), AtomId(2), BondAst::from_order(2)),
            ],
        );
    }

    #[rstest]
    #[case::no_incident(AtomId(3), Some(0))]
    #[case::single(AtomId(0), Some(1))]
    #[case::three_around_center(AtomId(1), Some(3))]
    #[case::double(AtomId(2), Some(2))]
    fn test_atom_view_valence(#[case] center: AtomId, #[case] expected: Option<u32>) {
        let molecule = mol!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#);
        assert_eq!(molecule.atom(center).valence(), expected);
    }

    #[rstest]
    fn test_atom_view_valence_undetermined() {
        let molecule = mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "*"]]}"#);
        assert_eq!(molecule.atom(AtomId(0)).valence(), None);
    }

    #[rstest]
    #[case::with_constraint(Some(AtomConstraint::valence(4)), ValueAst::Lit(4))]
    #[case::absent(None, ValueAst::Undetermined)]
    fn test_atom_view_valence_constraint(
        #[case] constraint: Option<AtomConstraint>,
        #[case] expected: ValueAst,
    ) {
        let mut atom = AtomAst::from_element(Element::C);
        if let Some(c) = constraint {
            atom.constraints.add(c);
        }
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(molecule.atom(AtomId(0)).constraints().valence(), expected);
    }

    #[rstest]
    #[case::donor(AtomId(0), Some(1))]
    #[case::acceptor(AtomId(1), Some(0))]
    fn test_atom_view_donated_pairs(#[case] atom: AtomId, #[case] expected: Option<u32>) {
        let molecule = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(molecule.atom(atom).donated_pairs(), expected);
    }

    #[rstest]
    fn test_atom_view_donated_pairs_constraint() {
        let mut atom = AtomAst::from_element(Element::N);
        atom.constraints.add(AtomConstraint::donated_pairs(1));
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().donated_pairs(),
            ValueAst::Lit(1),
        );
    }

    #[rstest]
    #[case::donor(AtomId(0), Some(0))]
    #[case::acceptor(AtomId(1), Some(1))]
    fn test_atom_view_accepted_pairs(#[case] atom: AtomId, #[case] expected: Option<u32>) {
        let molecule = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(molecule.atom(atom).accepted_pairs(), expected);
    }

    #[rstest]
    fn test_atom_view_accepted_pairs_constraint() {
        let mut atom = AtomAst::from_element(Element::C);
        atom.constraints.add(AtomConstraint::accepted_pairs(2));
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().accepted_pairs(),
            ValueAst::Lit(2),
        );
    }

    #[rstest]
    #[case::lit(ValueAst::Lit(2), Some(2))]
    #[case::undetermined(ValueAst::Undetermined, None)]
    fn test_atom_view_aromatic_valence(
        #[case] entry: ValueAst,
        #[case] expected: Option<u32>,
    ) {
        let molecule = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            vec![],
            vec![(
                vec![AtomId(0), AtomId(1)],
                AromaticSystemAst::new(vec![entry, ValueAst::Lit(1)]),
            )],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(molecule.atom(AtomId(0)).aromatic_valence(), expected);
    }

    #[rstest]
    fn test_atom_view_aromatic_valence_not_in_system() {
        let molecule = mol!(r#"{:atoms ["C"] :bonds []}"#);
        assert_eq!(molecule.atom(AtomId(0)).aromatic_valence(), Some(0));
    }

    #[rstest]
    #[case::in_system(AtomId(0), true)]
    #[case::not_in_system(AtomId(3), false)]
    fn test_atom_view_is_in_aromatic_system(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
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
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().aromatic_valence(),
            AromaticValenceAst::Aromatic(ValueAst::Lit(1)),
        );
    }

    #[rstest]
    #[case::single_bond(vec![(vec![AtomId(0), AtomId(1)], vec![ValueAst::Lit(2), ValueAst::Lit(2)])], Some(2))]
    #[case::two_bonds(
        vec![
            (vec![AtomId(0), AtomId(1)], vec![ValueAst::Lit(2), ValueAst::Lit(2)]),
            (vec![AtomId(0), AtomId(2)], vec![ValueAst::Lit(1), ValueAst::Lit(1)]),
        ],
        Some(3),
    )]
    #[case::undetermined_aborts(
        vec![(vec![AtomId(0), AtomId(1)], vec![ValueAst::Undetermined, ValueAst::Lit(2)])],
        None,
    )]
    fn test_atom_view_multicenter_valence(
        #[case] bonds: Vec<(Vec<AtomId>, Vec<ValueAst>)>,
        #[case] expected: Option<u32>,
    ) {
        let multicenter: Vec<_> = bonds
            .into_iter()
            .map(|(parts, electrons)| (parts, MulticenterBondAst::new(electrons)))
            .collect();
        let molecule = MoleculeAst::from_parts(
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
        assert_eq!(molecule.atom(AtomId(0)).multicenter_valence(), expected);
    }

    #[rstest]
    fn test_atom_view_multicenter_valence_constraint() {
        let mut atom = AtomAst::from_element(Element::C);
        atom.constraints.add(AtomConstraint::multicenter_valence(
            MulticenterValenceAst::Multicenter(ValueAst::Lit(2)),
        ));
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().multicenter_valence(),
            MulticenterValenceAst::Multicenter(ValueAst::Lit(2)),
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
            vec![BondId(0), BondId(1), BondId(2)],
        );
    }

    #[rstest]
    fn test_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(BondId, [AtomId; 2], BondAst)> = molecule
            .bonds()
            .iter()
            .map(|v| (v.id, v.atoms(), v.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondId(0), [AtomId(0), AtomId(1)], BondAst::from_order(1)),
                (BondId(1), [AtomId(1), AtomId(2)], BondAst::from_order(2)),
                (BondId(2), [AtomId(2), AtomId(3)], BondAst::from_order(1)),
            ],
        );
    }

    #[rstest]
    fn test_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.bonds().get(BondId(1));
        assert_eq!(view.id, BondId(1));
        assert_eq!(view.atoms(), [AtomId(1), AtomId(2)]);
        assert_eq!(*view.ast, BondAst::from_order(2));
    }

    #[rstest]
    fn test_bond_views_index(molecule: MoleculeAst) {
        let bond: &BondAst = &molecule.bonds()[BondId(1)];
        assert_eq!(*bond, BondAst::from_order(2));
    }

    // --- BondView ---

    #[rstest]
    fn test_bond_view_atoms(molecule: MoleculeAst) {
        assert_eq!(molecule.bond(BondId(1)).atoms(), [AtomId(1), AtomId(2)]);
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
            vec![DativeBondId(0)],
        );
    }

    #[rstest]
    fn test_dative_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(DativeBondId, AtomId, DativeBondAst)> = molecule
            .dative_bonds()
            .iter()
            .map(|v| (v.id, v.acceptor, v.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                DativeBondId(0),
                AtomId(3),
                DativeBondAst::from_order(1).with_acceptor_slot(1),
            )],
        );
    }

    #[rstest]
    fn test_dative_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.dative_bonds().get(DativeBondId(0));
        assert_eq!(view.id, DativeBondId(0));
        assert_eq!(view.acceptor, AtomId(3));
    }

    #[rstest]
    fn test_dative_bond_views_index(molecule: MoleculeAst) {
        let dative: &DativeBondAst = &molecule.dative_bonds()[DativeBondId(0)];
        assert_eq!(dative.order, ValueAst::Lit(1));
    }

    // --- DativeBondView ---

    #[rstest]
    fn test_dative_bond_view_atoms(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .dative_bond(DativeBondId(0))
                .atoms()
                .collect::<Vec<_>>(),
            vec![AtomId(2), AtomId(3)],
        );
    }

    #[rstest]
    fn test_dative_bond_view_donors(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .dative_bond(DativeBondId(0))
                .donors()
                .collect::<Vec<_>>(),
            vec![AtomId(2)],
        );
    }

    #[rstest]
    fn test_dative_bond_view_acceptor(molecule: MoleculeAst) {
        assert_eq!(molecule.dative_bond(DativeBondId(0)).acceptor, AtomId(3));
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
            vec![AromaticSystemId(0)],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(AromaticSystemId, Vec<AtomId>)> = molecule
            .aromatic_systems()
            .iter()
            .map(|v| (v.id, v.atoms().collect()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                AromaticSystemId(0),
                vec![AtomId(0), AtomId(1), AtomId(2)]
            )],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_get(molecule: MoleculeAst) {
        let view = molecule.aromatic_systems().get(AromaticSystemId(0));
        assert_eq!(view.id, AromaticSystemId(0));
        assert_eq!(
            view.atoms().collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_index(molecule: MoleculeAst) {
        let _: &AromaticSystemAst = &molecule.aromatic_systems()[AromaticSystemId(0)];
    }

    // --- AromaticSystemView ---

    #[rstest]
    fn test_aromatic_system_view_atoms(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .atoms()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_aromatic_system_view_bonds(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .bonds()
                .collect::<Vec<_>>(),
            vec![BondId(0), BondId(1)],
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
            vec![MulticenterBondId(0)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(MulticenterBondId, Vec<AtomId>)> = molecule
            .multicenter_bonds()
            .iter()
            .map(|v| (v.id, v.atoms().collect()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                MulticenterBondId(0),
                vec![AtomId(0), AtomId(1), AtomId(2)],
            )],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.multicenter_bonds().get(MulticenterBondId(0));
        assert_eq!(view.id, MulticenterBondId(0));
        assert_eq!(
            view.atoms().collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_index(molecule: MoleculeAst) {
        let _: &MulticenterBondAst = &molecule.multicenter_bonds()[MulticenterBondId(0)];
    }

    // --- MulticenterBondView ---

    #[rstest]
    fn test_multicenter_bond_view_atoms(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .multicenter_bond(MulticenterBondId(0))
                .atoms()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
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
            vec![NoncovalentBondId(0)],
        );
    }

    #[rstest]
    fn test_noncovalent_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(NoncovalentBondId, [AtomId; 2], NoncovalentBondAst)> = molecule
            .noncovalent_bonds()
            .iter()
            .map(|v| (v.id, v.atoms(), v.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                NoncovalentBondId(0),
                [AtomId(0), AtomId(3)],
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
        );
    }

    #[rstest]
    fn test_noncovalent_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.noncovalent_bonds().get(NoncovalentBondId(0));
        assert_eq!(view.id, NoncovalentBondId(0));
        assert_eq!(view.atoms(), [AtomId(0), AtomId(3)]);
    }

    #[rstest]
    fn test_noncovalent_bond_views_index(molecule: MoleculeAst) {
        let _: &NoncovalentBondAst = &molecule.noncovalent_bonds()[NoncovalentBondId(0)];
    }

    // --- NoncovalentBondView ---

    #[rstest]
    fn test_noncovalent_bond_view_atoms(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .noncovalent_bond(NoncovalentBondId(0))
                .atoms(),
            [AtomId(0), AtomId(3)],
        );
    }

    // --- NeighborView ---

    #[rstest]
    fn test_neighbor_view_fields(molecule: MoleculeAst) {
        let collected: Vec<(BondId, AtomId, BondAst)> = molecule
            .neighbors(AtomId(2))
            .map(|n| (n.bond, n.atom, n.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondId(1), AtomId(1), BondAst::from_order(2)),
                (BondId(2), AtomId(3), BondAst::from_order(1)),
            ],
        );
    }
}
