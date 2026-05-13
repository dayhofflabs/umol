//! Molecule structural AST.

use std::collections::HashSet;
use std::iter;
use std::ops::Index;
use std::sync::Arc;

pub use builder::MoleculeBuilder;
use umol_graph_core::{FixedRelationSet, Graph, NodeId, RelationId, VarRelationSet};

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::constraint::{Constraint, Constraints};
use super::dative::DativeBondAst;
use super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;
use super::rings::{self, RingFamily, RingSet};
use super::subgraph::MoleculeSubgraph;
use super::views::{
    AromaticSystemView, AromaticSystemViews, AtomView, AtomViewMut, AtomViews, BondView,
    BondViewMut, BondViews, DativeBondView, DativeBondViews, GraphView, MulticenterBondView,
    MulticenterBondViews, NeighborView, NoncovalentBondView, NoncovalentBondViews,
};

mod builder;
mod rewrite;

/// Molecule AST: structural representation of a molecule (ground or pattern).
///
/// Topology and per-atom/bond data are `Arc`-shared (copy-on-write). The AST
/// itself only allows attribute mutation (`atom_mut`, `bond_mut`); structural
/// edits go through `MoleculeBuilder` via [`MoleculeAst::edit`].
///
/// Carries a single-slot last-request-wins ring cache (`rings(family,
/// max_ring_size)`). Caching is sound because topology is invariant across
/// every in-place mutation path (only attribute fields change). Structural
/// edits go through the builder, which produces a fresh `MoleculeAst` with
/// an empty cache.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeAst {
    graph: Graph,
    atoms: Arc<Vec<AtomAst>>,
    bonds: Arc<Vec<BondAst>>,
    dative_bonds: Arc<VarRelationSet<DativeBondAst>>,
    aromatic_systems: Arc<VarRelationSet<AromaticSystemAst>>,
    multicenter_bonds: Arc<VarRelationSet<MulticenterBondAst>>,
    noncovalent_bonds: Arc<FixedRelationSet<NoncovalentBondAst, 2>>,
    constraints: Constraints,
    /// Single-slot ring cache; never participates in `PartialEq`/`Hash`.
    ring_cache: RingCache,
}

#[derive(Clone, Debug, Default)]
struct RingCache(Option<Box<RingCacheEntry>>);

#[derive(Clone, Debug)]
struct RingCacheEntry {
    family: RingFamily,
    max_ring_size: usize,
    rings: RingSet,
}

impl PartialEq for RingCache {
    fn eq(&self, _other: &Self) -> bool {
        true
    }
}

impl Eq for RingCache {}

impl Default for MoleculeAst {
    fn default() -> Self {
        Self {
            graph: Graph::default(),
            atoms: Arc::new(Vec::new()),
            bonds: Arc::new(Vec::new()),
            dative_bonds: Arc::new(VarRelationSet::default()),
            aromatic_systems: Arc::new(VarRelationSet::default()),
            multicenter_bonds: Arc::new(VarRelationSet::default()),
            noncovalent_bonds: Arc::new(FixedRelationSet::default()),
            constraints: Constraints::new(),
            ring_cache: RingCache::default(),
        }
    }
}

impl AsRef<MoleculeAst> for MoleculeAst {
    fn as_ref(&self) -> &MoleculeAst {
        self
    }
}

impl MoleculeAst {
    /// Empty molecule: zero atoms, zero bonds, zero relations, zero
    /// constraints. Mirrors `Vec::new()` / `HashMap::new()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Common shape: atoms plus pairwise bonds, no relations or constraints.
    pub fn from_atoms_and_bonds(
        atoms: Vec<AtomAst>,
        bonds: Vec<(AtomId, AtomId, BondAst)>,
    ) -> Self {
        Self::from_parts(
            atoms,
            bonds,
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        )
    }

    /// Start an empty `MoleculeBuilder` for fluent / programmatic
    /// construction. Use [`MoleculeAst::edit`] to start from an existing
    /// molecule.
    pub fn builder() -> MoleculeBuilder {
        Self::new().edit()
    }

    /// Full structural constructor: every entity-type vector is supplied
    /// directly. The escape hatch when the molecule has relations or
    /// molecule-level constraints; tests covering all entity types route
    /// through here.
    #[allow(clippy::too_many_arguments)]
    pub fn from_parts(
        atoms: Vec<AtomAst>,
        bonds: Vec<(AtomId, AtomId, BondAst)>,
        dative: Vec<(Vec<AtomId>, AtomId, DativeBondAst)>,
        aromatic: Vec<(Vec<AtomId>, AromaticSystemAst)>,
        multicenter: Vec<(Vec<AtomId>, MulticenterBondAst)>,
        noncovalent: Vec<(AtomId, AtomId, NoncovalentBondAst)>,
        constraints: Constraints,
    ) -> Self {
        let node_count = atoms.len();
        let edges: Vec<[u32; 2]> = bonds.iter().map(|(s, t, _)| [s.0, t.0]).collect();
        let bond_data: Vec<BondAst> = bonds.into_iter().map(|(_, _, d)| d).collect();
        let graph = Graph::new(node_count, &edges);

        let dative_bonds = VarRelationSet::new(
            dative
                .into_iter()
                .map(|(donors, acceptor, mut d)| {
                    let acceptor_node = NodeId::from(acceptor);
                    let mut participants: Vec<NodeId> = donors
                        .into_iter()
                        .map(NodeId::from)
                        .chain(iter::once(acceptor_node))
                        .collect();
                    participants.sort_unstable();
                    let slot = participants
                        .iter()
                        .position(|&n| n == acceptor_node)
                        .expect("acceptor must appear in participants");
                    d.acceptor_slot = slot as u8;
                    (participants, d)
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
            ring_cache: RingCache::default(),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn from_arcs(
        graph: Graph,
        atoms: Arc<Vec<AtomAst>>,
        bonds: Arc<Vec<BondAst>>,
        dative_bonds: Arc<VarRelationSet<DativeBondAst>>,
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
            ring_cache: RingCache::default(),
        }
    }

    /// AtomId/BondId-typed adapter exposing the pure-graph algorithms.
    pub fn graph(&self) -> GraphView<'_> {
        GraphView::new(&self.graph)
    }

    /// Raw underlying graph with `NodeId` / `EdgeId` types. Escape hatch
    /// for code that needs the graph-core API directly; use [`Self::graph`]
    /// for AtomId/BondId-typed access.
    pub fn raw_graph(&self) -> &Graph {
        &self.graph
    }

    pub fn neighbors(&self, atom: AtomId) -> impl Iterator<Item = NeighborView<'_>> {
        let bonds = &self.bonds;
        self.graph
            .neighbors(NodeId::from(atom))
            .iter()
            .map(move |n| {
                NeighborView::new(
                    BondId::from(n.edge),
                    AtomId::from(n.node),
                    &bonds[n.edge.index()],
                    self,
                )
            })
    }

    pub fn atoms(&self) -> AtomViews<'_> {
        AtomViews::new(self, &self.atoms)
    }

    pub fn atom(&self, idx: AtomId) -> AtomView<'_> {
        self.atoms().get(idx)
    }

    pub fn bonds(&self) -> BondViews<'_> {
        BondViews::new(self, &self.bonds)
    }

    pub fn bond(&self, idx: BondId) -> BondView<'_> {
        self.bonds().get(idx)
    }

    pub fn dative_bonds(&self) -> DativeBondViews<'_> {
        DativeBondViews::new(self, &self.dative_bonds)
    }

    pub fn dative_bond(&self, idx: DativeBondId) -> DativeBondView<'_> {
        self.dative_bonds().get(idx)
    }

    pub fn aromatic_systems(&self) -> AromaticSystemViews<'_> {
        AromaticSystemViews::new(self, &self.aromatic_systems)
    }

    pub fn aromatic_system(&self, idx: AromaticSystemId) -> AromaticSystemView<'_> {
        self.aromatic_systems().get(idx)
    }

    pub fn multicenter_bonds(&self) -> MulticenterBondViews<'_> {
        MulticenterBondViews::new(self, &self.multicenter_bonds)
    }

    pub fn multicenter_bond(&self, idx: MulticenterBondId) -> MulticenterBondView<'_> {
        self.multicenter_bonds().get(idx)
    }

    pub fn noncovalent_bonds(&self) -> NoncovalentBondViews<'_> {
        NoncovalentBondViews::new(self, &self.noncovalent_bonds)
    }

    pub fn noncovalent_bond(&self, idx: NoncovalentBondId) -> NoncovalentBondView<'_> {
        self.noncovalent_bonds().get(idx)
    }

    pub fn connecting_bond(&self, a: AtomId, b: AtomId) -> Option<BondId> {
        self.graph
            .find_edge(NodeId::from(a), NodeId::from(b))
            .map(BondId::from)
    }

    pub fn connecting_dative_bond(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<DativeBondId> {
        let atoms: HashSet<AtomId> = atoms.into_iter().collect();
        let &first = atoms.iter().next()?;
        self.dative_bonds_incident(first).find(|&idx| {
            let v = self.dative_bond(idx);
            let parts: HashSet<AtomId> = v.atoms().collect();
            parts == atoms
        })
    }

    pub fn connecting_noncovalent_bond(
        &self,
        a: AtomId,
        b: AtomId,
    ) -> Option<NoncovalentBondId> {
        self.noncovalent_bonds_incident(a).find(|&idx| {
            let v = self.noncovalent_bond(idx);
            (v.atoms()[0] == a && v.atoms()[1] == b) || (v.atoms()[0] == b && v.atoms()[1] == a)
        })
    }

    pub fn connecting_aromatic_system(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<AromaticSystemId> {
        let atoms: HashSet<AtomId> = atoms.into_iter().collect();
        let &first = atoms.iter().next()?;
        self.aromatic_systems_incident(first).find(|&idx| {
            let v = self.aromatic_system(idx);
            let parts: HashSet<AtomId> = v.atoms().collect();
            parts == atoms
        })
    }

    pub fn connecting_multicenter_bond(
        &self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Option<MulticenterBondId> {
        let atoms: HashSet<AtomId> = atoms.into_iter().collect();
        let &first = atoms.iter().next()?;
        self.multicenter_bonds_incident(first).find(|&idx| {
            let v = self.multicenter_bond(idx);
            let parts: HashSet<AtomId> = v.atoms().collect();
            parts == atoms
        })
    }

    pub fn dative_bonds_incident(&self, atom: AtomId) -> impl Iterator<Item = DativeBondId> + '_ {
        self.dative_bonds
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| DativeBondId::from(rid))
    }

    pub fn aromatic_systems_incident(
        &self,
        atom: AtomId,
    ) -> impl Iterator<Item = AromaticSystemId> + '_ {
        self.aromatic_systems
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| AromaticSystemId::from(rid))
    }

    pub fn multicenter_bonds_incident(
        &self,
        atom: AtomId,
    ) -> impl Iterator<Item = MulticenterBondId> + '_ {
        self.multicenter_bonds
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| MulticenterBondId::from(rid))
    }

    pub fn noncovalent_bonds_incident(
        &self,
        atom: AtomId,
    ) -> impl Iterator<Item = NoncovalentBondId> + '_ {
        self.noncovalent_bonds
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| NoncovalentBondId::from(rid))
    }

    pub fn induced_subgraph(&self, atoms: &[AtomId]) -> MoleculeSubgraph {
        let keep: HashSet<AtomId> = atoms.iter().copied().collect();
        let remove_atoms: Vec<AtomId> = (0..self.atoms().count())
            .map(AtomId::from)
            .filter(|a| !keep.contains(a))
            .collect();
        let remove_bonds: Vec<BondId> = self
            .bonds()
            .iter()
            .filter(|b| !keep.contains(&b.atoms()[0]) || !keep.contains(&b.atoms()[1]))
            .map(|b| b.id)
            .collect();
        let mut builder = self.edit();
        let remap = builder.remove(&remove_atoms, &remove_bonds);
        let ast = builder.build();

        let atom_map: Vec<AtomId> = (0..self.atoms().count())
            .map(AtomId::from)
            .filter(|&a| remap.atom(a).is_some())
            .collect();
        let bond_map: Vec<BondId> = (0..self.bonds().count())
            .map(BondId::from)
            .filter(|&b| remap.bond(b).is_some())
            .collect();
        let dative_bond_map: Vec<DativeBondId> = (0..self.dative_bonds().count())
            .map(DativeBondId::from)
            .filter(|&d| remap.dative_bond(d).is_some())
            .collect();
        let aromatic_system_map: Vec<AromaticSystemId> = (0..self.aromatic_systems().count())
            .map(AromaticSystemId::from)
            .filter(|&a| remap.aromatic_system(a).is_some())
            .collect();
        let multicenter_bond_map: Vec<MulticenterBondId> = (0..self.multicenter_bonds().count())
            .map(MulticenterBondId::from)
            .filter(|&m| remap.multicenter_bond(m).is_some())
            .collect();
        let noncovalent_bond_map: Vec<NoncovalentBondId> = (0..self.noncovalent_bonds().count())
            .map(NoncovalentBondId::from)
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

    pub fn induced_bonds(&self, atoms: &[AtomId]) -> Vec<BondId> {
        let mut nodes: Vec<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        nodes.sort_unstable();
        self.graph
            .induced_edges(&nodes)
            .map(BondId::from)
            .collect()
    }

    pub fn induced_dative_bonds(&self, atoms: &[AtomId]) -> Vec<DativeBondId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.dative_bonds
            .relation_ids()
            .filter(|&rid| {
                self.dative_bonds
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(DativeBondId::from)
            .collect()
    }

    pub fn induced_aromatic_systems(&self, atoms: &[AtomId]) -> Vec<AromaticSystemId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.aromatic_systems
            .relation_ids()
            .filter(|&rid| {
                self.aromatic_systems
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(AromaticSystemId::from)
            .collect()
    }

    pub fn induced_multicenter_bonds(&self, atoms: &[AtomId]) -> Vec<MulticenterBondId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.multicenter_bonds
            .relation_ids()
            .filter(|&rid| {
                self.multicenter_bonds
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(MulticenterBondId::from)
            .collect()
    }

    pub fn induced_noncovalent_bonds(&self, atoms: &[AtomId]) -> Vec<NoncovalentBondId> {
        let set: HashSet<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        self.noncovalent_bonds
            .relation_ids()
            .filter(|&rid| {
                self.noncovalent_bonds
                    .participants(rid)
                    .iter()
                    .all(|p| set.contains(p))
            })
            .map(NoncovalentBondId::from)
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

    /// Cached ring enumeration. Single-slot last-request-wins: a call with
    /// different `(family, max_ring_size)` than the previous one replaces the
    /// stored result. Topology is invariant across in-place mutation paths,
    /// so the cache stays valid for the molecule's lifetime; structural
    /// edits go through the builder, which produces a fresh `MoleculeAst`.
    pub fn rings(&mut self, family: RingFamily, max_ring_size: usize) -> &RingSet {
        let stale = match self.ring_cache.0.as_deref() {
            Some(entry) => entry.family != family || entry.max_ring_size != max_ring_size,
            None => true,
        };
        if stale {
            let rings = rings::enumerate_rings(&self.graph, family, max_ring_size, |_| true);
            self.ring_cache.0 = Some(Box::new(RingCacheEntry {
                family,
                max_ring_size,
                rings,
            }));
        }
        &self.ring_cache.0.as_deref().unwrap().rings
    }

    /// Uncached ring enumeration with a custom atom filter. The cache slot
    /// is not consulted or written; callers that want to amortize across
    /// queries should use [`MoleculeAst::rings`] and filter at use-time.
    pub fn enumerate_rings(
        &self,
        family: RingFamily,
        max_ring_size: usize,
        atom_filter: impl Fn(AtomId) -> bool,
    ) -> RingSet {
        rings::enumerate_rings(&self.graph, family, max_ring_size, atom_filter)
    }

    pub fn atom_mut(&mut self, idx: AtomId) -> AtomViewMut<'_> {
        let ast = &mut Arc::make_mut(&mut self.atoms)[idx.index()];
        AtomViewMut { id: idx, ast }
    }

    pub fn atoms_mut(&mut self) -> impl Iterator<Item = &mut AtomAst> {
        Arc::make_mut(&mut self.atoms).iter_mut()
    }

    pub fn bond_mut(&mut self, idx: BondId) -> BondViewMut<'_> {
        let [s, t] = self.graph.edge_endpoints(idx.into());
        let data = &mut Arc::make_mut(&mut self.bonds)[idx.index()];
        BondViewMut::new(idx, [AtomId::from(s), AtomId::from(t)], data)
    }

    pub fn bonds_mut(&mut self) -> impl Iterator<Item = &mut BondAst> {
        Arc::make_mut(&mut self.bonds).iter_mut()
    }

    pub fn dative_bond_mut(&mut self, idx: DativeBondId) -> &mut DativeBondAst {
        Arc::make_mut(&mut self.dative_bonds).data_mut(RelationId::from(idx))
    }

    pub fn dative_bonds_mut(&mut self) -> impl Iterator<Item = &mut DativeBondAst> {
        Arc::make_mut(&mut self.dative_bonds).data_iter_mut()
    }

    pub fn aromatic_system_mut(&mut self, idx: AromaticSystemId) -> &mut AromaticSystemAst {
        Arc::make_mut(&mut self.aromatic_systems).data_mut(RelationId::from(idx))
    }

    pub fn aromatic_systems_mut(&mut self) -> impl Iterator<Item = &mut AromaticSystemAst> {
        Arc::make_mut(&mut self.aromatic_systems).data_iter_mut()
    }

    pub fn multicenter_bond_mut(&mut self, idx: MulticenterBondId) -> &mut MulticenterBondAst {
        Arc::make_mut(&mut self.multicenter_bonds).data_mut(RelationId::from(idx))
    }

    pub fn multicenter_bonds_mut(&mut self) -> impl Iterator<Item = &mut MulticenterBondAst> {
        Arc::make_mut(&mut self.multicenter_bonds).data_iter_mut()
    }

    pub fn noncovalent_bond_mut(&mut self, idx: NoncovalentBondId) -> &mut NoncovalentBondAst {
        Arc::make_mut(&mut self.noncovalent_bonds).data_mut(RelationId::from(idx))
    }

    pub fn noncovalent_bonds_mut(&mut self) -> impl Iterator<Item = &mut NoncovalentBondAst> {
        Arc::make_mut(&mut self.noncovalent_bonds).data_iter_mut()
    }

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
    /// `MulticenterBond` / `NoncovalentBond` entries. The order of inserted
    /// entries in `self.constraints` is unspecified.
    pub fn lift_constraints(&mut self) {
        let atom_count = self.atoms().count();
        let bond_count = self.bonds().count();
        let dative_count = self.dative_bonds().count();
        let aromatic_count = self.aromatic_systems().count();
        let multicenter_count = self.multicenter_bonds().count();
        let noncovalent_count = self.noncovalent_bonds().count();

        let mut additions: Vec<Constraint> = Vec::new();
        for i in 0..atom_count {
            let idx = AtomId::from(i);
            for c in self.atom_mut(idx).ast.constraints.take() {
                additions.push(Constraint::Atom(idx, c));
            }
        }
        for i in 0..bond_count {
            let idx = BondId::from(i);
            for c in self.bond_mut(idx).ast.constraints.take() {
                additions.push(Constraint::Bond(idx, c));
            }
        }
        for i in 0..dative_count {
            let idx = DativeBondId::from(i);
            for c in self.dative_bond_mut(idx).constraints.take() {
                additions.push(Constraint::DativeBond(idx, c));
            }
        }
        for i in 0..aromatic_count {
            let idx = AromaticSystemId::from(i);
            for c in self.aromatic_system_mut(idx).constraints.take() {
                additions.push(Constraint::AromaticSystem(idx, c));
            }
        }
        for i in 0..multicenter_count {
            let idx = MulticenterBondId::from(i);
            for c in self.multicenter_bond_mut(idx).constraints.take() {
                additions.push(Constraint::MulticenterBond(idx, c));
            }
        }
        for i in 0..noncovalent_count {
            let idx = NoncovalentBondId::from(i);
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
                    self.atom_mut(idx).ast.constraints.add(inner);
                }
                Constraint::Bond(idx, inner) => {
                    self.bond_mut(idx).ast.constraints.add(inner);
                }
                Constraint::DativeBond(idx, inner) => {
                    self.dative_bond_mut(idx).constraints.add(inner);
                }
                Constraint::AromaticSystem(idx, inner) => {
                    self.aromatic_system_mut(idx).constraints.add(inner);
                }
                Constraint::MulticenterBond(idx, inner) => {
                    self.multicenter_bond_mut(idx).constraints.add(inner);
                }
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

impl Index<AtomId> for MoleculeAst {
    type Output = AtomAst;
    fn index(&self, idx: AtomId) -> &AtomAst {
        &self.atoms[idx.index()]
    }
}

impl Index<BondId> for MoleculeAst {
    type Output = BondAst;
    fn index(&self, idx: BondId) -> &BondAst {
        &self.bonds[idx.index()]
    }
}

impl Index<DativeBondId> for MoleculeAst {
    type Output = DativeBondAst;
    fn index(&self, idx: DativeBondId) -> &DativeBondAst {
        self.dative_bonds.data(RelationId::from(idx))
    }
}

impl Index<AromaticSystemId> for MoleculeAst {
    type Output = AromaticSystemAst;
    fn index(&self, idx: AromaticSystemId) -> &AromaticSystemAst {
        self.aromatic_systems.data(RelationId::from(idx))
    }
}

impl Index<MulticenterBondId> for MoleculeAst {
    type Output = MulticenterBondAst;
    fn index(&self, idx: MulticenterBondId) -> &MulticenterBondAst {
        self.multicenter_bonds.data(RelationId::from(idx))
    }
}

impl Index<NoncovalentBondId> for MoleculeAst {
    type Output = NoncovalentBondAst;
    fn index(&self, idx: NoncovalentBondId) -> &NoncovalentBondAst {
        self.noncovalent_bonds.data(RelationId::from(idx))
    }
}

#[cfg(test)]
mod tests;
