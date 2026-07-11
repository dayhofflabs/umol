//! Molecule structural AST.

use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::iter;
use std::ops::Index;
use std::sync::{Arc, OnceLock};

pub use build::MoleculeBuilder;
pub use editor::MoleculeEditor;
pub use fragment::{Fragment, Port, PortArg};
pub use spec::{AtomArg, MoleculeSpec, MoleculeSpecTerm};
use umol_graph_core::{
    Correspondence, EdgeId, FixedRelationSet, FixedVarBirelationSet, Graph, NodeId, Ordered,
    RelationId, RelationParticipant, Remapping, UnionFind, Unordered, VarRelationSet,
};

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::constraint::{Constraint, Constraints, MoleculeConstraint, RelationalConstraint};
use super::correspondence::MoleculeCorrespondence;
use super::dative::DativeBondAst;
use super::edit::{AtomHandle, BondHandle, Edit};
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigand;
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;
use super::remap::IdRemapping;
use super::ring::{RingFamily, RingSet};
use super::stereo::{StereoAtomAst, StereoBondAst};
use super::traits::Lattice;
use super::view::{
    AromaticSystemView, AromaticSystemViews, AtomView, AtomViewMut, AtomViews, BondView,
    BondViewMut, BondViews, DativeBondView, DativeBondViews, GraphView, MulticenterBondView,
    MulticenterBondViews, NeighborView, NoncovalentBondView, NoncovalentBondViews, StereoAtomView,
    StereoAtomViews, StereoBondView, StereoBondViews,
};

mod build;
mod editor;
mod fragment;
mod pushout;
pub mod spec;
pub(super) mod transact;

/// Molecule AST: structural representation of a molecule (ground or pattern).
///
/// Topology and per-atom/bond data are `Arc`-shared (copy-on-write). The AST
/// itself only allows attribute mutation (`atom_mut`, `bond_mut`); structural
/// edits go through `MoleculeEditor` via [`MoleculeAst::edit`].
///
/// Carries a single-slot canonical-rings cache (`OnceLock<RingSet>`) populated
/// lazily on the first call to [`MoleculeAst::rings`]. The cache stores
/// Vismara relevant cycles up to max ring size 22; non-canonical enumeration
/// goes through [`MoleculeAst::rings_with`], which is uncached and returns
/// owned. Topology is invariant across in-place attribute mutation, so the
/// cache remains valid for the molecule's lifetime; structural edits go
/// through the builder, which produces a fresh `MoleculeAst` with an empty
/// cache. The cache slot is excluded from `PartialEq` / `Hash` so identity
/// is independent of cache state.
#[derive(Debug, Default)]
pub struct MoleculeAst {
    graph: Graph,
    atoms: Arc<Vec<AtomAst>>,
    bonds: Arc<Vec<BondAst>>,
    dative_bonds: Arc<FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondAst>>,
    aromatic_systems: Arc<VarRelationSet<NodeId, Unordered, AromaticSystemAst>>,
    multicenter_bonds: Arc<VarRelationSet<NodeId, Unordered, MulticenterBondAst>>,
    noncovalent_bonds: Arc<FixedRelationSet<NodeId, Unordered, NoncovalentBondAst, 2>>,
    stereo_atoms:
        Arc<FixedVarBirelationSet<NodeId, Ordered, 1, StereoLigand, Ordered, StereoAtomAst>>,
    stereo_bonds:
        Arc<FixedVarBirelationSet<EdgeId, Ordered, 1, StereoLigand, Ordered, StereoBondAst>>,
    constraints: Constraints,
    rings_cache: OnceLock<RingSet>,
}

impl Clone for MoleculeAst {
    fn clone(&self) -> Self {
        Self {
            graph: self.graph.clone(),
            atoms: self.atoms.clone(),
            bonds: self.bonds.clone(),
            dative_bonds: self.dative_bonds.clone(),
            aromatic_systems: self.aromatic_systems.clone(),
            multicenter_bonds: self.multicenter_bonds.clone(),
            noncovalent_bonds: self.noncovalent_bonds.clone(),
            stereo_atoms: self.stereo_atoms.clone(),
            stereo_bonds: self.stereo_bonds.clone(),
            constraints: self.constraints.clone(),
            rings_cache: OnceLock::new(),
        }
    }
}

impl PartialEq for MoleculeAst {
    fn eq(&self, other: &Self) -> bool {
        self.graph == other.graph
            && self.atoms == other.atoms
            && self.bonds == other.bonds
            && self.dative_bonds == other.dative_bonds
            && self.aromatic_systems == other.aromatic_systems
            && self.multicenter_bonds == other.multicenter_bonds
            && self.noncovalent_bonds == other.noncovalent_bonds
            && self.stereo_atoms == other.stereo_atoms
            && self.stereo_bonds == other.stereo_bonds
            && self.constraints == other.constraints
    }
}

impl Eq for MoleculeAst {}

impl AsRef<MoleculeAst> for MoleculeAst {
    fn as_ref(&self) -> &MoleculeAst {
        self
    }
}

/// Constructor input for [`MoleculeAst::from_parts`].
#[derive(Debug, Default, Clone)]
pub struct MoleculeParts {
    pub atoms: Vec<AtomAst>,
    pub bonds: Vec<(AtomId, AtomId, BondAst)>,
    pub dative: Vec<(Vec<AtomId>, AtomId, DativeBondAst)>,
    pub aromatic: Vec<(Vec<AtomId>, AromaticSystemAst)>,
    pub multicenter: Vec<(Vec<AtomId>, MulticenterBondAst)>,
    pub noncovalent: Vec<(AtomId, AtomId, NoncovalentBondAst)>,
    pub stereo_atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)>,
    pub stereo_bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondAst)>,
    pub constraints: Constraints,
}

impl MoleculeAst {
    /// Empty molecule: zero atoms, zero bonds, zero relations, zero
    /// constraints. Mirrors `Vec::new()` / `HashMap::new()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Start an empty `MoleculeEditor` for fluent / programmatic
    /// construction. Use [`MoleculeAst::edit`] to start from an existing
    /// molecule.
    pub fn builder() -> MoleculeBuilder {
        MoleculeBuilder::new()
    }

    /// Full structural constructor from a flat [`MoleculeParts`]: every
    /// entity-type field is supplied directly. The topology-only case fills
    /// just `atoms` and `bonds`; relations and molecule-level constraints go
    /// in the remaining fields.
    pub fn from_parts(parts: MoleculeParts) -> Self {
        let MoleculeParts {
            atoms,
            bonds,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            stereo_atoms,
            stereo_bonds,
            constraints,
        } = parts;
        let node_count = atoms.len();
        let edges: Vec<[u32; 2]> = bonds
            .iter()
            .map(|(first, second, _)| [first.0, second.0])
            .collect();
        let bond_data: Vec<BondAst> = bonds.into_iter().map(|(_, _, d)| d).collect();
        let graph = Graph::new(node_count, &edges);

        let dative_bonds = FixedVarBirelationSet::new(
            dative
                .into_iter()
                .map(|(donors, acceptor, d)| {
                    (
                        [NodeId::from(acceptor)],
                        donors.into_iter().map(NodeId::from).collect(),
                        d,
                    )
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

        let stereo_atoms = FixedVarBirelationSet::new(
            stereo_atoms
                .into_iter()
                .map(|(site, ligands, d)| ([NodeId::from(site)], ligands, d))
                .collect(),
        );

        let stereo_bonds = FixedVarBirelationSet::new(
            stereo_bonds
                .into_iter()
                .map(|(site, ligands, d)| ([EdgeId::from(site)], ligands, d))
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
            stereo_atoms: Arc::new(stereo_atoms),
            stereo_bonds: Arc::new(stereo_bonds),
            constraints,
            rings_cache: OnceLock::new(),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn from_arcs(
        graph: Graph,
        atoms: Arc<Vec<AtomAst>>,
        bonds: Arc<Vec<BondAst>>,
        dative_bonds: Arc<
            FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondAst>,
        >,
        aromatic_systems: Arc<VarRelationSet<NodeId, Unordered, AromaticSystemAst>>,
        multicenter_bonds: Arc<VarRelationSet<NodeId, Unordered, MulticenterBondAst>>,
        noncovalent_bonds: Arc<FixedRelationSet<NodeId, Unordered, NoncovalentBondAst, 2>>,
        stereo_atoms: Arc<
            FixedVarBirelationSet<NodeId, Ordered, 1, StereoLigand, Ordered, StereoAtomAst>,
        >,
        stereo_bonds: Arc<
            FixedVarBirelationSet<EdgeId, Ordered, 1, StereoLigand, Ordered, StereoBondAst>,
        >,
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
            stereo_atoms,
            stereo_bonds,
            constraints,
            rings_cache: OnceLock::new(),
        }
    }

    /// AtomId/BondId-typed adapter exposing the pure-graph algorithms.
    pub fn graph(&self) -> GraphView<'_> {
        GraphView::new(&self.graph)
    }

    /// Raw underlying graph with `NodeId` / `EdgeId` types. Escape hatch
    /// for code that needs the graph-core API directly; use [`Self::graph`]
    /// for AtomId/BondId-typed access.
    #[inline]
    pub fn raw_graph(&self) -> &Graph {
        &self.graph
    }

    /// Neighbors of `atom`, ordered by ascending neighbor atom id.
    pub fn neighbors(&self, atom: AtomId) -> impl Iterator<Item = NeighborView<'_>> {
        self.graph
            .neighbors(NodeId::from(atom))
            .iter()
            .map(move |n| NeighborView::new(AtomId::from(n.node), BondId::from(n.edge), self))
    }

    // Slices have to store node/edge ids, otherwise need owned vectors.
    pub fn atoms(&self) -> AtomViews<'_> {
        AtomViews::new(self, &self.atoms)
    }

    /// View of the atom with `id`.
    ///
    /// Panics if `id` is not an atom in this molecule. Use
    /// `self.atoms().get(id)` for checked lookup.
    pub fn atom(&self, id: AtomId) -> AtomView<'_> {
        self.atoms()
            .get(id)
            .expect("atom id must refer to an atom in this molecule")
    }

    pub fn bonds(&self) -> BondViews<'_> {
        BondViews::new(self, &self.bonds)
    }

    /// View of the bond with `id`.
    ///
    /// Panics if `id` is not a bond in this molecule. Use
    /// `self.bonds().get(id)` for checked lookup.
    pub fn bond(&self, id: BondId) -> BondView<'_> {
        self.bonds()
            .get(id)
            .expect("bond id must refer to a bond in this molecule")
    }

    pub fn dative_bonds(&self) -> DativeBondViews<'_> {
        DativeBondViews::new(self, &self.dative_bonds)
    }

    /// View of the dative bond with `id`.
    ///
    /// Panics if `id` is not a dative bond in this molecule. Use
    /// `self.dative_bonds().get(id)` for checked lookup.
    pub fn dative_bond(&self, id: DativeBondId) -> DativeBondView<'_> {
        self.dative_bonds()
            .get(id)
            .expect("dative bond id must refer to a dative bond in this molecule")
    }

    pub fn aromatic_systems(&self) -> AromaticSystemViews<'_> {
        AromaticSystemViews::new(self, &self.aromatic_systems)
    }

    /// View of the aromatic system with `id`.
    ///
    /// Panics if `id` is not an aromatic system in this molecule. Use
    /// `self.aromatic_systems().get(id)` for checked lookup.
    pub fn aromatic_system(&self, id: AromaticSystemId) -> AromaticSystemView<'_> {
        self.aromatic_systems()
            .get(id)
            .expect("aromatic system id must refer to an aromatic system in this molecule")
    }

    pub fn multicenter_bonds(&self) -> MulticenterBondViews<'_> {
        MulticenterBondViews::new(self, &self.multicenter_bonds)
    }

    /// View of the multicenter bond with `id`.
    ///
    /// Panics if `id` is not a multicenter bond in this molecule. Use
    /// `self.multicenter_bonds().get(id)` for checked lookup.
    pub fn multicenter_bond(&self, id: MulticenterBondId) -> MulticenterBondView<'_> {
        self.multicenter_bonds()
            .get(id)
            .expect("multicenter bond id must refer to a multicenter bond in this molecule")
    }

    pub fn noncovalent_bonds(&self) -> NoncovalentBondViews<'_> {
        NoncovalentBondViews::new(self, &self.noncovalent_bonds)
    }

    /// View of the noncovalent bond with `id`.
    ///
    /// Panics if `id` is not a noncovalent bond in this molecule. Use
    /// `self.noncovalent_bonds().get(id)` for checked lookup.
    pub fn noncovalent_bond(&self, id: NoncovalentBondId) -> NoncovalentBondView<'_> {
        self.noncovalent_bonds()
            .get(id)
            .expect("noncovalent bond id must refer to a noncovalent bond in this molecule")
    }

    pub fn stereo_atoms(&self) -> StereoAtomViews<'_> {
        StereoAtomViews::new(self, &self.stereo_atoms)
    }

    /// View of the stereo atom with `id`.
    ///
    /// Panics if `id` is not a stereo atom in this molecule. Use
    /// `self.stereo_atoms().get(id)` for checked lookup.
    pub fn stereo_atom(&self, id: StereoAtomId) -> StereoAtomView<'_> {
        self.stereo_atoms()
            .get(id)
            .expect("stereo atom id must refer to a stereo atom in this molecule")
    }

    pub fn stereo_bonds(&self) -> StereoBondViews<'_> {
        StereoBondViews::new(self, &self.stereo_bonds)
    }

    /// View of the stereo bond with `id`.
    ///
    /// Panics if `id` is not a stereo bond in this molecule. Use
    /// `self.stereo_bonds().get(id)` for checked lookup.
    pub fn stereo_bond(&self, id: StereoBondId) -> StereoBondView<'_> {
        self.stereo_bonds()
            .get(id)
            .expect("stereo bond id must refer to a stereo bond in this molecule")
    }

    /// The subgraph induced by `atoms` (deduplicated, host order preserved), as an injective
    /// sub→host [`MoleculeCorrespondence`]: sub entity `i` maps to its host id. `extract` / `edits`
    /// materialize it. A bond/overlay is included iff all its participants are in `atoms`.
    pub fn induced_subgraph(&self, atoms: &[AtomId]) -> MoleculeCorrespondence {
        let mut host_atoms: Vec<AtomId> = Vec::with_capacity(atoms.len());
        let mut atom_set: HashSet<AtomId> = HashSet::with_capacity(atoms.len());
        for &a in atoms {
            if atom_set.insert(a) {
                host_atoms.push(a);
            }
        }

        let host_bonds: Vec<BondId> = self
            .bonds()
            .iter()
            .filter(|b| {
                let [a, b_end] = b.atom_ids();
                atom_set.contains(&a) && atom_set.contains(&b_end)
            })
            .map(|b| b.id)
            .collect();
        let host_dative_bonds: Vec<DativeBondId> = self
            .dative_bonds()
            .iter()
            .filter(|v| v.atom_ids().all(|a| atom_set.contains(&a)))
            .map(|v| v.id)
            .collect();
        let host_aromatic_systems: Vec<AromaticSystemId> = self
            .aromatic_systems()
            .iter()
            .filter(|v| v.atom_ids().all(|a| atom_set.contains(&a)))
            .map(|v| v.id)
            .collect();
        let host_multicenter_bonds: Vec<MulticenterBondId> = self
            .multicenter_bonds()
            .iter()
            .filter(|v| v.atom_ids().all(|a| atom_set.contains(&a)))
            .map(|v| v.id)
            .collect();
        let host_noncovalent_bonds: Vec<NoncovalentBondId> = self
            .noncovalent_bonds()
            .iter()
            .filter(|v| {
                let [a, b] = v.atom_ids();
                atom_set.contains(&a) && atom_set.contains(&b)
            })
            .map(|v| v.id)
            .collect();
        let host_stereo_atoms: Vec<StereoAtomId> = self
            .stereo_atoms()
            .iter()
            .filter(|v| v.atom_ids().all(|a| atom_set.contains(&a)))
            .map(|v| v.id)
            .collect();
        let host_stereo_bonds: Vec<StereoBondId> = self
            .stereo_bonds()
            .iter()
            .filter(|v| v.atom_ids().all(|a| atom_set.contains(&a)))
            .map(|v| v.id)
            .collect();

        let atom_images: Vec<NodeId> = host_atoms.iter().map(|&a| NodeId::from(a)).collect();
        MoleculeCorrespondence::new(
            Correspondence::from_images(&atom_images, self.atoms().count()),
            Correspondence::from_images(&host_bonds, self.bonds().count()),
            Correspondence::from_images(&host_dative_bonds, self.dative_bonds().count()),
            Correspondence::from_images(&host_aromatic_systems, self.aromatic_systems().count()),
            Correspondence::from_images(&host_multicenter_bonds, self.multicenter_bonds().count()),
            Correspondence::from_images(&host_noncovalent_bonds, self.noncovalent_bonds().count()),
            Correspondence::from_images(&host_stereo_atoms, self.stereo_atoms().count()),
            Correspondence::from_images(&host_stereo_bonds, self.stereo_bonds().count()),
        )
    }

    /// Materialize an induced-subgraph correspondence `sub` (over `self` as host) as a standalone
    /// molecule: drop every host atom/bond absent from `sub`. Host order preserved, gaps compacted;
    /// overlay drops cascade through the builder.
    pub fn extract(&self, sub: &MoleculeCorrespondence) -> MoleculeAst {
        let kept: HashSet<AtomId> = sub
            .atoms()
            .mates()
            .iter()
            .map(|&(_, host)| AtomId::from(host))
            .collect();
        let remove_atoms: Vec<AtomId> = (0..self.atoms().count())
            .map(AtomId::from)
            .filter(|a| !kept.contains(a))
            .collect();
        let remove_bonds: Vec<BondId> = self
            .bonds()
            .iter()
            .filter(|b| {
                let [a, b_end] = b.atom_ids();
                !kept.contains(&a) || !kept.contains(&b_end)
            })
            .map(|b| b.id)
            .collect();
        let mut builder = self.edit();
        builder.remove(&remove_atoms, &remove_bonds);
        builder.build()
    }

    /// Edits transforming `self` into the extracted subgraph `sub`: one `RemoveTopology` over the
    /// host atoms/bonds not in `sub` (empty when `sub` covers the whole molecule).
    pub fn edits(&self, sub: &MoleculeCorrespondence) -> Vec<Edit> {
        let kept: HashSet<AtomId> = sub
            .atoms()
            .mates()
            .iter()
            .map(|&(_, host)| AtomId::from(host))
            .collect();
        let kept_bonds: HashSet<BondId> =
            sub.bonds().mates().iter().map(|&(_, host)| host).collect();
        let removed_atoms: Vec<AtomHandle> = (0..self.atoms().count())
            .map(AtomId::from)
            .filter(|a| !kept.contains(a))
            .map(AtomHandle::Id)
            .collect();
        let removed_bonds: Vec<BondHandle> = (0..self.bonds().count())
            .map(BondId::from)
            .filter(|b| !kept_bonds.contains(b))
            .map(BondHandle::Id)
            .collect();
        if removed_atoms.is_empty() && removed_bonds.is_empty() {
            return Vec::new();
        }
        vec![Edit::RemoveTopology {
            atoms: removed_atoms,
            bonds: removed_bonds,
        }]
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
            && self
                .stereo_atoms
                .relation_ids()
                .all(|id| self.stereo_atoms.data(id).is_ground())
            && self
                .stereo_bonds
                .relation_ids()
                .all(|id| self.stereo_bonds.data(id).is_ground())
    }

    /// Canonical ring set: Vismara relevant cycles up to max ring size 22,
    /// applied to every atom. Cached in a single-slot `OnceLock` populated
    /// lazily on first call; subsequent calls return the same borrow.
    pub fn rings(&self) -> &RingSet {
        self.rings_cache
            .get_or_init(|| RingSet::enumerate(RingFamily::Relevant, 22, |_| true, &self.graph))
    }

    /// Ring enumeration with caller-specified family, maximum size, and
    /// atom filter. Uncached; each call recomputes.
    pub fn rings_with(
        &self,
        family: RingFamily,
        max_ring_size: usize,
        atom_filter: impl Fn(AtomId) -> bool,
    ) -> RingSet {
        RingSet::enumerate(family, max_ring_size, atom_filter, &self.graph)
    }

    pub fn atom_mut(&mut self, id: AtomId) -> AtomViewMut<'_> {
        let ast = &mut Arc::make_mut(&mut self.atoms)[id.index()];
        AtomViewMut { id, ast }
    }

    pub fn atoms_mut(&mut self) -> impl Iterator<Item = &mut AtomAst> {
        Arc::make_mut(&mut self.atoms).iter_mut()
    }

    pub fn bond_mut(&mut self, id: BondId) -> BondViewMut<'_> {
        let [s, t] = self.graph.edge_endpoints(id.into());
        let data = &mut Arc::make_mut(&mut self.bonds)[id.index()];
        BondViewMut::new(id, [AtomId::from(s), AtomId::from(t)], data)
    }

    pub fn bonds_mut(&mut self) -> impl Iterator<Item = &mut BondAst> {
        Arc::make_mut(&mut self.bonds).iter_mut()
    }

    pub fn dative_bond_mut(&mut self, id: DativeBondId) -> &mut DativeBondAst {
        Arc::make_mut(&mut self.dative_bonds).data_mut(RelationId::from(id))
    }

    pub fn dative_bonds_mut(&mut self) -> impl Iterator<Item = &mut DativeBondAst> {
        Arc::make_mut(&mut self.dative_bonds).data_iter_mut()
    }

    pub fn aromatic_system_mut(&mut self, id: AromaticSystemId) -> &mut AromaticSystemAst {
        Arc::make_mut(&mut self.aromatic_systems).data_mut(RelationId::from(id))
    }

    pub fn aromatic_systems_mut(&mut self) -> impl Iterator<Item = &mut AromaticSystemAst> {
        Arc::make_mut(&mut self.aromatic_systems).data_iter_mut()
    }

    pub fn multicenter_bond_mut(&mut self, id: MulticenterBondId) -> &mut MulticenterBondAst {
        Arc::make_mut(&mut self.multicenter_bonds).data_mut(RelationId::from(id))
    }

    pub fn multicenter_bonds_mut(&mut self) -> impl Iterator<Item = &mut MulticenterBondAst> {
        Arc::make_mut(&mut self.multicenter_bonds).data_iter_mut()
    }

    pub fn noncovalent_bond_mut(&mut self, id: NoncovalentBondId) -> &mut NoncovalentBondAst {
        Arc::make_mut(&mut self.noncovalent_bonds).data_mut(RelationId::from(id))
    }

    pub fn noncovalent_bonds_mut(&mut self) -> impl Iterator<Item = &mut NoncovalentBondAst> {
        Arc::make_mut(&mut self.noncovalent_bonds).data_iter_mut()
    }

    pub fn stereo_atom_mut(&mut self, id: StereoAtomId) -> &mut StereoAtomAst {
        Arc::make_mut(&mut self.stereo_atoms).data_mut(RelationId::from(id))
    }

    pub fn stereo_atoms_mut(&mut self) -> impl Iterator<Item = &mut StereoAtomAst> {
        Arc::make_mut(&mut self.stereo_atoms).data_iter_mut()
    }

    pub fn stereo_bond_mut(&mut self, id: StereoBondId) -> &mut StereoBondAst {
        Arc::make_mut(&mut self.stereo_bonds).data_mut(RelationId::from(id))
    }

    pub fn stereo_bonds_mut(&mut self) -> impl Iterator<Item = &mut StereoBondAst> {
        Arc::make_mut(&mut self.stereo_bonds).data_iter_mut()
    }

    pub fn constraints(&self) -> &Constraints {
        &self.constraints
    }

    pub fn constraints_mut(&mut self) -> &mut Constraints {
        &mut self.constraints
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// True if the molecule-scope `Constraints` tree is non-empty.
    /// Per-entity constraint stores are not consulted.
    pub fn has_constraints(&self) -> bool {
        !self.constraints.is_empty()
    }

    pub fn has_dative_bonds(&self) -> bool {
        self.dative_bonds.relation_count() > 0
    }

    pub fn has_aromatic_systems(&self) -> bool {
        self.aromatic_systems.relation_count() > 0
    }

    pub fn has_multicenter_bonds(&self) -> bool {
        self.multicenter_bonds.relation_count() > 0
    }

    pub fn has_noncovalent_bonds(&self) -> bool {
        self.noncovalent_bonds.relation_count() > 0
    }

    pub fn has_stereo_atoms(&self) -> bool {
        self.stereo_atoms.relation_count() > 0
    }

    pub fn has_stereo_bonds(&self) -> bool {
        self.stereo_bonds.relation_count() > 0
    }

    /// Drain every entity's inline `constraints` store into `self.constraints`
    /// as `Constraint::Atom` / `Bond` / `DativeBond` / `AromaticSystem` /
    /// `MulticenterBond` / `NoncovalentBond` / `StereoAtom` / `StereoBond`
    /// entries. The order of inserted
    /// entries in `self.constraints` is unspecified.
    pub fn lift_constraints(&mut self) {
        let atom_count = self.atoms().count();
        let bond_count = self.bonds().count();
        let dative_count = self.dative_bonds().count();
        let aromatic_count = self.aromatic_systems().count();
        let multicenter_count = self.multicenter_bonds().count();
        let noncovalent_count = self.noncovalent_bonds().count();
        let stereo_atom_count = self.stereo_atoms().count();
        let stereo_bond_count = self.stereo_bonds().count();

        let mut additions: Vec<Constraint> = Vec::new();
        for i in 0..atom_count {
            let id = AtomId::from(i);
            for c in self.atom_mut(id).ast.constraints.take() {
                additions.push(Constraint::Atom(id, c));
            }
        }
        for i in 0..bond_count {
            let id = BondId::from(i);
            for c in self.bond_mut(id).ast.constraints.take() {
                additions.push(Constraint::Bond(id, c));
            }
        }
        for i in 0..dative_count {
            let id = DativeBondId::from(i);
            for c in self.dative_bond_mut(id).constraints.take() {
                additions.push(Constraint::DativeBond(id, c));
            }
        }
        for i in 0..aromatic_count {
            let id = AromaticSystemId::from(i);
            for c in self.aromatic_system_mut(id).constraints.take() {
                additions.push(Constraint::AromaticSystem(id, c));
            }
        }
        for i in 0..multicenter_count {
            let id = MulticenterBondId::from(i);
            for c in self.multicenter_bond_mut(id).constraints.take() {
                additions.push(Constraint::MulticenterBond(id, c));
            }
        }
        for i in 0..noncovalent_count {
            let id = NoncovalentBondId::from(i);
            for c in self.noncovalent_bond_mut(id).constraints.take() {
                additions.push(Constraint::NoncovalentBond(id, c));
            }
        }
        for i in 0..stereo_atom_count {
            let id = StereoAtomId::from(i);
            let kind = self
                .stereo_atom_mut(id)
                .configuration
                .kind()
                .expect("molecule stereo atom has a concrete kind");
            for c in self.stereo_atom_mut(id).constraints.take() {
                additions.push(Constraint::StereoAtom(id, kind, c));
            }
        }
        for i in 0..stereo_bond_count {
            let id = StereoBondId::from(i);
            let kind = self
                .stereo_bond_mut(id)
                .configuration
                .kind()
                .expect("molecule stereo bond has a concrete kind");
            for c in self.stereo_bond_mut(id).constraints.take() {
                additions.push(Constraint::StereoBond(id, kind, c));
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
    pub fn inline_constraints(&mut self) {
        let entries = self.constraints.take();
        let mut leftover: Vec<Constraint> = Vec::new();
        for c in entries {
            match c {
                Constraint::Atom(id, inner) => {
                    self.atom_mut(id).ast.constraints.set(inner);
                }
                Constraint::Bond(id, inner) => {
                    self.bond_mut(id).ast.constraints.set(inner);
                }
                Constraint::DativeBond(id, inner) => {
                    self.dative_bond_mut(id).constraints.set(inner);
                }
                Constraint::AromaticSystem(id, inner) => {
                    self.aromatic_system_mut(id).constraints.set(inner);
                }
                Constraint::MulticenterBond(id, inner) => {
                    self.multicenter_bond_mut(id).constraints.set(inner);
                }
                Constraint::NoncovalentBond(_, inner) => match inner {},
                // The carried kind is dropped here; kind/degree consistency
                // against the element is the C4 validator's job.
                Constraint::StereoAtom(id, _kind, inner) => {
                    self.stereo_atom_mut(id).constraints.set(inner);
                }
                Constraint::StereoBond(id, _kind, inner) => {
                    self.stereo_bond_mut(id).constraints.set(inner);
                }
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

    pub fn edit(&self) -> MoleculeEditor {
        MoleculeEditor::from_parts(
            self.graph.clone(),
            Arc::clone(&self.atoms),
            Arc::clone(&self.bonds),
            Arc::clone(&self.dative_bonds),
            Arc::clone(&self.aromatic_systems),
            Arc::clone(&self.multicenter_bonds),
            Arc::clone(&self.noncovalent_bonds),
            Arc::clone(&self.stereo_atoms),
            Arc::clone(&self.stereo_bonds),
            self.constraints.clone(),
        )
    }

    /// Append `other` disjointly, in place: `self` keeps its ids (the prefix) and `other`'s
    /// entities follow, each id offset by `self`'s per-family count. Returns the
    /// `MoleculeCorrespondence` mapping `other`'s ids to their new ids in the union. Pure
    /// renumbering — no gluing, no chemistry.
    pub fn join_from(&mut self, other: &MoleculeAst) -> MoleculeCorrespondence {
        let atom_offset = self.atoms().count();
        let bond_offset = self.bonds().count();
        let dative_offset = self.dative_bonds().count();
        let aromatic_offset = self.aromatic_systems().count();
        let multicenter_offset = self.multicenter_bonds().count();
        let noncovalent_offset = self.noncovalent_bonds().count();
        let stereo_atom_count = self.stereo_atoms().count();
        let stereo_bond_count = self.stereo_bonds().count();
        let shift = |id: AtomId| AtomId(id.0 + atom_offset as u32);

        let mut editor = self.edit();
        for atom in other.atoms().iter() {
            editor.add_atom(atom.ast.clone());
        }
        for bond in other.bonds().iter() {
            let [a, b] = bond.atom_ids();
            editor.add_bond(shift(a), shift(b), bond.ast.clone());
        }
        for dative in other.dative_bonds().iter() {
            let donors: Vec<AtomId> = dative.donors().map(|d| shift(d.id)).collect();
            editor.add_dative_bond(donors, shift(dative.acceptor_id()), dative.ast.clone());
        }
        for system in other.aromatic_systems().iter() {
            let atoms: Vec<AtomId> = system.atom_ids().map(shift).collect();
            editor.add_aromatic_system(atoms, system.ast.clone());
        }
        for bond in other.multicenter_bonds().iter() {
            let atoms: Vec<AtomId> = bond.atom_ids().map(shift).collect();
            editor.add_multicenter_bond(atoms, bond.ast.clone());
        }
        for bond in other.noncovalent_bonds().iter() {
            let [a, b] = bond.atom_ids();
            editor.add_noncovalent_bond([shift(a), shift(b)], bond.ast.clone());
        }
        let ligand_remapping = Remapping::new(
            (0..other.atoms().count())
                .map(|k| NodeId((k + atom_offset) as u32))
                .collect(),
            (0..other.bonds().count())
                .map(|k| EdgeId((k + bond_offset) as u32))
                .collect(),
        );
        for rid in other.stereo_atoms.relation_ids() {
            let site = shift(AtomId::from(other.stereo_atoms.participants_1(rid)[0]));
            let ligands: Vec<StereoLigand> = other
                .stereo_atoms
                .participants_2(rid)
                .iter()
                .map(|ligand| ligand.remap(&ligand_remapping))
                .collect();
            editor.add_stereo_atom(site, ligands, other.stereo_atoms.data(rid).clone());
        }
        for rid in other.stereo_bonds.relation_ids() {
            let site = BondId(other.stereo_bonds.participants_1(rid)[0].0 + bond_offset as u32);
            let ligands: Vec<StereoLigand> = other
                .stereo_bonds
                .participants_2(rid)
                .iter()
                .map(|ligand| ligand.remap(&ligand_remapping))
                .collect();
            editor.add_stereo_bond(site, ligands, other.stereo_bonds.data(rid).clone());
        }
        if !other.constraints.is_empty() {
            let remapping = IdRemapping::new(
                offset_map(atom_offset, other.atoms().count()),
                offset_map(bond_offset, other.bonds().count()),
                offset_map(dative_offset, other.dative_bonds().count()),
                offset_map(aromatic_offset, other.aromatic_systems().count()),
                offset_map(multicenter_offset, other.multicenter_bonds().count()),
                offset_map(noncovalent_offset, other.noncovalent_bonds().count()),
                offset_map(stereo_atom_count, other.stereo_atoms().count()),
                offset_map(stereo_bond_count, other.stereo_bonds().count()),
            );
            for constraint in other.constraints.iter() {
                editor
                    .constraints_mut()
                    .push(constraint.clone().remap(&remapping));
            }
        }
        *self = editor.build();

        MoleculeCorrespondence::new(
            offset_correspondence(atom_offset, other.atoms().count()),
            offset_correspondence(bond_offset, other.bonds().count()),
            offset_correspondence(dative_offset, other.dative_bonds().count()),
            offset_correspondence(aromatic_offset, other.aromatic_systems().count()),
            offset_correspondence(multicenter_offset, other.multicenter_bonds().count()),
            offset_correspondence(noncovalent_offset, other.noncovalent_bonds().count()),
            offset_correspondence(stereo_atom_count, other.stereo_atoms().count()),
            offset_correspondence(stereo_bond_count, other.stereo_bonds().count()),
        )
    }

    /// Disjoint union of `self` and `other` as a fresh molecule, with the correspondence mapping
    /// `other`'s ids into it. `self.clone()` then [`join_from`](Self::join_from).
    pub fn join(&self, other: &MoleculeAst) -> (MoleculeAst, MoleculeCorrespondence) {
        let mut union = self.clone();
        let correspondence = union.join_from(other);
        (union, correspondence)
    }

    /// Decompose into connected components — a conservative partition where every relation keeps its
    /// atoms in one component (a spanning overlay prevents the split rather than being dropped). Each
    /// component is a fresh, compactly-renumbered `MoleculeAst` paired with the
    /// `MoleculeCorrespondence` mapping its ids back to `self`. Components are ordered by their lowest
    /// original atom.
    pub fn split(&self) -> Vec<(MoleculeAst, MoleculeCorrespondence)> {
        let atom_count = self.atoms().count();
        let mut uf = UnionFind::new(atom_count);
        for bond in self.bonds().iter() {
            let [a, b] = bond.atom_ids();
            uf.union(a.index(), b.index());
        }
        for dative in self.dative_bonds().iter() {
            union_participants(&mut uf, dative.atom_ids());
        }
        for system in self.aromatic_systems().iter() {
            union_participants(&mut uf, system.atom_ids());
        }
        for bond in self.multicenter_bonds().iter() {
            union_participants(&mut uf, bond.atom_ids());
        }
        for bond in self.noncovalent_bonds().iter() {
            let [a, b] = bond.atom_ids();
            uf.union(a.index(), b.index());
        }
        for rid in self.stereo_atoms.relation_ids() {
            let site = AtomId::from(self.stereo_atoms.participants_1(rid)[0]);
            for ligand in self.stereo_atoms.participants_2(rid) {
                uf.union(site.index(), ligand.atom_id.index());
            }
        }
        for rid in self.stereo_bonds.relation_ids() {
            let [a, b] = self
                .bond(BondId(self.stereo_bonds.participants_1(rid)[0].0))
                .atom_ids();
            uf.union(a.index(), b.index());
            for ligand in self.stereo_bonds.participants_2(rid) {
                uf.union(a.index(), ligand.atom_id.index());
            }
        }
        for constraint in self.constraints.iter() {
            union_participants(&mut uf, self.constraint_atoms(constraint));
        }

        let mut atom_component = vec![0usize; atom_count];
        let mut atom_compact = vec![0u32; atom_count];
        let mut component_atoms: Vec<Vec<AtomId>> = Vec::new();
        let mut index_of_root: HashMap<usize, usize> = HashMap::new();
        for i in 0..atom_count {
            let root = uf.find(i);
            let component = *index_of_root.entry(root).or_insert_with(|| {
                component_atoms.push(Vec::new());
                component_atoms.len() - 1
            });
            atom_component[i] = component;
            atom_compact[i] = component_atoms[component].len() as u32;
            component_atoms[component].push(AtomId(i as u32));
        }
        let compact = |a: AtomId| AtomId(atom_compact[a.index()]);
        let component_of = |a: AtomId| atom_component[a.index()];
        let compaction = Remapping::new(
            (0..atom_count).map(|i| NodeId(atom_compact[i])).collect(),
            Vec::new(),
        );

        component_atoms
            .iter()
            .enumerate()
            .map(|(component, atoms)| {
                let mut editor = MoleculeAst::new().edit();
                for atom in atoms {
                    editor.add_atom(self[*atom].clone());
                }
                let mut bond_compact: HashMap<BondId, BondId> = HashMap::new();
                for (original, bond) in self.bonds().iter().enumerate() {
                    let [a, b] = bond.atom_ids();
                    if component_of(a) == component {
                        let new_bond = editor.add_bond(compact(a), compact(b), bond.ast.clone());
                        bond_compact.insert(BondId(original as u32), new_bond);
                    }
                }
                for dative in self.dative_bonds().iter() {
                    if component_of(dative.acceptor_id()) == component {
                        let donors = dative.donors().map(|d| compact(d.id)).collect();
                        editor.add_dative_bond(
                            donors,
                            compact(dative.acceptor_id()),
                            dative.ast.clone(),
                        );
                    }
                }
                for system in self.aromatic_systems().iter() {
                    let members: Vec<AtomId> = system.atom_ids().collect();
                    if members
                        .first()
                        .is_some_and(|a| component_of(*a) == component)
                    {
                        editor.add_aromatic_system(
                            members.iter().map(|a| compact(*a)).collect(),
                            system.ast.clone(),
                        );
                    }
                }
                for bond in self.multicenter_bonds().iter() {
                    let members: Vec<AtomId> = bond.atom_ids().collect();
                    if members
                        .first()
                        .is_some_and(|a| component_of(*a) == component)
                    {
                        editor.add_multicenter_bond(
                            members.iter().map(|a| compact(*a)).collect(),
                            bond.ast.clone(),
                        );
                    }
                }
                for bond in self.noncovalent_bonds().iter() {
                    let [a, b] = bond.atom_ids();
                    if component_of(a) == component {
                        editor.add_noncovalent_bond([compact(a), compact(b)], bond.ast.clone());
                    }
                }
                for rid in self.stereo_atoms.relation_ids() {
                    let site = AtomId::from(self.stereo_atoms.participants_1(rid)[0]);
                    if component_of(site) == component {
                        let ligands: Vec<StereoLigand> = self
                            .stereo_atoms
                            .participants_2(rid)
                            .iter()
                            .map(|ligand| ligand.remap(&compaction))
                            .collect();
                        editor.add_stereo_atom(
                            compact(site),
                            ligands,
                            self.stereo_atoms.data(rid).clone(),
                        );
                    }
                }
                for rid in self.stereo_bonds.relation_ids() {
                    let bond = BondId(self.stereo_bonds.participants_1(rid)[0].0);
                    let [a, _] = self.bond(bond).atom_ids();
                    if component_of(a) == component {
                        let ligands: Vec<StereoLigand> = self
                            .stereo_bonds
                            .participants_2(rid)
                            .iter()
                            .map(|ligand| ligand.remap(&compaction))
                            .collect();
                        editor.add_stereo_bond(
                            bond_compact[&bond],
                            ligands,
                            self.stereo_bonds.data(rid).clone(),
                        );
                    }
                }
                let entities = editor.build();

                let atom_images: Vec<NodeId> = atoms.iter().map(|a| NodeId(a.0)).collect();
                let atom_correspondence = Correspondence::from_images(&atom_images, atom_count);
                let correspondence =
                    MoleculeCorrespondence::induce(&entities, self, atom_correspondence);

                // route each constraint to the component holding its atoms (conservative union
                // guarantees they share one), remapped to compact ids
                let remapping = idremapping_from_correspondence(&correspondence);
                let mut editor = entities.edit();
                for constraint in self.constraints.iter() {
                    if self
                        .constraint_atoms(constraint)
                        .first()
                        .is_some_and(|a| component_of(*a) == component)
                    {
                        editor
                            .constraints_mut()
                            .push(constraint.clone().remap(&remapping));
                    }
                }
                (editor.build(), correspondence)
            })
            .collect()
    }

    /// Every atom a constraint transitively references — the conservative binding set `split` unions
    /// and the connectivity validator checks (mirrors the ids each `Constraint` arm carries, resolved
    /// to atoms; a `None`-scoped molecule constraint binds the whole molecule).
    pub fn constraint_atoms(&self, constraint: &Constraint) -> Vec<AtomId> {
        let all_atoms = || {
            (0..self.atoms().count() as u32)
                .map(AtomId)
                .collect::<Vec<_>>()
        };
        match constraint {
            Constraint::Atom(id, _) => vec![*id],
            Constraint::Bond(id, _) => self.bond(*id).atom_ids().to_vec(),
            Constraint::DativeBond(id, _) => self.dative_bond(*id).atom_ids().collect(),
            Constraint::AromaticSystem(id, _) => self.aromatic_system(*id).atom_ids().collect(),
            Constraint::MulticenterBond(id, _) => self.multicenter_bond(*id).atom_ids().collect(),
            Constraint::NoncovalentBond(id, _) => self.noncovalent_bond(*id).atom_ids().to_vec(),
            Constraint::StereoAtom(id, _, _) => {
                let view = self.stereo_atom(*id);
                iter::once(view.site_id())
                    .chain(view.ligands().map(|ligand| ligand.atom_id()))
                    .collect()
            }
            Constraint::StereoBond(id, _, _) => {
                let view = self.stereo_bond(*id);
                let [a, b] = self.bond(view.site_id()).atom_ids();
                [a, b]
                    .into_iter()
                    .chain(view.ligands().map(|ligand| ligand.atom_id()))
                    .collect()
            }
            Constraint::Relational(relational) => self.relational_constraint_atoms(relational),
            Constraint::Molecule(molecule) => match molecule {
                MoleculeConstraint::ChargeSum { atoms, .. }
                | MoleculeConstraint::SpinSum { atoms, .. }
                | MoleculeConstraint::Connected { atoms } => {
                    atoms.clone().unwrap_or_else(all_atoms)
                }
                MoleculeConstraint::BondOrderSum { bonds, .. } => match bonds {
                    Some(bonds) => bonds
                        .iter()
                        .flat_map(|b| self.bond(*b).atom_ids())
                        .collect(),
                    None => all_atoms(),
                },
                MoleculeConstraint::SubPattern { anchor, .. } => {
                    anchor.atoms().iter().flat_map(|&(a, b)| [a, b]).collect()
                }
            },
            Constraint::And(constraints) | Constraint::Or(constraints) => constraints
                .iter()
                .flat_map(|c| self.constraint_atoms(c))
                .collect(),
            Constraint::Not(constraint) => self.constraint_atoms(constraint),
        }
    }

    /// A relational constraint's atoms: its primary entity's atoms plus any explicit atom operands.
    fn relational_constraint_atoms(&self, constraint: &RelationalConstraint) -> Vec<AtomId> {
        let dative = |id, extra: &[AtomId]| {
            self.dative_bond(id)
                .atom_ids()
                .chain(extra.iter().copied())
                .collect::<Vec<_>>()
        };
        let aromatic = |id, extra: &[AtomId]| {
            self.aromatic_system(id)
                .atom_ids()
                .chain(extra.iter().copied())
                .collect::<Vec<_>>()
        };
        let multicenter = |id, extra: &[AtomId]| {
            self.multicenter_bond(id)
                .atom_ids()
                .chain(extra.iter().copied())
                .collect::<Vec<_>>()
        };
        let noncovalent = |id, extra: &[AtomId]| {
            self.noncovalent_bond(id)
                .atom_ids()
                .into_iter()
                .chain(extra.iter().copied())
                .collect::<Vec<_>>()
        };
        let stereo_atom = |id, extra: &[AtomId]| {
            let view = self.stereo_atom(id);
            iter::once(view.site_id())
                .chain(view.ligands().map(|ligand| ligand.atom_id()))
                .chain(extra.iter().copied())
                .collect::<Vec<_>>()
        };
        let stereo_bond = |id, extra: &[AtomId]| {
            let view = self.stereo_bond(id);
            let [a, b] = self.bond(view.site_id()).atom_ids();
            [a, b]
                .into_iter()
                .chain(view.ligands().map(|ligand| ligand.atom_id()))
                .chain(extra.iter().copied())
                .collect::<Vec<_>>()
        };
        match constraint {
            RelationalConstraint::DativeBondDonors { bond, atoms } => dative(*bond, atoms),
            RelationalConstraint::DativeBondDonor { bond, atom } => dative(*bond, &[*atom]),
            RelationalConstraint::DativeBondContainsAllDonors { bond, atoms } => {
                dative(*bond, atoms)
            }
            RelationalConstraint::DativeBondAllDonors { bond, .. } => dative(*bond, &[]),
            RelationalConstraint::DativeBondAnyDonor { bond, .. } => dative(*bond, &[]),
            RelationalConstraint::DativeBondAcceptor { bond, atom } => dative(*bond, &[*atom]),
            RelationalConstraint::DativeBondAcceptorSatisfies { bond, .. } => dative(*bond, &[]),
            RelationalConstraint::DativeBondParallels {
                dative: id,
                parallel,
            } => dative(*id, &self.bond(*parallel).atom_ids()),
            RelationalConstraint::AromaticSystemAtoms { system, atoms } => aromatic(*system, atoms),
            RelationalConstraint::AromaticSystemContains { system, atom } => {
                aromatic(*system, &[*atom])
            }
            RelationalConstraint::AromaticSystemContainsAll { system, atoms } => {
                aromatic(*system, atoms)
            }
            RelationalConstraint::AromaticSystemAllAtoms { system, .. } => aromatic(*system, &[]),
            RelationalConstraint::AromaticSystemAnyAtom { system, .. } => aromatic(*system, &[]),
            RelationalConstraint::MulticenterBondAtoms { bond, atoms } => multicenter(*bond, atoms),
            RelationalConstraint::MulticenterBondContains { bond, atom } => {
                multicenter(*bond, &[*atom])
            }
            RelationalConstraint::MulticenterBondContainsAll { bond, atoms } => {
                multicenter(*bond, atoms)
            }
            RelationalConstraint::MulticenterBondAllAtoms { bond, .. } => multicenter(*bond, &[]),
            RelationalConstraint::MulticenterBondAnyAtom { bond, .. } => multicenter(*bond, &[]),
            RelationalConstraint::NoncovalentBondEnds { bond, atoms } => noncovalent(*bond, atoms),
            RelationalConstraint::NoncovalentBondContains { bond, atom } => {
                noncovalent(*bond, &[*atom])
            }
            RelationalConstraint::NoncovalentBondEndsSatisfy { bond, .. } => {
                noncovalent(*bond, &[])
            }
            RelationalConstraint::StereoAtomSite {
                stereo_atom: id,
                atom,
            } => stereo_atom(*id, &[*atom]),
            RelationalConstraint::StereoAtomContains {
                stereo_atom: id,
                atom,
            } => stereo_atom(*id, &[*atom]),
            RelationalConstraint::StereoAtomLigands {
                stereo_atom: id,
                atoms,
            } => stereo_atom(*id, atoms),
            RelationalConstraint::StereoAtomAllLigands {
                stereo_atom: id, ..
            } => stereo_atom(*id, &[]),
            RelationalConstraint::StereoAtomAnyLigand {
                stereo_atom: id, ..
            } => stereo_atom(*id, &[]),
            RelationalConstraint::StereoBondSite {
                stereo_bond: id,
                bond,
            } => stereo_bond(*id, &self.bond(*bond).atom_ids()),
            RelationalConstraint::StereoBondContains {
                stereo_bond: id,
                atom,
            } => stereo_bond(*id, &[*atom]),
            RelationalConstraint::StereoBondLigands {
                stereo_bond: id,
                atoms,
            } => stereo_bond(*id, atoms),
            RelationalConstraint::StereoBondAllLigands {
                stereo_bond: id, ..
            } => stereo_bond(*id, &[]),
            RelationalConstraint::StereoBondAnyLigand {
                stereo_bond: id, ..
            } => stereo_bond(*id, &[]),
        }
    }
}

/// The correspondence mapping a family's `other` ids (`0..other_count`) to their union ids
/// (`offset..offset+other_count`) — the per-family offset a disjoint `join` induces.
fn offset_correspondence<Id: Copy + Ord + From<usize>>(
    offset: usize,
    other_count: usize,
) -> Correspondence<Id> {
    let images: Vec<Id> = (0..other_count).map(|k| Id::from(offset + k)).collect();
    Correspondence::from_images(&images, offset + other_count)
}

/// The per-family offset map (`other` ids → their union ids) that `join_from` feeds to
/// `Constraints::remap`.
fn offset_map<Id: Copy + Eq + Hash + From<usize>>(offset: usize, count: usize) -> HashMap<Id, Id> {
    (0..count)
        .map(|k| (Id::from(k), Id::from(offset + k)))
        .collect()
}

/// Union every atom of a relation into one component (all participants share the first's set).
fn union_participants(uf: &mut UnionFind, atoms: impl IntoIterator<Item = AtomId>) {
    let mut atoms = atoms.into_iter();
    if let Some(first) = atoms.next() {
        for atom in atoms {
            uf.union(first.index(), atom.index());
        }
    }
}

/// The per-family `original → compact` remapping a `split` component induces, read off its
/// `component → original` correspondence.
fn idremapping_from_correspondence(correspondence: &MoleculeCorrespondence) -> IdRemapping {
    IdRemapping::new(
        correspondence
            .atoms()
            .mates()
            .iter()
            .map(|&(compact, original)| (AtomId::from(original), AtomId::from(compact)))
            .collect(),
        correspondence
            .bonds()
            .mates()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .dative_bonds()
            .mates()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .aromatic_systems()
            .mates()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .multicenter_bonds()
            .mates()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .noncovalent_bonds()
            .mates()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .stereo_atoms()
            .mates()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
        correspondence
            .stereo_bonds()
            .mates()
            .iter()
            .map(|&(compact, original)| (original, compact))
            .collect(),
    )
}

impl Index<AtomId> for MoleculeAst {
    type Output = AtomAst;
    fn index(&self, id: AtomId) -> &AtomAst {
        &self.atoms[id.index()]
    }
}

impl Index<BondId> for MoleculeAst {
    type Output = BondAst;
    fn index(&self, id: BondId) -> &BondAst {
        &self.bonds[id.index()]
    }
}

impl Index<DativeBondId> for MoleculeAst {
    type Output = DativeBondAst;
    fn index(&self, id: DativeBondId) -> &DativeBondAst {
        self.dative_bonds.data(RelationId::from(id))
    }
}

impl Index<AromaticSystemId> for MoleculeAst {
    type Output = AromaticSystemAst;
    fn index(&self, id: AromaticSystemId) -> &AromaticSystemAst {
        self.aromatic_systems.data(RelationId::from(id))
    }
}

impl Index<MulticenterBondId> for MoleculeAst {
    type Output = MulticenterBondAst;
    fn index(&self, id: MulticenterBondId) -> &MulticenterBondAst {
        self.multicenter_bonds.data(RelationId::from(id))
    }
}

impl Index<NoncovalentBondId> for MoleculeAst {
    type Output = NoncovalentBondAst;
    fn index(&self, id: NoncovalentBondId) -> &NoncovalentBondAst {
        self.noncovalent_bonds.data(RelationId::from(id))
    }
}

impl Index<StereoAtomId> for MoleculeAst {
    type Output = StereoAtomAst;
    fn index(&self, id: StereoAtomId) -> &StereoAtomAst {
        self.stereo_atoms.data(RelationId::from(id))
    }
}

impl Index<StereoBondId> for MoleculeAst {
    type Output = StereoBondAst;
    fn index(&self, id: StereoBondId) -> &StereoBondAst {
        self.stereo_bonds.data(RelationId::from(id))
    }
}

#[cfg(test)]
mod tests;
