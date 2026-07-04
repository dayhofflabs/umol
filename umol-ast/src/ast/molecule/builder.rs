//! Structural editing for `MoleculeAst`. The AST itself only allows attribute
//! mutation; structural change (add atoms/bonds/relations, remove anything)
//! goes through `MoleculeBuilder`.
//!
//! Storage is lazy: each Arc-shared field stays shared until first write,
//! at which point only that field decomposes to a mutable form. `build`
//! re-wraps everything in `Arc`, reusing untouched shared data.

use std::collections::HashSet;
use std::mem;
use std::sync::Arc;

use umol_graph_core::{
    compact_edge_vec, compact_node_vec, Compaction, EdgeId, FactorOrdering, FixedRelationSet,
    FixedVarBirelationSet, Graph, NodeId, Ordered, RelationId, RelationParticipant, Unordered,
    VarRelationSet,
};

use super::super::aromatic::AromaticSystemAst;
use super::super::atom::AtomAst;
use super::super::bond::BondAst;
use super::super::constraint::{Constraint, Constraints};
use super::super::dative::DativeBondAst;
use super::super::edit::{
    AddedAromaticSystem, AddedAtom, AddedBond, AddedDativeBond, AddedMulticenterBond,
    AddedNoncovalentBond, AddedStereoAtom, AddedStereoBond, CascadedConstraints,
    RemovedAromaticSystem, RemovedAtom, RemovedBond, RemovedDativeBond, RemovedMulticenterBond,
    RemovedNoncovalentBond, RemovedOverlays, RemovedStereoAtom, RemovedStereoBond,
};
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::ligand::StereoLigand;
use super::super::multicenter::MulticenterBondAst;
use super::super::noncovalent::NoncovalentBondAst;
use super::super::remap::{IdCompaction, UndoCompaction};
use super::super::stereo::{StereoAtomAst, StereoBondAst};
use super::super::view::{
    AromaticSystemBuilderView, AromaticSystemBuilderViewMut, AtomBuilderView, AtomBuilderViewMut,
    BondBuilderView, BondBuilderViewMut, DativeBondBuilderView, DativeBondBuilderViewMut,
    MulticenterBondBuilderView, MulticenterBondBuilderViewMut, NoncovalentBondBuilderView,
    NoncovalentBondBuilderViewMut, StereoAtomBuilderView, StereoAtomBuilderViewMut,
    StereoBondBuilderView, StereoBondBuilderViewMut,
};
use super::MoleculeAst;

#[derive(Clone)]
enum FixedSetStorage<P, O, D, const N: usize> {
    Shared(Arc<FixedRelationSet<P, O, D, N>>),
    Mutable(Vec<([P; N], D)>),
}

impl<P, O, D, const N: usize> FixedSetStorage<P, O, D, N>
where
    P: RelationParticipant,
    O: FactorOrdering,
    D: Clone,
{
    fn push(&mut self, participants: [P; N], data: D) -> u32 {
        self.materialize();
        let FixedSetStorage::Mutable(vec) = self else {
            unreachable!()
        };
        let id = vec.len() as u32;
        vec.push((participants, data));
        id
    }

    fn materialize(&mut self) {
        if let FixedSetStorage::Shared(arc) = self {
            let entries: Vec<([P; N], D)> = (0..arc.relation_count())
                .map(|i| {
                    let rid = RelationId(i as u32);
                    (*arc.participants(rid), arc.data(rid).clone())
                })
                .collect();
            *self = FixedSetStorage::Mutable(entries);
        }
    }

    fn into_arc(self) -> Arc<FixedRelationSet<P, O, D, N>> {
        match self {
            FixedSetStorage::Shared(arc) => arc,
            FixedSetStorage::Mutable(vec) => Arc::new(FixedRelationSet::new(vec)),
        }
    }

    fn relation_count(&self) -> usize {
        match self {
            FixedSetStorage::Shared(arc) => arc.relation_count(),
            FixedSetStorage::Mutable(vec) => vec.len(),
        }
    }

    fn participants(&self, i: usize) -> [P; N] {
        match self {
            FixedSetStorage::Shared(arc) => *arc.participants(RelationId(i as u32)),
            FixedSetStorage::Mutable(vec) => vec[i].0,
        }
    }

    fn apply_compaction(self, compaction: &Compaction) -> Self {
        match self {
            FixedSetStorage::Shared(arc) => {
                FixedSetStorage::Shared(Arc::new(arc.apply_compaction(compaction)))
            }
            FixedSetStorage::Mutable(vec) => {
                let compacted: Vec<([P; N], D)> = vec
                    .into_iter()
                    .filter_map(|(mut participants, d)| {
                        for slot in &mut participants {
                            *slot = (*slot).compact(compaction)?;
                        }
                        Some((participants, d))
                    })
                    .collect();
                FixedSetStorage::Mutable(compacted)
            }
        }
    }

    fn remove_relations(&mut self, ids: &[RelationId]) {
        if ids.is_empty() {
            return;
        }
        self.materialize();
        let FixedSetStorage::Mutable(vec) = self else {
            unreachable!()
        };
        let remove: HashSet<RelationId> = ids.iter().copied().collect();
        let mut dst = 0usize;
        for src in 0..vec.len() {
            if !remove.contains(&RelationId(src as u32)) {
                vec.swap(dst, src);
                dst += 1;
            }
        }
        vec.truncate(dst);
    }

    fn entries(&self) -> Vec<([P; N], D)> {
        match self {
            FixedSetStorage::Shared(arc) => (0..arc.relation_count())
                .map(|i| {
                    let rid = RelationId(i as u32);
                    (*arc.participants(rid), arc.data(rid).clone())
                })
                .collect(),
            FixedSetStorage::Mutable(vec) => vec.clone(),
        }
    }
}

#[derive(Clone)]
enum VarSetStorage<P, O, D> {
    Shared(Arc<VarRelationSet<P, O, D>>),
    Mutable(Vec<(Vec<P>, D)>),
}

impl<P, O, D> VarSetStorage<P, O, D>
where
    P: RelationParticipant,
    O: FactorOrdering,
    D: Clone,
{
    fn push(&mut self, participants: Vec<P>, data: D) -> u32 {
        self.materialize();
        let VarSetStorage::Mutable(vec) = self else {
            unreachable!()
        };
        let id = vec.len() as u32;
        vec.push((participants, data));
        id
    }

    fn materialize(&mut self) {
        if let VarSetStorage::Shared(arc) = self {
            let entries: Vec<(Vec<P>, D)> = (0..arc.relation_count())
                .map(|i| {
                    let rid = RelationId(i as u32);
                    (arc.participants(rid).to_vec(), arc.data(rid).clone())
                })
                .collect();
            *self = VarSetStorage::Mutable(entries);
        }
    }

    fn into_arc(self) -> Arc<VarRelationSet<P, O, D>> {
        match self {
            VarSetStorage::Shared(arc) => arc,
            VarSetStorage::Mutable(vec) => Arc::new(VarRelationSet::new(vec)),
        }
    }

    fn relation_count(&self) -> usize {
        match self {
            VarSetStorage::Shared(arc) => arc.relation_count(),
            VarSetStorage::Mutable(vec) => vec.len(),
        }
    }

    fn participants(&self, i: usize) -> Vec<P> {
        match self {
            VarSetStorage::Shared(arc) => arc.participants(RelationId(i as u32)).to_vec(),
            VarSetStorage::Mutable(vec) => vec[i].0.clone(),
        }
    }

    fn apply_compaction(self, compaction: &Compaction) -> Self {
        match self {
            VarSetStorage::Shared(arc) => {
                VarSetStorage::Shared(Arc::new(arc.apply_compaction(compaction)))
            }
            VarSetStorage::Mutable(vec) => {
                let compacted: Vec<(Vec<P>, D)> = vec
                    .into_iter()
                    .filter_map(|(participants, d)| {
                        let mapped: Option<Vec<P>> = participants
                            .into_iter()
                            .map(|p| p.compact(compaction))
                            .collect();
                        mapped.map(|p| (p, d))
                    })
                    .collect();
                VarSetStorage::Mutable(compacted)
            }
        }
    }

    fn remove_relations(&mut self, ids: &[RelationId]) {
        if ids.is_empty() {
            return;
        }
        self.materialize();
        let VarSetStorage::Mutable(vec) = self else {
            unreachable!()
        };
        let remove: HashSet<RelationId> = ids.iter().copied().collect();
        let mut dst = 0usize;
        for src in 0..vec.len() {
            if !remove.contains(&RelationId(src as u32)) {
                vec.swap(dst, src);
                dst += 1;
            }
        }
        vec.truncate(dst);
    }

    fn entries(&self) -> Vec<(Vec<P>, D)> {
        match self {
            VarSetStorage::Shared(arc) => (0..arc.relation_count())
                .map(|i| {
                    let rid = RelationId(i as u32);
                    (arc.participants(rid).to_vec(), arc.data(rid).clone())
                })
                .collect(),
            VarSetStorage::Mutable(vec) => vec.clone(),
        }
    }
}

/// Builder storage for a fixed-arity-first-factor / var-second-factor birelation:
/// shared until first mutation, then a `Vec` of entries.
#[derive(Clone)]
enum FixedVarSetStorage<L1, O1, const N1: usize, L2, O2, D> {
    Shared(Arc<FixedVarBirelationSet<L1, O1, N1, L2, O2, D>>),
    Mutable(Vec<([L1; N1], Vec<L2>, D)>),
}

impl<L1, O1, const N1: usize, L2, O2, D> FixedVarSetStorage<L1, O1, N1, L2, O2, D>
where
    L1: RelationParticipant,
    O1: FactorOrdering,
    L2: RelationParticipant,
    O2: FactorOrdering,
    D: Clone,
{
    fn push(&mut self, participants_1: [L1; N1], participants_2: Vec<L2>, data: D) -> u32 {
        self.materialize();
        let FixedVarSetStorage::Mutable(vec) = self else {
            unreachable!()
        };
        let id = vec.len() as u32;
        vec.push((participants_1, participants_2, data));
        id
    }

    fn materialize(&mut self) {
        if let FixedVarSetStorage::Shared(arc) = self {
            *self = FixedVarSetStorage::Mutable(arc_entries(arc));
        }
    }

    fn into_arc(self) -> Arc<FixedVarBirelationSet<L1, O1, N1, L2, O2, D>> {
        match self {
            FixedVarSetStorage::Shared(arc) => arc,
            FixedVarSetStorage::Mutable(vec) => Arc::new(FixedVarBirelationSet::new(vec)),
        }
    }

    fn relation_count(&self) -> usize {
        match self {
            FixedVarSetStorage::Shared(arc) => arc.relation_count(),
            FixedVarSetStorage::Mutable(vec) => vec.len(),
        }
    }

    fn participants_1(&self, i: usize) -> [L1; N1] {
        match self {
            FixedVarSetStorage::Shared(arc) => *arc.participants_1(RelationId(i as u32)),
            FixedVarSetStorage::Mutable(vec) => vec[i].0,
        }
    }

    fn participants_2(&self, i: usize) -> Vec<L2> {
        match self {
            FixedVarSetStorage::Shared(arc) => arc.participants_2(RelationId(i as u32)).to_vec(),
            FixedVarSetStorage::Mutable(vec) => vec[i].1.clone(),
        }
    }

    fn apply_compaction(self, compaction: &Compaction) -> Self {
        match self {
            FixedVarSetStorage::Shared(arc) => {
                FixedVarSetStorage::Shared(Arc::new(arc.apply_compaction(compaction)))
            }
            FixedVarSetStorage::Mutable(vec) => {
                let compacted: Vec<([L1; N1], Vec<L2>, D)> = vec
                    .into_iter()
                    .filter_map(|(mut participants_1, participants_2, d)| {
                        for slot in &mut participants_1 {
                            *slot = (*slot).compact(compaction)?;
                        }
                        let participants_2: Option<Vec<L2>> = participants_2
                            .into_iter()
                            .map(|p| p.compact(compaction))
                            .collect();
                        Some((participants_1, participants_2?, d))
                    })
                    .collect();
                FixedVarSetStorage::Mutable(compacted)
            }
        }
    }

    fn remove_relations(&mut self, ids: &[RelationId]) {
        if ids.is_empty() {
            return;
        }
        self.materialize();
        let FixedVarSetStorage::Mutable(vec) = self else {
            unreachable!()
        };
        let remove: HashSet<RelationId> = ids.iter().copied().collect();
        let mut dst = 0usize;
        for src in 0..vec.len() {
            if !remove.contains(&RelationId(src as u32)) {
                vec.swap(dst, src);
                dst += 1;
            }
        }
        vec.truncate(dst);
    }

    fn entries(&self) -> Vec<([L1; N1], Vec<L2>, D)> {
        match self {
            FixedVarSetStorage::Shared(arc) => arc_entries(arc),
            FixedVarSetStorage::Mutable(vec) => vec.clone(),
        }
    }
}

fn arc_entries<L1, O1, const N1: usize, L2, O2, D>(
    arc: &FixedVarBirelationSet<L1, O1, N1, L2, O2, D>,
) -> Vec<([L1; N1], Vec<L2>, D)>
where
    L1: RelationParticipant,
    O1: FactorOrdering,
    L2: RelationParticipant,
    O2: FactorOrdering,
    D: Clone,
{
    (0..arc.relation_count())
        .map(|i| {
            let rid = RelationId(i as u32);
            (
                *arc.participants_1(rid),
                arc.participants_2(rid).to_vec(),
                arc.data(rid).clone(),
            )
        })
        .collect()
}

/// Indices of birelations whose first or second factor maps to `None` under
/// `compaction` (i.e. dropped by the structural removal).
fn birelation_removed<L1, O1, const N1: usize, L2, O2, D>(
    storage: &FixedVarSetStorage<L1, O1, N1, L2, O2, D>,
    compaction: &Compaction,
) -> Vec<RelationId>
where
    L1: RelationParticipant,
    O1: FactorOrdering,
    L2: RelationParticipant,
    O2: FactorOrdering,
    D: Clone,
{
    let mut removed = Vec::new();
    for i in 0..storage.relation_count() {
        let f1_gone = storage
            .participants_1(i)
            .iter()
            .any(|&p| p.compact(compaction).is_none());
        let f2_gone = storage
            .participants_2(i)
            .iter()
            .any(|&p| p.compact(compaction).is_none());
        if f1_gone || f2_gone {
            removed.push(RelationId(i as u32));
        }
    }
    removed
}

/// Un-map a surviving birelation's factors back to the pre-removal coordinate
/// system during rollback.
fn restore_birelation_participants<L1, const N1: usize, L2>(
    participants_1: [L1; N1],
    participants_2: Vec<L2>,
    undo_compaction: &UndoCompaction,
) -> ([L1; N1], Vec<L2>)
where
    L1: RelationParticipant,
    L2: RelationParticipant,
{
    let graph = undo_compaction.forward().graph();
    (
        participants_1.map(|p| p.uncompact(graph)),
        participants_2
            .into_iter()
            .map(|p| p.uncompact(graph))
            .collect(),
    )
}

fn fixed_relation_removed<P, O, D, const N: usize>(
    storage: &FixedSetStorage<P, O, D, N>,
    compaction: &Compaction,
) -> Vec<RelationId>
where
    P: RelationParticipant,
    O: FactorOrdering,
    D: Clone,
{
    let mut removed = Vec::new();
    for i in 0..storage.relation_count() {
        if storage
            .participants(i)
            .iter()
            .any(|&p| p.compact(compaction).is_none())
        {
            removed.push(RelationId(i as u32));
        }
    }
    removed
}

fn var_relation_removed<P, O, D>(
    storage: &VarSetStorage<P, O, D>,
    compaction: &Compaction,
) -> Vec<RelationId>
where
    P: RelationParticipant,
    O: FactorOrdering,
    D: Clone,
{
    let mut removed = Vec::new();
    for i in 0..storage.relation_count() {
        if storage
            .participants(i)
            .iter()
            .any(|&p| p.compact(compaction).is_none())
        {
            removed.push(RelationId(i as u32));
        }
    }
    removed
}

fn restore_var_participants<P: RelationParticipant>(
    parts: Vec<P>,
    undo_compaction: &UndoCompaction,
) -> Vec<P> {
    let remapping = undo_compaction.forward().graph();
    parts.into_iter().map(|p| p.uncompact(remapping)).collect()
}

fn restore_fixed_participants<P: RelationParticipant, const N: usize>(
    parts: [P; N],
    undo_compaction: &UndoCompaction,
) -> [P; N] {
    let remapping = undo_compaction.forward().graph();
    parts.map(|p| p.uncompact(remapping))
}

/// Mutable builder for a `MoleculeAst`. Accumulates atoms, bonds, and
/// relations (dative, aromatic, multicenter, noncovalent), then finalizes
/// into an immutable `MoleculeAst`. Supports incremental removal with
/// index remapping via `remove`.
#[derive(Clone)]
pub struct MoleculeBuilder {
    graph: Graph,
    atoms: Arc<Vec<AtomAst>>,
    bonds: Arc<Vec<BondAst>>,
    dative_bonds: FixedVarSetStorage<NodeId, Ordered, 1, NodeId, Unordered, DativeBondAst>,
    aromatic_systems: VarSetStorage<NodeId, Unordered, AromaticSystemAst>,
    multicenter_bonds: VarSetStorage<NodeId, Unordered, MulticenterBondAst>,
    noncovalent_bonds: FixedSetStorage<NodeId, Unordered, NoncovalentBondAst, 2>,
    stereo_atoms: FixedVarSetStorage<NodeId, Ordered, 1, StereoLigand, Ordered, StereoAtomAst>,
    stereo_bonds: FixedVarSetStorage<EdgeId, Ordered, 1, StereoLigand, Ordered, StereoBondAst>,
    constraints: Constraints,
}

impl MoleculeBuilder {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn from_parts(
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
            dative_bonds: FixedVarSetStorage::Shared(dative_bonds),
            aromatic_systems: VarSetStorage::Shared(aromatic_systems),
            multicenter_bonds: VarSetStorage::Shared(multicenter_bonds),
            noncovalent_bonds: FixedSetStorage::Shared(noncovalent_bonds),
            stereo_atoms: FixedVarSetStorage::Shared(stereo_atoms),
            stereo_bonds: FixedVarSetStorage::Shared(stereo_bonds),
            constraints,
        }
    }

    /// Append an atom directly to the builder.
    ///
    /// This is a low-level, non-transactional construction primitive. Use
    /// `transact` for checked atomic edits or `transact_unchecked` for trusted
    /// generated edit batches.
    pub fn add_atom(&mut self, atom: AtomAst) -> AtomId {
        let id = self.graph.add_node();
        Arc::make_mut(&mut self.atoms).push(atom);
        AtomId::from(id)
    }

    /// Append a localized bond directly to the builder.
    ///
    /// This is a low-level, non-transactional construction primitive. It
    /// assumes `first` and `second` are valid atom ids in the current dense layout.
    pub fn add_bond(&mut self, first: AtomId, second: AtomId, bond: BondAst) -> BondId {
        let id = self
            .graph
            .add_edge(NodeId::from(first), NodeId::from(second));
        Arc::make_mut(&mut self.bonds).push(bond);
        BondId::from(id)
    }

    /// Append a dative-bond overlay directly to the builder. The acceptor is
    /// factor 1; the donors are factor 2 (sorted by the `Unordered`
    /// canonicalization).
    pub fn add_dative_bond(
        &mut self,
        donors: Vec<AtomId>,
        acceptor: AtomId,
        bond: DativeBondAst,
    ) -> DativeBondId {
        let donors: Vec<NodeId> = donors.into_iter().map(NodeId::from).collect();
        DativeBondId(
            self.dative_bonds
                .push([NodeId::from(acceptor)], donors, bond),
        )
    }

    /// Append an aromatic-system overlay directly to the builder.
    pub fn add_aromatic_system(
        &mut self,
        atoms: Vec<AtomId>,
        data: AromaticSystemAst,
    ) -> AromaticSystemId {
        let nodes: Vec<NodeId> = atoms.into_iter().map(NodeId::from).collect();
        let i = self.aromatic_systems.push(nodes, data);
        AromaticSystemId(i)
    }

    /// Append a multicenter-bond overlay directly to the builder.
    pub fn add_multicenter_bond(
        &mut self,
        atoms: Vec<AtomId>,
        data: MulticenterBondAst,
    ) -> MulticenterBondId {
        let nodes: Vec<NodeId> = atoms.into_iter().map(NodeId::from).collect();
        let i = self.multicenter_bonds.push(nodes, data);
        MulticenterBondId(i)
    }

    /// Append a noncovalent-bond overlay directly to the builder.
    pub fn add_noncovalent_bond(
        &mut self,
        ends: [AtomId; 2],
        bond: NoncovalentBondAst,
    ) -> NoncovalentBondId {
        let i = self
            .noncovalent_bonds
            .push([NodeId::from(ends[0]), NodeId::from(ends[1])], bond);
        NoncovalentBondId(i)
    }

    /// Append a stereo-atom overlay directly to the builder.
    pub fn add_stereo_atom(
        &mut self,
        site: AtomId,
        ligands: Vec<StereoLigand>,
        ast: StereoAtomAst,
    ) -> StereoAtomId {
        StereoAtomId(self.stereo_atoms.push([NodeId::from(site)], ligands, ast))
    }

    /// Append a stereo-bond overlay directly to the builder.
    pub fn add_stereo_bond(
        &mut self,
        site: BondId,
        ligands: Vec<StereoLigand>,
        ast: StereoBondAst,
    ) -> StereoBondId {
        StereoBondId(self.stereo_bonds.push([EdgeId::from(site)], ligands, ast))
    }

    /// Add a molecule-level constraint (molecule-scope predicate or
    /// combinator). Unconditional per-entity constraints belong inline on the
    /// entity — use `atom_mut(id).constraints.set(c)` etc.
    pub fn push_constraint(&mut self, c: Constraint) {
        self.constraints.push(c);
    }

    // -- Attribute access -----------------------------------------------------
    //
    // Mutable views edit entity data in place. Structural add/remove stays on
    // the builder itself because dense removal can compact many unrelated ids.

    pub fn atom(&self, id: AtomId) -> AtomBuilderView<'_> {
        AtomBuilderView {
            id,
            ast: &self.atoms[id.index()],
        }
    }

    pub fn atom_mut(&mut self, id: AtomId) -> AtomBuilderViewMut<'_> {
        let ast = &mut Arc::make_mut(&mut self.atoms)[id.index()];
        AtomBuilderViewMut { id, ast }
    }

    pub fn bond(&self, id: BondId) -> BondBuilderView<'_> {
        let endpoints = self.graph.edge_endpoints(EdgeId::from(id));
        let atoms = [AtomId::from(endpoints[0]), AtomId::from(endpoints[1])];
        BondBuilderView {
            id,
            ast: &self.bonds[id.index()],
            atoms,
        }
    }

    pub fn bond_mut(&mut self, id: BondId) -> BondBuilderViewMut<'_> {
        let endpoints = self.graph.edge_endpoints(EdgeId::from(id));
        let atoms = [AtomId::from(endpoints[0]), AtomId::from(endpoints[1])];
        let ast = &mut Arc::make_mut(&mut self.bonds)[id.index()];
        BondBuilderViewMut { id, ast, atoms }
    }

    pub fn dative_bond(&self, id: DativeBondId) -> DativeBondBuilderView<'_> {
        match &self.dative_bonds {
            FixedVarSetStorage::Shared(arc) => {
                let rid = RelationId(id.0);
                DativeBondBuilderView {
                    id,
                    ast: arc.data(rid),
                    acceptor_id: AtomId::from(arc.participants_1(rid)[0]),
                    donors: arc.participants_2(rid),
                }
            }
            FixedVarSetStorage::Mutable(vec) => {
                let entry = &vec[id.index()];
                DativeBondBuilderView {
                    id,
                    ast: &entry.2,
                    acceptor_id: AtomId::from(entry.0[0]),
                    donors: &entry.1,
                }
            }
        }
    }

    pub fn dative_bond_mut(&mut self, id: DativeBondId) -> DativeBondBuilderViewMut<'_> {
        self.dative_bonds.materialize();
        let FixedVarSetStorage::Mutable(vec) = &mut self.dative_bonds else {
            unreachable!()
        };
        let entry = &mut vec[id.index()];
        let acceptor_id = AtomId::from(entry.0[0]);
        DativeBondBuilderViewMut {
            id,
            acceptor_id,
            donors: &entry.1,
            ast: &mut entry.2,
        }
    }

    pub fn aromatic_system(&self, id: AromaticSystemId) -> AromaticSystemBuilderView<'_> {
        match &self.aromatic_systems {
            VarSetStorage::Shared(arc) => {
                let rid = RelationId(id.0);
                AromaticSystemBuilderView {
                    id,
                    ast: arc.data(rid),
                    atoms: arc.participants(rid),
                }
            }
            VarSetStorage::Mutable(vec) => {
                let entry = &vec[id.index()];
                AromaticSystemBuilderView {
                    id,
                    ast: &entry.1,
                    atoms: &entry.0,
                }
            }
        }
    }

    pub fn aromatic_system_mut(
        &mut self,
        id: AromaticSystemId,
    ) -> AromaticSystemBuilderViewMut<'_> {
        self.aromatic_systems.materialize();
        let VarSetStorage::Mutable(vec) = &mut self.aromatic_systems else {
            unreachable!()
        };
        let entry = &mut vec[id.index()];
        AromaticSystemBuilderViewMut {
            id,
            atoms: &entry.0,
            ast: &mut entry.1,
        }
    }

    pub fn multicenter_bond(&self, id: MulticenterBondId) -> MulticenterBondBuilderView<'_> {
        match &self.multicenter_bonds {
            VarSetStorage::Shared(arc) => {
                let rid = RelationId(id.0);
                MulticenterBondBuilderView {
                    id,
                    ast: arc.data(rid),
                    atoms: arc.participants(rid),
                }
            }
            VarSetStorage::Mutable(vec) => {
                let entry = &vec[id.index()];
                MulticenterBondBuilderView {
                    id,
                    ast: &entry.1,
                    atoms: &entry.0,
                }
            }
        }
    }

    pub fn multicenter_bond_mut(
        &mut self,
        id: MulticenterBondId,
    ) -> MulticenterBondBuilderViewMut<'_> {
        self.multicenter_bonds.materialize();
        let VarSetStorage::Mutable(vec) = &mut self.multicenter_bonds else {
            unreachable!()
        };
        let entry = &mut vec[id.index()];
        MulticenterBondBuilderViewMut {
            id,
            atoms: &entry.0,
            ast: &mut entry.1,
        }
    }

    pub fn noncovalent_bond(&self, id: NoncovalentBondId) -> NoncovalentBondBuilderView<'_> {
        match &self.noncovalent_bonds {
            FixedSetStorage::Shared(arc) => {
                let rid = RelationId(id.0);
                let parts = arc.participants(rid);
                NoncovalentBondBuilderView {
                    id,
                    ast: arc.data(rid),
                    atoms: [AtomId::from(parts[0]), AtomId::from(parts[1])],
                }
            }
            FixedSetStorage::Mutable(vec) => {
                let entry = &vec[id.index()];
                NoncovalentBondBuilderView {
                    id,
                    ast: &entry.1,
                    atoms: [AtomId::from(entry.0[0]), AtomId::from(entry.0[1])],
                }
            }
        }
    }

    pub fn noncovalent_bond_mut(
        &mut self,
        id: NoncovalentBondId,
    ) -> NoncovalentBondBuilderViewMut<'_> {
        self.noncovalent_bonds.materialize();
        let FixedSetStorage::Mutable(vec) = &mut self.noncovalent_bonds else {
            unreachable!()
        };
        let entry = &mut vec[id.index()];
        let atoms = [AtomId::from(entry.0[0]), AtomId::from(entry.0[1])];
        NoncovalentBondBuilderViewMut {
            id,
            ast: &mut entry.1,
            atoms,
        }
    }

    pub fn stereo_atom(&self, id: StereoAtomId) -> StereoAtomBuilderView<'_> {
        match &self.stereo_atoms {
            FixedVarSetStorage::Shared(arc) => {
                let rid = RelationId(id.0);
                StereoAtomBuilderView {
                    id,
                    ast: arc.data(rid),
                    site: AtomId::from(arc.participants_1(rid)[0]),
                    ligands: arc.participants_2(rid),
                }
            }
            FixedVarSetStorage::Mutable(vec) => {
                let entry = &vec[id.index()];
                StereoAtomBuilderView {
                    id,
                    ast: &entry.2,
                    site: AtomId::from(entry.0[0]),
                    ligands: &entry.1,
                }
            }
        }
    }

    pub fn stereo_bond(&self, id: StereoBondId) -> StereoBondBuilderView<'_> {
        match &self.stereo_bonds {
            FixedVarSetStorage::Shared(arc) => {
                let rid = RelationId(id.0);
                StereoBondBuilderView {
                    id,
                    ast: arc.data(rid),
                    site: BondId::from(arc.participants_1(rid)[0]),
                    ligands: arc.participants_2(rid),
                }
            }
            FixedVarSetStorage::Mutable(vec) => {
                let entry = &vec[id.index()];
                StereoBondBuilderView {
                    id,
                    ast: &entry.2,
                    site: BondId::from(entry.0[0]),
                    ligands: &entry.1,
                }
            }
        }
    }

    pub fn stereo_atom_mut(&mut self, id: StereoAtomId) -> StereoAtomBuilderViewMut<'_> {
        self.stereo_atoms.materialize();
        let FixedVarSetStorage::Mutable(vec) = &mut self.stereo_atoms else {
            unreachable!()
        };
        let entry = &mut vec[id.index()];
        let site = AtomId::from(entry.0[0]);
        StereoAtomBuilderViewMut {
            id,
            ast: &mut entry.2,
            site,
            ligands: &entry.1,
        }
    }

    pub fn stereo_bond_mut(&mut self, id: StereoBondId) -> StereoBondBuilderViewMut<'_> {
        self.stereo_bonds.materialize();
        let FixedVarSetStorage::Mutable(vec) = &mut self.stereo_bonds else {
            unreachable!()
        };
        let entry = &mut vec[id.index()];
        let site = BondId::from(entry.0[0]);
        StereoBondBuilderViewMut {
            id,
            ast: &mut entry.2,
            site,
            ligands: &entry.1,
        }
    }

    pub fn constraints(&self) -> &Constraints {
        &self.constraints
    }

    pub fn constraints_mut(&mut self) -> &mut Constraints {
        &mut self.constraints
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    pub fn dative_bond_count(&self) -> usize {
        self.dative_bonds.relation_count()
    }

    pub fn aromatic_system_count(&self) -> usize {
        self.aromatic_systems.relation_count()
    }

    pub fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.relation_count()
    }

    pub fn noncovalent_bond_count(&self) -> usize {
        self.noncovalent_bonds.relation_count()
    }

    pub fn stereo_atom_count(&self) -> usize {
        self.stereo_atoms.relation_count()
    }

    pub fn stereo_bond_count(&self) -> usize {
        self.stereo_bonds.relation_count()
    }

    // -- Relation removal -----------------------------------------------------

    /// Remove dative-bond overlays directly from the builder.
    ///
    /// This is a low-level dense removal primitive. It compacts molecule-level
    /// constraints but does not build rollback data.
    pub fn remove_dative_bonds(&mut self, ids: &[DativeBondId]) {
        let raw: Vec<RelationId> = ids.iter().map(|&i| i.into()).collect();
        self.dative_bonds.remove_relations(&raw);
        let id_compaction = IdCompaction::new(
            Compaction::new(Vec::new(), Vec::new()),
            raw,
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        self.constraints.compact(&id_compaction);
    }

    /// Remove aromatic-system overlays directly from the builder.
    ///
    /// This is a low-level dense removal primitive. It compacts molecule-level
    /// constraints but does not build rollback data.
    pub fn remove_aromatic_systems(&mut self, ids: &[AromaticSystemId]) {
        let raw: Vec<RelationId> = ids.iter().map(|&i| i.into()).collect();
        self.aromatic_systems.remove_relations(&raw);
        let id_compaction = IdCompaction::new(
            Compaction::new(Vec::new(), Vec::new()),
            Vec::new(),
            raw,
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        self.constraints.compact(&id_compaction);
    }

    /// Remove multicenter-bond overlays directly from the builder.
    ///
    /// This is a low-level dense removal primitive. It compacts molecule-level
    /// constraints but does not build rollback data.
    pub fn remove_multicenter_bonds(&mut self, ids: &[MulticenterBondId]) {
        let raw: Vec<RelationId> = ids.iter().map(|&i| i.into()).collect();
        self.multicenter_bonds.remove_relations(&raw);
        let id_compaction = IdCompaction::new(
            Compaction::new(Vec::new(), Vec::new()),
            Vec::new(),
            Vec::new(),
            raw,
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        self.constraints.compact(&id_compaction);
    }

    /// Remove noncovalent-bond overlays directly from the builder.
    ///
    /// This is a low-level dense removal primitive. It compacts molecule-level
    /// constraints but does not build rollback data.
    pub fn remove_noncovalent_bonds(&mut self, ids: &[NoncovalentBondId]) {
        let raw: Vec<RelationId> = ids.iter().map(|&i| i.into()).collect();
        self.noncovalent_bonds.remove_relations(&raw);
        let id_compaction = IdCompaction::new(
            Compaction::new(Vec::new(), Vec::new()),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            raw,
            Vec::new(),
            Vec::new(),
        );
        self.constraints.compact(&id_compaction);
    }

    /// Remove stereo-atom overlays directly from the builder.
    ///
    /// Low-level dense removal primitive; compacts molecule-level constraints but
    /// does not build rollback data.
    pub fn remove_stereo_atoms(&mut self, ids: &[StereoAtomId]) {
        let raw: Vec<RelationId> = ids.iter().map(|&i| i.into()).collect();
        self.stereo_atoms.remove_relations(&raw);
        let id_compaction = IdCompaction::new(
            Compaction::new(Vec::new(), Vec::new()),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            raw,
            Vec::new(),
        );
        self.constraints.compact(&id_compaction);
    }

    /// Remove stereo-bond overlays directly from the builder.
    ///
    /// Low-level dense removal primitive; compacts molecule-level constraints but
    /// does not build rollback data.
    pub fn remove_stereo_bonds(&mut self, ids: &[StereoBondId]) {
        let raw: Vec<RelationId> = ids.iter().map(|&i| i.into()).collect();
        self.stereo_bonds.remove_relations(&raw);
        let id_compaction = IdCompaction::new(
            Compaction::new(Vec::new(), Vec::new()),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            raw,
        );
        self.constraints.compact(&id_compaction);
    }

    // -- Topological removal --------------------------------------------------

    /// Remove topology directly from the builder and return the forward compaction.
    ///
    /// This is the low-level dense topology-removal primitive. It removes the
    /// requested atoms and bonds, cascades relations whose participants were
    /// removed, compacts molecule-level constraints, and returns the forward
    /// `IdCompaction` for downstream id holders. It does not build rollback
    /// data; checked transactions capture the removed payloads before calling
    /// this method.
    pub fn remove(&mut self, atoms: &[AtomId], bonds: &[BondId]) -> IdCompaction {
        let nodes: Vec<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        let edges: Vec<EdgeId> = bonds.iter().map(|&b| EdgeId::from(b)).collect();
        let compaction = self.graph.remove_cascading(&nodes, &edges);

        let new_atoms = compact_node_vec(&compaction, &self.atoms);
        let new_bonds = compact_edge_vec(&compaction, &self.bonds);
        self.atoms = Arc::new(new_atoms);
        self.bonds = Arc::new(new_bonds);

        let removed_dative = birelation_removed(&self.dative_bonds, &compaction);
        let removed_aromatic = var_relation_removed(&self.aromatic_systems, &compaction);
        let removed_multicenter = var_relation_removed(&self.multicenter_bonds, &compaction);
        let removed_noncovalent = fixed_relation_removed(&self.noncovalent_bonds, &compaction);
        let removed_stereo_atoms = birelation_removed(&self.stereo_atoms, &compaction);
        let removed_stereo_bonds = birelation_removed(&self.stereo_bonds, &compaction);

        let dative = mem::replace(
            &mut self.dative_bonds,
            FixedVarSetStorage::Shared(Arc::new(FixedVarBirelationSet::default())),
        );
        self.dative_bonds = dative.apply_compaction(&compaction);

        let aromatic = mem::replace(
            &mut self.aromatic_systems,
            VarSetStorage::Shared(Arc::new(VarRelationSet::default())),
        );
        self.aromatic_systems = aromatic.apply_compaction(&compaction);

        let multicenter = mem::replace(
            &mut self.multicenter_bonds,
            VarSetStorage::Shared(Arc::new(VarRelationSet::default())),
        );
        self.multicenter_bonds = multicenter.apply_compaction(&compaction);

        let noncovalent = mem::replace(
            &mut self.noncovalent_bonds,
            FixedSetStorage::Shared(Arc::new(FixedRelationSet::default())),
        );
        self.noncovalent_bonds = noncovalent.apply_compaction(&compaction);

        // Forward-compact stereo overlays: a stereo element whose site or any
        // ligand atom/bond was removed drops out (cascade). The dropped ids
        // (computed above) feed `IdCompaction` so rollback (`restore_topology`)
        // can reinsert them.
        let stereo_atoms = mem::replace(
            &mut self.stereo_atoms,
            FixedVarSetStorage::Shared(Arc::new(FixedVarBirelationSet::default())),
        );
        self.stereo_atoms = stereo_atoms.apply_compaction(&compaction);
        let stereo_bonds = mem::replace(
            &mut self.stereo_bonds,
            FixedVarSetStorage::Shared(Arc::new(FixedVarBirelationSet::default())),
        );
        self.stereo_bonds = stereo_bonds.apply_compaction(&compaction);

        let id_compaction = IdCompaction::new(
            compaction,
            removed_dative,
            removed_aromatic,
            removed_multicenter,
            removed_noncovalent,
            removed_stereo_atoms,
            removed_stereo_bonds,
        );
        self.constraints.compact(&id_compaction);
        id_compaction
    }

    // -- Undo of additions ----------------------------------------------------

    pub(super) fn remove_added_topology(&mut self, atoms: &[AddedAtom], bonds: &[AddedBond]) {
        let atom_ids: Vec<AtomId> = atoms.iter().map(|a| a.id).collect();
        let bond_ids: Vec<BondId> = bonds.iter().map(|b| b.id).collect();
        self.remove(&atom_ids, &bond_ids);
    }

    pub(super) fn remove_added_dative_bond(&mut self, added: &AddedDativeBond) {
        self.remove_dative_bonds(&[added.id]);
    }

    pub(super) fn remove_added_aromatic_system(&mut self, added: &AddedAromaticSystem) {
        self.remove_aromatic_systems(&[added.id]);
    }

    pub(super) fn remove_added_multicenter_bond(&mut self, added: &AddedMulticenterBond) {
        self.remove_multicenter_bonds(&[added.id]);
    }

    pub(super) fn remove_added_noncovalent_bond(&mut self, added: &AddedNoncovalentBond) {
        self.remove_noncovalent_bonds(&[added.id]);
    }

    pub(super) fn remove_added_stereo_atom(&mut self, added: &AddedStereoAtom) {
        self.remove_stereo_atoms(&[added.id]);
    }

    pub(super) fn remove_added_stereo_bond(&mut self, added: &AddedStereoBond) {
        self.remove_stereo_bonds(&[added.id]);
    }

    // -- Undo of removals -----------------------------------------------------

    pub(super) fn restore_topology(
        &mut self,
        atoms: Vec<RemovedAtom>,
        bonds: Vec<RemovedBond>,
        overlays: RemovedOverlays,
        undo_compaction: &UndoCompaction,
        cascade: CascadedConstraints,
    ) {
        self.restore_atoms(atoms, undo_compaction);
        self.restore_bonds(bonds, undo_compaction);
        self.restore_dative_bonds(overlays.dative_bonds, undo_compaction);
        self.restore_aromatic_systems(overlays.aromatic_systems, undo_compaction);
        self.restore_multicenter_bonds(overlays.multicenter_bonds, undo_compaction);
        self.restore_noncovalent_bonds(overlays.noncovalent_bonds, undo_compaction);
        self.restore_stereo_atoms(overlays.stereo_atoms, undo_compaction);
        self.restore_stereo_bonds(overlays.stereo_bonds, undo_compaction);
        cascade.rollback_into(&mut self.constraints);
    }

    // -- Restore primitives ---------------------------------------------------

    fn restore_atoms(&mut self, removed: Vec<RemovedAtom>, undo_compaction: &UndoCompaction) {
        let mut next = vec![None; self.atoms.len() + removed.len()];
        for removed in removed {
            next[removed.id.index()] = Some(removed.ast);
        }
        for (idx, atom) in self.atoms.iter().cloned().enumerate() {
            let old = undo_compaction.uncompact_atom(AtomId(idx as u32));
            next[old.index()] = Some(atom);
        }
        self.atoms = Arc::new(next.into_iter().map(Option::unwrap).collect());
    }

    fn restore_bonds(&mut self, removed: Vec<RemovedBond>, undo_compaction: &UndoCompaction) {
        let mut old_endpoints: Vec<Option<[AtomId; 2]>> =
            vec![None; self.bonds.len() + removed.len()];
        let mut old_bonds: Vec<Option<BondAst>> = vec![None; self.bonds.len() + removed.len()];
        for removed in removed {
            old_endpoints[removed.id.index()] = Some(removed.endpoints);
            old_bonds[removed.id.index()] = Some(removed.ast);
        }
        for (idx, bond) in self.bonds.iter().cloned().enumerate() {
            let old_id = undo_compaction.uncompact_bond(BondId(idx as u32));
            let endpoints = self.graph.edge_endpoints(EdgeId(idx as u32));
            old_endpoints[old_id.index()] = Some([
                undo_compaction.uncompact_atom(AtomId::from(endpoints[0])),
                undo_compaction.uncompact_atom(AtomId::from(endpoints[1])),
            ]);
            old_bonds[old_id.index()] = Some(bond);
        }
        let endpoints: Vec<[u32; 2]> = old_endpoints
            .into_iter()
            .map(|e| {
                let e = e.unwrap();
                [e[0].0, e[1].0]
            })
            .collect();
        self.graph = Graph::new(self.atoms.len(), &endpoints);
        self.bonds = Arc::new(old_bonds.into_iter().map(Option::unwrap).collect());
    }

    pub(super) fn restore_dative_bonds(
        &mut self,
        removed: Vec<RemovedDativeBond>,
        undo_compaction: &UndoCompaction,
    ) {
        let current = self.dative_bonds.entries();
        let mut next = vec![None; current.len() + removed.len()];
        for (idx, (acceptor, donors, data)) in current.into_iter().enumerate() {
            let old_id = undo_compaction.uncompact_dative_bond(DativeBondId(idx as u32));
            let (acceptor, donors) =
                restore_birelation_participants(acceptor, donors, undo_compaction);
            next[old_id.index()] = Some((acceptor, donors, data));
        }
        for removed in removed {
            let (acceptor, donors) = removed
                .atoms
                .split_last()
                .expect("dative bond has an acceptor");
            next[removed.id.index()] = Some((
                [NodeId::from(*acceptor)],
                donors.iter().map(|&a| NodeId::from(a)).collect(),
                removed.ast,
            ));
        }
        self.dative_bonds =
            FixedVarSetStorage::Mutable(next.into_iter().map(Option::unwrap).collect());
    }

    pub(super) fn restore_aromatic_systems(
        &mut self,
        removed: Vec<RemovedAromaticSystem>,
        undo_compaction: &UndoCompaction,
    ) {
        let current = self.aromatic_systems.entries();
        let mut next = vec![None; current.len() + removed.len()];
        for (idx, (parts, data)) in current.into_iter().enumerate() {
            let old_id = undo_compaction.uncompact_aromatic_system(AromaticSystemId(idx as u32));
            next[old_id.index()] = Some((restore_var_participants(parts, undo_compaction), data));
        }
        for removed in removed {
            next[removed.id.index()] = Some((
                removed.atoms.into_iter().map(NodeId::from).collect(),
                removed.ast,
            ));
        }
        self.aromatic_systems =
            VarSetStorage::Mutable(next.into_iter().map(Option::unwrap).collect());
    }

    pub(super) fn restore_multicenter_bonds(
        &mut self,
        removed: Vec<RemovedMulticenterBond>,
        undo_compaction: &UndoCompaction,
    ) {
        let current = self.multicenter_bonds.entries();
        let mut next = vec![None; current.len() + removed.len()];
        for (idx, (parts, data)) in current.into_iter().enumerate() {
            let old_id = undo_compaction.uncompact_multicenter_bond(MulticenterBondId(idx as u32));
            next[old_id.index()] = Some((restore_var_participants(parts, undo_compaction), data));
        }
        for removed in removed {
            next[removed.id.index()] = Some((
                removed.atoms.into_iter().map(NodeId::from).collect(),
                removed.ast,
            ));
        }
        self.multicenter_bonds =
            VarSetStorage::Mutable(next.into_iter().map(Option::unwrap).collect());
    }

    pub(super) fn restore_noncovalent_bonds(
        &mut self,
        removed: Vec<RemovedNoncovalentBond>,
        undo_compaction: &UndoCompaction,
    ) {
        let current = self.noncovalent_bonds.entries();
        let mut next = vec![None; current.len() + removed.len()];
        for (idx, (parts, data)) in current.into_iter().enumerate() {
            let old_id = undo_compaction.uncompact_noncovalent_bond(NoncovalentBondId(idx as u32));
            next[old_id.index()] = Some((restore_fixed_participants(parts, undo_compaction), data));
        }
        for removed in removed {
            next[removed.id.index()] = Some((
                [
                    NodeId::from(removed.atoms[0]),
                    NodeId::from(removed.atoms[1]),
                ],
                removed.ast,
            ));
        }
        self.noncovalent_bonds =
            FixedSetStorage::Mutable(next.into_iter().map(Option::unwrap).collect());
    }

    pub(super) fn restore_stereo_atoms(
        &mut self,
        removed: Vec<RemovedStereoAtom>,
        undo_compaction: &UndoCompaction,
    ) {
        let current = self.stereo_atoms.entries();
        let mut next = vec![None; current.len() + removed.len()];
        for (idx, (site, ligands, data)) in current.into_iter().enumerate() {
            let old_id = undo_compaction.uncompact_stereo_atom(StereoAtomId(idx as u32));
            let (site, ligands) = restore_birelation_participants(site, ligands, undo_compaction);
            next[old_id.index()] = Some((site, ligands, data));
        }
        for removed in removed {
            next[removed.id.index()] =
                Some(([NodeId::from(removed.site)], removed.ligands, removed.ast));
        }
        self.stereo_atoms =
            FixedVarSetStorage::Mutable(next.into_iter().map(Option::unwrap).collect());
    }

    pub(super) fn restore_stereo_bonds(
        &mut self,
        removed: Vec<RemovedStereoBond>,
        undo_compaction: &UndoCompaction,
    ) {
        let current = self.stereo_bonds.entries();
        let mut next = vec![None; current.len() + removed.len()];
        for (idx, (site, ligands, data)) in current.into_iter().enumerate() {
            let old_id = undo_compaction.uncompact_stereo_bond(StereoBondId(idx as u32));
            let (site, ligands) = restore_birelation_participants(site, ligands, undo_compaction);
            next[old_id.index()] = Some((site, ligands, data));
        }
        for removed in removed {
            next[removed.id.index()] =
                Some(([EdgeId::from(removed.site)], removed.ligands, removed.ast));
        }
        self.stereo_bonds =
            FixedVarSetStorage::Mutable(next.into_iter().map(Option::unwrap).collect());
    }

    pub fn build(self) -> MoleculeAst {
        MoleculeAst::from_arcs(
            self.graph,
            self.atoms,
            self.bonds,
            self.dative_bonds.into_arc(),
            self.aromatic_systems.into_arc(),
            self.multicenter_bonds.into_arc(),
            self.noncovalent_bonds.into_arc(),
            self.stereo_atoms.into_arc(),
            self.stereo_bonds.into_arc(),
            self.constraints,
        )
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::AtomConstraint;
    use crate::ast::dative::DativeBondAst;
    use crate::ast::RemovedConstraint;
    use crate::mol;

    #[fixture]
    fn triatomic() -> MoleculeBuilder {
        let mut b = MoleculeAst::default().edit();
        b.add_atom(AtomAst::from_element(Element::C));
        b.add_atom(AtomAst::from_element(Element::N));
        b.add_atom(AtomAst::from_element(Element::O));
        b.add_bond(AtomId(0), AtomId(1), BondAst::from_order(1));
        b.add_bond(AtomId(1), AtomId(2), BondAst::from_order(2));
        b
    }

    #[rstest]
    fn test_molecule_builder_restore_topology(mut triatomic: MoleculeBuilder) {
        let dropped_constraint = Constraint::Atom(AtomId(1), AtomConstraint::degree(3));
        triatomic.push_constraint(dropped_constraint.clone());
        let expected = triatomic.clone().build();
        let removed_atoms = vec![RemovedAtom {
            id: AtomId(1),
            ast: triatomic.atom(AtomId(1)).ast.clone(),
        }];
        let removed_bonds = vec![
            RemovedBond {
                id: BondId(0),
                endpoints: triatomic.bond(BondId(0)).atoms,
                ast: triatomic.bond(BondId(0)).ast.clone(),
            },
            RemovedBond {
                id: BondId(1),
                endpoints: triatomic.bond(BondId(1)).atoms,
                ast: triatomic.bond(BondId(1)).ast.clone(),
            },
        ];

        let compaction = triatomic.remove(&[AtomId(1)], &[]);
        triatomic.restore_topology(
            removed_atoms,
            removed_bonds,
            RemovedOverlays::default(),
            &compaction.undo_compaction(),
            CascadedConstraints {
                removed: vec![RemovedConstraint {
                    position: 0,
                    constraint: dropped_constraint,
                }],
                modified: Vec::new(),
            },
        );

        assert_eq!(triatomic.build(), expected);
    }

    #[rstest]
    fn test_molecule_builder_remove_added_topology(mut triatomic: MoleculeBuilder) {
        let expected = triatomic.clone().build();
        let added_atom = AddedAtom {
            id: triatomic.add_atom(AtomAst::from_element(Element::F)),
            ast: AtomAst::from_element(Element::F),
        };
        let added_bond = AddedBond {
            id: triatomic.add_bond(AtomId(2), added_atom.id, BondAst::from_order(1)),
            endpoints: [AtomId(2), added_atom.id],
            ast: BondAst::from_order(1),
        };

        triatomic.remove_added_topology(&[added_atom], &[added_bond]);

        assert_eq!(triatomic.build(), expected);
    }

    #[rstest]
    fn test_molecule_builder_restore_dative_bond() {
        let mut b = MoleculeAst::default().edit();
        b.add_atom(AtomAst::from_element(Element::C));
        b.add_atom(AtomAst::from_element(Element::N));
        b.add_dative_bond(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1));
        b.add_dative_bond(vec![AtomId(1)], AtomId(0), DativeBondAst::from_order(2));
        let expected = b.clone().build();
        let view = b.dative_bond(DativeBondId(0));
        let removed = RemovedDativeBond {
            id: DativeBondId(0),
            atoms: view.atom_ids().collect(),
            ast: view.ast.clone(),
        };

        b.remove_dative_bonds(&[DativeBondId(0)]);
        let undo = IdCompaction::relations(
            vec![removed.id.into()],
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .undo_compaction();
        b.restore_dative_bonds(vec![removed], &undo);

        assert_eq!(b.build(), expected);
    }

    // `edit()` → `build()` reproduces the AST including both stereo overlays.
    #[rstest]
    fn test_molecule_builder_build() {
        let ast = mol!(
            r#"{:atoms ["C" "C" "C" "F" "Cl"]
                :bonds [[0 1 "1"] [1 2 "2"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 3 4] :type "Th1"}]
                :stereo-bonds [{:site 1 :ligands [0 2] :type "Ct1"}]}"#
        );
        assert_eq!(ast.edit().build(), ast);
    }

    // `remove` forward-compacts stereo-atom node refs: removing a non-participant
    // shifts the surviving site/ligand ids; removing the site drops the element.
    #[rstest]
    #[case::remaps_surviving(vec![AtomId(0)], vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]])]
    #[case::drops_on_site_removal(vec![AtomId(1)], vec![])]
    fn test_molecule_builder_remove_stereo_atom(
        #[case] remove_atoms: Vec<AtomId>,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        let ast = mol!(
            r#"{:atoms ["C" "C" "F" "Cl" "Br"]
                :bonds [[1 2 "1"] [1 3 "1"] [1 4 "1"]]
                :stereo-atoms [{:site 1 :ligands [2 3 4] :type "Th1"}]}"#
        );
        let mut builder = ast.edit();
        builder.remove(&remove_atoms, &[]);
        let surviving: Vec<Vec<AtomId>> = builder
            .build()
            .stereo_atoms()
            .iter()
            .map(|view| view.atom_ids().collect())
            .collect();
        assert_eq!(surviving, expected);
    }

    // `remove` forward-compacts the stereo-bond edge site: removing a non-site bond
    // shifts the surviving site; removing the site bond drops the element.
    #[rstest]
    #[case::remaps_surviving(vec![BondId(0)], vec![BondId(0)])]
    #[case::drops_on_site_removal(vec![BondId(1)], vec![])]
    fn test_molecule_builder_remove_stereo_bond(
        #[case] remove_bonds: Vec<BondId>,
        #[case] expected: Vec<BondId>,
    ) {
        let ast = mol!(
            r#"{:atoms ["C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
                :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"#
        );
        let mut builder = ast.edit();
        builder.remove(&[], &remove_bonds);
        let surviving: Vec<BondId> = builder
            .build()
            .stereo_bonds()
            .iter()
            .map(|view| view.site_id())
            .collect();
        assert_eq!(surviving, expected);
    }
}
