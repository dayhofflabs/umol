//! Structural editing for `Molecule`. The molecule itself only allows attribute
//! mutation; structural change (add atoms/bonds/relations, remove anything)
//! goes through `MoleculeEditor`.
//!
//! Storage is lazy: each Arc-shared field stays shared until first write,
//! at which point only that field decomposes to a mutable form. `build`
//! re-wraps everything in `Arc`, reusing untouched shared data.

use std::collections::HashSet;
use std::mem;
use std::sync::Arc;

use umol_graph_core::{
    Compaction, EdgeId, FixedRelationSet, FixedVarBirelationSet, Graph, GraphCompaction, NodeId,
    RelationId, RelationParticipant, VarRelationSet,
};
use umol_perm::{DynPermutation, Permutation};

use super::super::aromatic::{AromaticSystemForm, AromaticSystems};
use super::super::atom::AtomForm;
use super::super::bond::BondForm;
use super::super::compact::{MoleculeCompaction, UndoCompaction};
use super::super::constraint::{Constraint, Constraints};
use super::super::dative::{DativeBondForm, DativeBonds};
use super::super::edit::{
    AddedAromaticSystem, AddedAtom, AddedBond, AddedDativeBond, AddedMulticenterBond,
    AddedNoncovalentBond, AddedStereoAtom, AddedStereoBond, RemovedAromaticSystem, RemovedAtom,
    RemovedBond, RemovedDativeBond, RemovedMulticenterBond, RemovedNoncovalentBond,
    RemovedOverlays, RemovedStereoAtom, RemovedStereoBond,
};
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::ligand::StereoLigand;
use super::super::multicenter::{MulticenterBondForm, MulticenterBonds};
use super::super::noncovalent::{NoncovalentBondForm, NoncovalentBonds};
use super::super::stereo::{StereoAtomForm, StereoAtoms, StereoBondForm, StereoBonds};
use super::super::traits::{FrameTransport, Normalize};
use super::super::view::{
    AromaticSystemEditorView, AromaticSystemEditorViewMut, AtomEditorView, AtomEditorViewMut,
    BondEditorView, BondEditorViewMut, DativeBondEditorView, DativeBondEditorViewMut,
    MulticenterBondEditorView, MulticenterBondEditorViewMut, NoncovalentBondEditorView,
    NoncovalentBondEditorViewMut, StereoAtomEditorView, StereoAtomEditorViewMut,
    StereoBondEditorView, StereoBondEditorViewMut,
};
use super::{Molecule, MoleculeIntegrityError};

#[derive(Clone)]
enum FixedSetStorage<P, D, const N: usize> {
    Shared(Arc<FixedRelationSet<P, D, N>>),
    Mutable(Vec<([P; N], D)>),
}

impl<P, D, const N: usize> FixedSetStorage<P, D, N>
where
    P: RelationParticipant,
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
        if matches!(self, FixedSetStorage::Shared(_)) {
            let FixedSetStorage::Shared(arc) =
                mem::replace(self, FixedSetStorage::Mutable(Vec::new()))
            else {
                unreachable!()
            };
            let entries = match Arc::try_unwrap(arc) {
                Ok(relation_set) => relation_set.into_entries(),
                Err(arc) => (0..arc.count())
                    .map(|i| {
                        let id = RelationId(i as u32);
                        (*arc.participants(id), arc.data(id).clone())
                    })
                    .collect(),
            };
            *self = FixedSetStorage::Mutable(entries);
        }
    }

    fn into_arc(self) -> Arc<FixedRelationSet<P, D, N>> {
        match self {
            FixedSetStorage::Shared(arc) => arc,
            FixedSetStorage::Mutable(vec) => Arc::new(FixedRelationSet::new(vec)),
        }
    }

    fn count(&self) -> usize {
        match self {
            FixedSetStorage::Shared(arc) => arc.count(),
            FixedSetStorage::Mutable(vec) => vec.len(),
        }
    }

    fn participants(&self, i: usize) -> [P; N] {
        match self {
            FixedSetStorage::Shared(arc) => *arc.participants(RelationId(i as u32)),
            FixedSetStorage::Mutable(vec) => vec[i].0,
        }
    }

    fn data(&self, i: usize) -> D {
        match self {
            FixedSetStorage::Shared(arc) => arc.data(RelationId(i as u32)).clone(),
            FixedSetStorage::Mutable(vec) => vec[i].1.clone(),
        }
    }

    /// Whether relation `i` coincides with `query` — multiset equality of the participants.
    ///
    /// The known-id sibling of a coincidence search, and the same operation graph-core exposes as
    /// `is_coincident`. Mutable storage is not kept canonical, so both sides are sorted.
    fn is_coincident(&self, i: usize, query: &[P]) -> bool {
        let stored = self.participants(i);
        stored.len() == query.len() && {
            let mut stored_sorted = stored.to_vec();
            stored_sorted.sort_unstable();
            let mut query_sorted = query.to_vec();
            query_sorted.sort_unstable();
            stored_sorted == query_sorted
        }
    }

    fn compact(self, compaction: &GraphCompaction) -> (Self, Compaction<RelationId>) {
        match self {
            FixedSetStorage::Shared(arc) => {
                let (compacted, removed) = arc.tracked_compact(compaction);
                (FixedSetStorage::Shared(Arc::new(compacted)), removed)
            }
            FixedSetStorage::Mutable(vec) => {
                let source_count = vec.len();
                let mut removed = Vec::new();
                let mut compacted: Vec<([P; N], D)> = Vec::with_capacity(vec.len());
                for (index, (mut participants, d)) in vec.into_iter().enumerate() {
                    let survives = participants
                        .iter_mut()
                        .all(|participant| match (*participant).compact(compaction) {
                            Some(mapped) => {
                                *participant = mapped;
                                true
                            }
                            None => false,
                        });
                    if survives {
                        compacted.push((participants, d));
                    } else {
                        removed.push(RelationId(index as u32));
                    }
                }
                (
                    FixedSetStorage::Mutable(compacted),
                    Compaction::new(source_count, removed)
                        .expect("removed relations belong to the source set"),
                )
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
            FixedSetStorage::Shared(arc) => (0..arc.count())
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
enum VarSetStorage<P, D> {
    Shared(Arc<VarRelationSet<P, D>>),
    Mutable(Vec<(Vec<P>, D)>),
}

impl<P, D> VarSetStorage<P, D>
where
    P: RelationParticipant,
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
        if matches!(self, VarSetStorage::Shared(_)) {
            let VarSetStorage::Shared(arc) = mem::replace(self, VarSetStorage::Mutable(Vec::new()))
            else {
                unreachable!()
            };
            let entries = match Arc::try_unwrap(arc) {
                Ok(relation_set) => relation_set.into_entries(),
                Err(arc) => (0..arc.count())
                    .map(|i| {
                        let id = RelationId(i as u32);
                        (arc.participants(id).to_vec(), arc.data(id).clone())
                    })
                    .collect(),
            };
            *self = VarSetStorage::Mutable(entries);
        }
    }

    fn into_arc(self) -> Arc<VarRelationSet<P, D>> {
        match self {
            VarSetStorage::Shared(arc) => arc,
            VarSetStorage::Mutable(vec) => Arc::new(VarRelationSet::new(vec)),
        }
    }

    fn count(&self) -> usize {
        match self {
            VarSetStorage::Shared(arc) => arc.count(),
            VarSetStorage::Mutable(vec) => vec.len(),
        }
    }

    fn participants(&self, i: usize) -> Vec<P> {
        match self {
            VarSetStorage::Shared(arc) => arc.participants(RelationId(i as u32)).to_vec(),
            VarSetStorage::Mutable(vec) => vec[i].0.clone(),
        }
    }

    /// Whether relation `i` coincides with `query` — multiset equality of the participants.
    ///
    /// The known-id sibling of a coincidence search, and the same operation graph-core exposes as
    /// `is_coincident`. Mutable storage is not kept canonical, so both sides are sorted.
    fn is_coincident(&self, i: usize, query: &[P]) -> bool {
        let stored = self.participants(i);
        stored.len() == query.len() && {
            let mut stored_sorted = stored.to_vec();
            stored_sorted.sort_unstable();
            let mut query_sorted = query.to_vec();
            query_sorted.sort_unstable();
            stored_sorted == query_sorted
        }
    }

    fn data(&self, i: usize) -> D {
        match self {
            VarSetStorage::Shared(arc) => arc.data(RelationId(i as u32)).clone(),
            VarSetStorage::Mutable(vec) => vec[i].1.clone(),
        }
    }

    fn compact(self, compaction: &GraphCompaction) -> (Self, Compaction<RelationId>) {
        match self {
            VarSetStorage::Shared(arc) => {
                let (compacted, removed) = arc.tracked_compact(compaction);
                (VarSetStorage::Shared(Arc::new(compacted)), removed)
            }
            VarSetStorage::Mutable(vec) => {
                let source_count = vec.len();
                let mut removed = Vec::new();
                let mut compacted: Vec<(Vec<P>, D)> = Vec::with_capacity(vec.len());
                for (index, (participants, d)) in vec.into_iter().enumerate() {
                    let mapped: Option<Vec<P>> = participants
                        .into_iter()
                        .map(|p| p.compact(compaction))
                        .collect();
                    match mapped {
                        Some(participants) => compacted.push((participants, d)),
                        None => removed.push(RelationId(index as u32)),
                    }
                }
                (
                    VarSetStorage::Mutable(compacted),
                    Compaction::new(source_count, removed)
                        .expect("removed relations belong to the source set"),
                )
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
            VarSetStorage::Shared(arc) => (0..arc.count())
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
enum FixedVarSetStorage<L1, const N1: usize, L2, D> {
    Shared(Arc<FixedVarBirelationSet<L1, N1, L2, D>>),
    Mutable(Vec<([L1; N1], Vec<L2>, D)>),
}

impl<L1, const N1: usize, L2, D> FixedVarSetStorage<L1, N1, L2, D>
where
    L1: RelationParticipant,
    L2: RelationParticipant,
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
        if matches!(self, FixedVarSetStorage::Shared(_)) {
            let FixedVarSetStorage::Shared(arc) =
                mem::replace(self, FixedVarSetStorage::Mutable(Vec::new()))
            else {
                unreachable!()
            };
            let entries = match Arc::try_unwrap(arc) {
                Ok(relation_set) => relation_set.into_entries(),
                Err(arc) => arc_entries(&arc),
            };
            *self = FixedVarSetStorage::Mutable(entries);
        }
    }

    fn into_arc(self) -> Arc<FixedVarBirelationSet<L1, N1, L2, D>> {
        match self {
            FixedVarSetStorage::Shared(arc) => arc,
            FixedVarSetStorage::Mutable(vec) => Arc::new(FixedVarBirelationSet::new(vec)),
        }
    }

    fn count(&self) -> usize {
        match self {
            FixedVarSetStorage::Shared(arc) => arc.count(),
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

    fn data(&self, i: usize) -> D {
        match self {
            FixedVarSetStorage::Shared(arc) => arc.data(RelationId(i as u32)).clone(),
            FixedVarSetStorage::Mutable(vec) => vec[i].2.clone(),
        }
    }

    /// Per-factor permutations reindexing `(query_1, query_2)` into relation `i`'s stored participant
    /// order (`σ[k]` = the position in the query of the participant equal to `stored[k]`), or `None`
    /// when either factor's sets differ. Direct alignment — mutable storage is not kept canonical.
    #[allow(clippy::type_complexity)]
    /// Whether relation `i` coincides with `query_1` / `query_2` — multiset equality of each factor.
    ///
    /// The known-id sibling of a coincidence search, and the same operation graph-core exposes as
    /// `is_coincident`. Mutable storage is not kept canonical, so both sides are sorted.
    fn is_coincident(&self, i: usize, query_1: &[L1], query_2: &[L2]) -> bool {
        let (stored_1, stored_2) = (self.participants_1(i), self.participants_2(i));
        stored_1.len() == query_1.len() && stored_2.len() == query_2.len() && {
            let mut s1 = stored_1.to_vec();
            s1.sort_unstable();
            let mut q1 = query_1.to_vec();
            q1.sort_unstable();
            let mut s2 = stored_2.to_vec();
            s2.sort_unstable();
            let mut q2 = query_2.to_vec();
            q2.sort_unstable();
            s1 == q1 && s2 == q2
        }
    }

    fn compact(self, compaction: &GraphCompaction) -> (Self, Compaction<RelationId>) {
        match self {
            FixedVarSetStorage::Shared(arc) => {
                let (compacted, removed) = arc.tracked_compact(compaction);
                (FixedVarSetStorage::Shared(Arc::new(compacted)), removed)
            }
            FixedVarSetStorage::Mutable(vec) => {
                let source_count = vec.len();
                let mut removed = Vec::new();
                let mut compacted: Vec<([L1; N1], Vec<L2>, D)> = Vec::with_capacity(vec.len());
                for (index, (mut participants_1, participants_2, d)) in vec.into_iter().enumerate()
                {
                    let f1 = participants_1.iter_mut().all(|participant| {
                        match (*participant).compact(compaction) {
                            Some(mapped) => {
                                *participant = mapped;
                                true
                            }
                            None => false,
                        }
                    });
                    let f2: Option<Vec<L2>> = participants_2
                        .into_iter()
                        .map(|p| p.compact(compaction))
                        .collect();
                    match (f1, f2) {
                        (true, Some(participants_2)) => {
                            compacted.push((participants_1, participants_2, d))
                        }
                        _ => removed.push(RelationId(index as u32)),
                    }
                }
                (
                    FixedVarSetStorage::Mutable(compacted),
                    Compaction::new(source_count, removed)
                        .expect("removed relations belong to the source set"),
                )
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

fn arc_entries<L1, const N1: usize, L2, D>(
    arc: &FixedVarBirelationSet<L1, N1, L2, D>,
) -> Vec<([L1; N1], Vec<L2>, D)>
where
    L1: RelationParticipant,
    L2: RelationParticipant,
    D: Clone,
{
    (0..arc.count())
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

/// Mutable editor for a `Molecule`. Accumulates atoms, bonds, and
/// relations (dative, aromatic, multicenter, noncovalent), then finalizes
/// into an immutable `Molecule`. Supports incremental removal with
/// index remapping via `remove`.
#[derive(Clone)]
pub struct MoleculeEditor {
    graph: Graph,
    atoms: Arc<Vec<AtomForm>>,
    bonds: Arc<Vec<BondForm>>,
    dative_bonds: FixedVarSetStorage<NodeId, 1, NodeId, DativeBondForm>,
    aromatic_systems: VarSetStorage<NodeId, AromaticSystemForm>,
    multicenter_bonds: VarSetStorage<NodeId, MulticenterBondForm>,
    noncovalent_bonds: FixedSetStorage<NodeId, NoncovalentBondForm, 2>,
    stereo_atoms: FixedVarSetStorage<NodeId, 1, StereoLigand, StereoAtomForm>,
    stereo_bonds: FixedVarSetStorage<EdgeId, 1, StereoLigand, StereoBondForm>,
    constraints: Constraints,
}

impl MoleculeEditor {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn from_parts(
        graph: Graph,
        atoms: Arc<Vec<AtomForm>>,
        bonds: Arc<Vec<BondForm>>,
        dative_bonds: DativeBonds,
        aromatic_systems: AromaticSystems,
        multicenter_bonds: MulticenterBonds,
        noncovalent_bonds: NoncovalentBonds,
        stereo_atoms: StereoAtoms,
        stereo_bonds: StereoBonds,
        constraints: Constraints,
    ) -> Self {
        Self {
            graph,
            atoms,
            bonds,
            dative_bonds: FixedVarSetStorage::Shared(dative_bonds.into_arc()),
            aromatic_systems: VarSetStorage::Shared(aromatic_systems.into_arc()),
            multicenter_bonds: VarSetStorage::Shared(multicenter_bonds.into_arc()),
            noncovalent_bonds: FixedSetStorage::Shared(noncovalent_bonds.into_arc()),
            stereo_atoms: FixedVarSetStorage::Shared(stereo_atoms.into_arc()),
            stereo_bonds: FixedVarSetStorage::Shared(stereo_bonds.into_arc()),
            constraints,
        }
    }

    /// Append an atom directly to the editor.
    ///
    /// This is a low-level, non-transactional construction primitive. Use `transact` for checked
    /// atomic edits with rollback or `apply` for consuming application without an undo journal.
    pub fn add_atom(&mut self, atom: AtomForm) -> AtomId {
        let id = self.graph.add_node();
        Arc::make_mut(&mut self.atoms).push(atom);
        AtomId::from(id)
    }

    /// Append a localized bond directly to the editor.
    ///
    /// This is a low-level, non-transactional construction primitive. It
    /// assumes `first` and `second` are valid atom ids in the current dense layout.
    pub fn add_bond(&mut self, first: AtomId, second: AtomId, bond: BondForm) -> BondId {
        let id = self
            .graph
            .add_edge(NodeId::from(first), NodeId::from(second));
        Arc::make_mut(&mut self.bonds).push(bond);
        BondId::from(id)
    }

    /// Append a dative-bond overlay directly to the editor. The acceptor is
    /// factor 1; the donors are factor 2 (sorted by the `Unordered`
    /// canonicalization).
    pub fn add_dative_bond(
        &mut self,
        donors: Vec<AtomId>,
        acceptor: AtomId,
        bond: DativeBondForm,
    ) -> DativeBondId {
        let donors: Vec<NodeId> = donors.into_iter().map(NodeId::from).collect();
        DativeBondId(
            self.dative_bonds
                .push([NodeId::from(acceptor)], donors, bond),
        )
    }

    /// Append an aromatic-system overlay directly to the editor.
    pub fn add_aromatic_system(
        &mut self,
        atoms: Vec<AtomId>,
        data: AromaticSystemForm,
    ) -> AromaticSystemId {
        let nodes: Vec<NodeId> = atoms.into_iter().map(NodeId::from).collect();
        let i = self.aromatic_systems.push(nodes, data);
        AromaticSystemId(i)
    }

    /// Append a multicenter-bond overlay directly to the editor.
    pub fn add_multicenter_bond(
        &mut self,
        atoms: Vec<AtomId>,
        data: MulticenterBondForm,
    ) -> MulticenterBondId {
        let nodes: Vec<NodeId> = atoms.into_iter().map(NodeId::from).collect();
        let i = self.multicenter_bonds.push(nodes, data);
        MulticenterBondId(i)
    }

    /// Append a noncovalent-bond overlay directly to the editor.
    pub fn add_noncovalent_bond(
        &mut self,
        ends: [AtomId; 2],
        bond: NoncovalentBondForm,
    ) -> NoncovalentBondId {
        let i = self
            .noncovalent_bonds
            .push([NodeId::from(ends[0]), NodeId::from(ends[1])], bond);
        NoncovalentBondId(i)
    }

    /// Append a stereo-atom overlay directly to the editor.
    pub fn add_stereo_atom(
        &mut self,
        site: AtomId,
        ligands: Vec<StereoLigand>,
        attributes: StereoAtomForm,
    ) -> StereoAtomId {
        StereoAtomId(
            self.stereo_atoms
                .push([NodeId::from(site)], ligands, attributes),
        )
    }

    /// Append a stereo-bond overlay directly to the editor.
    pub fn add_stereo_bond(
        &mut self,
        site: BondId,
        ligands: Vec<StereoLigand>,
        attributes: StereoBondForm,
    ) -> StereoBondId {
        StereoBondId(
            self.stereo_bonds
                .push([EdgeId::from(site)], ligands, attributes),
        )
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
    // the editor itself because dense removal can compact many unrelated ids.

    pub fn atom(&self, id: AtomId) -> AtomEditorView<'_> {
        AtomEditorView {
            id,
            attributes: &self.atoms[id.index()],
        }
    }

    pub fn atom_mut(&mut self, id: AtomId) -> AtomEditorViewMut<'_> {
        let attributes = &mut Arc::make_mut(&mut self.atoms)[id.index()];
        AtomEditorViewMut { id, attributes }
    }

    pub fn bond(&self, id: BondId) -> BondEditorView<'_> {
        let endpoints = self.graph.edge_endpoints(EdgeId::from(id));
        let atoms = [AtomId::from(endpoints[0]), AtomId::from(endpoints[1])];
        BondEditorView {
            id,
            attributes: &self.bonds[id.index()],
            atoms,
        }
    }

    pub fn bond_mut(&mut self, id: BondId) -> BondEditorViewMut<'_> {
        let endpoints = self.graph.edge_endpoints(EdgeId::from(id));
        let atoms = [AtomId::from(endpoints[0]), AtomId::from(endpoints[1])];
        let attributes = &mut Arc::make_mut(&mut self.bonds)[id.index()];
        BondEditorViewMut {
            id,
            attributes,
            atoms,
        }
    }

    pub fn dative_bond(&self, id: DativeBondId) -> DativeBondEditorView<'_> {
        match &self.dative_bonds {
            FixedVarSetStorage::Shared(arc) => {
                let rid = RelationId(id.0);
                DativeBondEditorView::new(
                    id,
                    arc.participants_2(rid),
                    AtomId::from(arc.participants_1(rid)[0]),
                    arc.data(rid),
                )
            }
            FixedVarSetStorage::Mutable(vec) => {
                let entry = &vec[id.index()];
                DativeBondEditorView::new(id, &entry.1, AtomId::from(entry.0[0]), &entry.2)
            }
        }
    }

    pub fn dative_bond_mut(&mut self, id: DativeBondId) -> DativeBondEditorViewMut<'_> {
        self.dative_bonds.materialize();
        let FixedVarSetStorage::Mutable(vec) = &mut self.dative_bonds else {
            unreachable!()
        };
        let entry = &mut vec[id.index()];
        let acceptor = AtomId::from(entry.0[0]);
        DativeBondEditorViewMut::new(id, &entry.1, acceptor, &mut entry.2)
    }

    pub fn aromatic_system(&self, id: AromaticSystemId) -> AromaticSystemEditorView<'_> {
        match &self.aromatic_systems {
            VarSetStorage::Shared(arc) => {
                let rid = RelationId(id.0);
                AromaticSystemEditorView::new(id, arc.participants(rid), arc.data(rid))
            }
            VarSetStorage::Mutable(vec) => {
                let entry = &vec[id.index()];
                AromaticSystemEditorView::new(id, &entry.0, &entry.1)
            }
        }
    }

    pub fn aromatic_system_mut(&mut self, id: AromaticSystemId) -> AromaticSystemEditorViewMut<'_> {
        self.aromatic_systems.materialize();
        let VarSetStorage::Mutable(vec) = &mut self.aromatic_systems else {
            unreachable!()
        };
        let entry = &mut vec[id.index()];
        AromaticSystemEditorViewMut::new(id, &entry.0, &mut entry.1)
    }

    pub fn multicenter_bond(&self, id: MulticenterBondId) -> MulticenterBondEditorView<'_> {
        match &self.multicenter_bonds {
            VarSetStorage::Shared(arc) => {
                let rid = RelationId(id.0);
                MulticenterBondEditorView::new(id, arc.participants(rid), arc.data(rid))
            }
            VarSetStorage::Mutable(vec) => {
                let entry = &vec[id.index()];
                MulticenterBondEditorView::new(id, &entry.0, &entry.1)
            }
        }
    }

    pub fn multicenter_bond_mut(
        &mut self,
        id: MulticenterBondId,
    ) -> MulticenterBondEditorViewMut<'_> {
        self.multicenter_bonds.materialize();
        let VarSetStorage::Mutable(vec) = &mut self.multicenter_bonds else {
            unreachable!()
        };
        let entry = &mut vec[id.index()];
        MulticenterBondEditorViewMut::new(id, &entry.0, &mut entry.1)
    }

    pub fn noncovalent_bond(&self, id: NoncovalentBondId) -> NoncovalentBondEditorView<'_> {
        match &self.noncovalent_bonds {
            FixedSetStorage::Shared(arc) => {
                let rid = RelationId(id.0);
                let parts = arc.participants(rid);
                NoncovalentBondEditorView {
                    id,
                    attributes: arc.data(rid),
                    atoms: [AtomId::from(parts[0]), AtomId::from(parts[1])],
                }
            }
            FixedSetStorage::Mutable(vec) => {
                let entry = &vec[id.index()];
                NoncovalentBondEditorView {
                    id,
                    attributes: &entry.1,
                    atoms: [AtomId::from(entry.0[0]), AtomId::from(entry.0[1])],
                }
            }
        }
    }

    pub fn noncovalent_bond_mut(
        &mut self,
        id: NoncovalentBondId,
    ) -> NoncovalentBondEditorViewMut<'_> {
        self.noncovalent_bonds.materialize();
        let FixedSetStorage::Mutable(vec) = &mut self.noncovalent_bonds else {
            unreachable!()
        };
        let entry = &mut vec[id.index()];
        let atoms = [AtomId::from(entry.0[0]), AtomId::from(entry.0[1])];
        NoncovalentBondEditorViewMut {
            id,
            attributes: &mut entry.1,
            atoms,
        }
    }

    pub fn stereo_atom(&self, id: StereoAtomId) -> StereoAtomEditorView<'_> {
        match &self.stereo_atoms {
            FixedVarSetStorage::Shared(arc) => {
                let rid = RelationId(id.0);
                StereoAtomEditorView {
                    id,
                    attributes: arc.data(rid),
                    site: AtomId::from(arc.participants_1(rid)[0]),
                    ligands: arc.participants_2(rid),
                }
            }
            FixedVarSetStorage::Mutable(vec) => {
                let entry = &vec[id.index()];
                StereoAtomEditorView {
                    id,
                    attributes: &entry.2,
                    site: AtomId::from(entry.0[0]),
                    ligands: &entry.1,
                }
            }
        }
    }

    pub fn stereo_bond(&self, id: StereoBondId) -> StereoBondEditorView<'_> {
        match &self.stereo_bonds {
            FixedVarSetStorage::Shared(arc) => {
                let rid = RelationId(id.0);
                StereoBondEditorView {
                    id,
                    attributes: arc.data(rid),
                    site: BondId::from(arc.participants_1(rid)[0]),
                    ligands: arc.participants_2(rid),
                }
            }
            FixedVarSetStorage::Mutable(vec) => {
                let entry = &vec[id.index()];
                StereoBondEditorView {
                    id,
                    attributes: &entry.2,
                    site: BondId::from(entry.0[0]),
                    ligands: &entry.1,
                }
            }
        }
    }

    /// `true` iff noncovalent bond `id` structurally equals `(atoms, attributes)` — participants (unordered)
    /// and `attributes` up to normal form, `attributes` reindexed into the stored participant frame.
    pub(crate) fn noncovalent_bond_equiv(
        &self,
        id: NoncovalentBondId,
        atoms: [AtomId; 2],
        attributes: &NoncovalentBondForm,
    ) -> bool {
        let stored: Vec<AtomId> = self
            .noncovalent_bonds
            .participants(id.index())
            .iter()
            .map(|&node| AtomId::from(node))
            .collect();
        self.noncovalent_bonds
            .is_coincident(id.index(), &atoms.map(NodeId::from))
            && DynPermutation::between(&atoms, &stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| {
                    restated.normalized_eq(&self.noncovalent_bonds.data(id.index()))
                })
    }

    /// `true` iff aromatic system `id` structurally equals `(atoms, attributes)`.
    pub(crate) fn aromatic_system_equiv(
        &self,
        id: AromaticSystemId,
        atoms: &[AtomId],
        attributes: &AromaticSystemForm,
    ) -> bool {
        let stored: Vec<AtomId> = self
            .aromatic_systems
            .participants(id.index())
            .iter()
            .map(|&node| AtomId::from(node))
            .collect();
        self.aromatic_systems.is_coincident(
            id.index(),
            &atoms.iter().map(|&a| NodeId::from(a)).collect::<Vec<_>>(),
        ) && DynPermutation::between(atoms, &stored)
            .and_then(|action| attributes.clone().reframe_by(&action))
            .is_some_and(|restated| restated.normalized_eq(&self.aromatic_systems.data(id.index())))
    }

    /// `true` iff multicenter bond `id` structurally equals `(atoms, attributes)`.
    pub(crate) fn multicenter_bond_equiv(
        &self,
        id: MulticenterBondId,
        atoms: &[AtomId],
        attributes: &MulticenterBondForm,
    ) -> bool {
        let stored: Vec<AtomId> = self
            .multicenter_bonds
            .participants(id.index())
            .iter()
            .map(|&node| AtomId::from(node))
            .collect();
        self.multicenter_bonds.is_coincident(
            id.index(),
            &atoms.iter().map(|&a| NodeId::from(a)).collect::<Vec<_>>(),
        ) && DynPermutation::between(atoms, &stored)
            .and_then(|action| attributes.clone().reframe_by(&action))
            .is_some_and(|restated| {
                restated.normalized_eq(&self.multicenter_bonds.data(id.index()))
            })
    }

    /// `true` iff dative bond `id` structurally equals `(acceptor, donors, attributes)` — the acceptor
    /// (ordered, single) and donors (unordered) factors and `attributes` up to normal form.
    pub(crate) fn dative_bond_equiv(
        &self,
        id: DativeBondId,
        acceptor: AtomId,
        donors: &[AtomId],
        attributes: &DativeBondForm,
    ) -> bool {
        let stored: Vec<AtomId> = self
            .dative_bonds
            .participants_2(id.index())
            .iter()
            .map(|&node| AtomId::from(node))
            .collect();
        let donor_nodes: Vec<NodeId> = donors.iter().map(|&a| NodeId::from(a)).collect();
        self.dative_bonds
            .is_coincident(id.index(), &[NodeId::from(acceptor)], &donor_nodes)
            && DynPermutation::between(donors, &stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(&self.dative_bonds.data(id.index())))
    }

    /// `true` iff stereo atom `id` structurally equals `(site, ligands, attributes)`.
    pub(crate) fn stereo_atom_equiv(
        &self,
        id: StereoAtomId,
        site: AtomId,
        ligands: &[StereoLigand],
        attributes: &StereoAtomForm,
    ) -> bool {
        let stored = self.stereo_atoms.participants_2(id.index());
        AtomId::from(self.stereo_atoms.participants_1(id.index())[0]) == site
            && Permutation::between(ligands, &stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(&self.stereo_atoms.data(id.index())))
    }

    /// `true` iff stereo bond `id` structurally equals `(site, ligands, attributes)`.
    pub(crate) fn stereo_bond_equiv(
        &self,
        id: StereoBondId,
        site: BondId,
        ligands: &[StereoLigand],
        attributes: &StereoBondForm,
    ) -> bool {
        let stored = self.stereo_bonds.participants_2(id.index());
        BondId::from(self.stereo_bonds.participants_1(id.index())[0]) == site
            && Permutation::between(ligands, &stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(&self.stereo_bonds.data(id.index())))
    }

    pub fn stereo_atom_mut(&mut self, id: StereoAtomId) -> StereoAtomEditorViewMut<'_> {
        self.stereo_atoms.materialize();
        let FixedVarSetStorage::Mutable(vec) = &mut self.stereo_atoms else {
            unreachable!()
        };
        let entry = &mut vec[id.index()];
        let site = AtomId::from(entry.0[0]);
        StereoAtomEditorViewMut {
            id,
            attributes: &mut entry.2,
            site,
            ligands: &entry.1,
        }
    }

    pub fn stereo_bond_mut(&mut self, id: StereoBondId) -> StereoBondEditorViewMut<'_> {
        self.stereo_bonds.materialize();
        let FixedVarSetStorage::Mutable(vec) = &mut self.stereo_bonds else {
            unreachable!()
        };
        let entry = &mut vec[id.index()];
        let site = BondId::from(entry.0[0]);
        StereoBondEditorViewMut {
            id,
            attributes: &mut entry.2,
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
        self.dative_bonds.count()
    }

    pub fn aromatic_system_count(&self) -> usize {
        self.aromatic_systems.count()
    }

    pub fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.count()
    }

    pub fn noncovalent_bond_count(&self) -> usize {
        self.noncovalent_bonds.count()
    }

    pub fn stereo_atom_count(&self) -> usize {
        self.stereo_atoms.count()
    }

    pub fn stereo_bond_count(&self) -> usize {
        self.stereo_bonds.count()
    }

    // -- Relation removal -----------------------------------------------------

    /// Remove dative-bond overlays directly from the editor.
    ///
    /// This is a low-level dense removal primitive. It compacts molecule-level
    /// constraints but does not build rollback data.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn remove_dative_bonds(&mut self, ids: &[DativeBondId]) {
        self.tracked_remove_dative_bonds(ids);
    }

    /// Remove overlays and return the source-to-result compaction for all entity kinds.
    ///
    /// Leaves the same state as [`Self::remove_dative_bonds`]; unchanged families retain their counts.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn tracked_remove_dative_bonds(&mut self, ids: &[DativeBondId]) -> MoleculeCompaction {
        let raw: Vec<RelationId> = ids.iter().map(|&i| i.into()).collect();
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(self.atom_count()),
                Compaction::identity(self.bond_count()),
            ),
            Compaction::new(self.dative_bond_count(), ids.to_vec())
                .expect("removed entities belong to the source table"),
            Compaction::identity(self.aromatic_system_count()),
            Compaction::identity(self.multicenter_bond_count()),
            Compaction::identity(self.noncovalent_bond_count()),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.dative_bonds.remove_relations(&raw);
        self.constraints.compact(&compaction);
        compaction
    }

    /// Remove aromatic-system overlays directly from the editor.
    ///
    /// This is a low-level dense removal primitive. It compacts molecule-level
    /// constraints but does not build rollback data.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn remove_aromatic_systems(&mut self, ids: &[AromaticSystemId]) {
        self.tracked_remove_aromatic_systems(ids);
    }

    /// Remove overlays and return the source-to-result compaction for all entity kinds.
    ///
    /// Leaves the same state as [`Self::remove_aromatic_systems`]; unchanged families retain their counts.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn tracked_remove_aromatic_systems(
        &mut self,
        ids: &[AromaticSystemId],
    ) -> MoleculeCompaction {
        let raw: Vec<RelationId> = ids.iter().map(|&i| i.into()).collect();
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(self.atom_count()),
                Compaction::identity(self.bond_count()),
            ),
            Compaction::identity(self.dative_bond_count()),
            Compaction::new(self.aromatic_system_count(), ids.to_vec())
                .expect("removed entities belong to the source table"),
            Compaction::identity(self.multicenter_bond_count()),
            Compaction::identity(self.noncovalent_bond_count()),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.aromatic_systems.remove_relations(&raw);
        self.constraints.compact(&compaction);
        compaction
    }

    /// Remove multicenter-bond overlays directly from the editor.
    ///
    /// This is a low-level dense removal primitive. It compacts molecule-level
    /// constraints but does not build rollback data.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn remove_multicenter_bonds(&mut self, ids: &[MulticenterBondId]) {
        self.tracked_remove_multicenter_bonds(ids);
    }

    /// Remove overlays and return the source-to-result compaction for all entity kinds.
    ///
    /// Leaves the same state as [`Self::remove_multicenter_bonds`]; unchanged families retain their counts.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn tracked_remove_multicenter_bonds(
        &mut self,
        ids: &[MulticenterBondId],
    ) -> MoleculeCompaction {
        let raw: Vec<RelationId> = ids.iter().map(|&i| i.into()).collect();
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(self.atom_count()),
                Compaction::identity(self.bond_count()),
            ),
            Compaction::identity(self.dative_bond_count()),
            Compaction::identity(self.aromatic_system_count()),
            Compaction::new(self.multicenter_bond_count(), ids.to_vec())
                .expect("removed entities belong to the source table"),
            Compaction::identity(self.noncovalent_bond_count()),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.multicenter_bonds.remove_relations(&raw);
        self.constraints.compact(&compaction);
        compaction
    }

    /// Remove noncovalent-bond overlays directly from the editor.
    ///
    /// This is a low-level dense removal primitive. It compacts molecule-level
    /// constraints but does not build rollback data.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn remove_noncovalent_bonds(&mut self, ids: &[NoncovalentBondId]) {
        self.tracked_remove_noncovalent_bonds(ids);
    }

    /// Remove overlays and return the source-to-result compaction for all entity kinds.
    ///
    /// Leaves the same state as [`Self::remove_noncovalent_bonds`]; unchanged families retain their counts.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn tracked_remove_noncovalent_bonds(
        &mut self,
        ids: &[NoncovalentBondId],
    ) -> MoleculeCompaction {
        let raw: Vec<RelationId> = ids.iter().map(|&i| i.into()).collect();
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(self.atom_count()),
                Compaction::identity(self.bond_count()),
            ),
            Compaction::identity(self.dative_bond_count()),
            Compaction::identity(self.aromatic_system_count()),
            Compaction::identity(self.multicenter_bond_count()),
            Compaction::new(self.noncovalent_bond_count(), ids.to_vec())
                .expect("removed entities belong to the source table"),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.noncovalent_bonds.remove_relations(&raw);
        self.constraints.compact(&compaction);
        compaction
    }

    /// Remove stereo-atom overlays directly from the editor.
    ///
    /// Low-level dense removal primitive; compacts molecule-level constraints but
    /// does not build rollback data.
    pub fn remove_stereo_atoms(&mut self, ids: &[StereoAtomId]) {
        self.tracked_remove_stereo_atoms(ids);
    }

    /// Remove overlays and return the source-to-result compaction for all entity kinds.
    ///
    /// Leaves the same state as [`Self::remove_stereo_atoms`]; unchanged families retain their counts.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn tracked_remove_stereo_atoms(&mut self, ids: &[StereoAtomId]) -> MoleculeCompaction {
        let raw: Vec<RelationId> = ids.iter().map(|&i| i.into()).collect();
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(self.atom_count()),
                Compaction::identity(self.bond_count()),
            ),
            Compaction::identity(self.dative_bond_count()),
            Compaction::identity(self.aromatic_system_count()),
            Compaction::identity(self.multicenter_bond_count()),
            Compaction::identity(self.noncovalent_bond_count()),
            Compaction::new(self.stereo_atom_count(), ids.to_vec())
                .expect("removed entities belong to the source table"),
            Compaction::identity(self.stereo_bond_count()),
        );
        self.stereo_atoms.remove_relations(&raw);
        self.constraints.compact(&compaction);
        compaction
    }

    /// Remove stereo-bond overlays directly from the editor.
    ///
    /// Low-level dense removal primitive; compacts molecule-level constraints but
    /// does not build rollback data.
    pub fn remove_stereo_bonds(&mut self, ids: &[StereoBondId]) {
        self.tracked_remove_stereo_bonds(ids);
    }

    /// Remove overlays and return the source-to-result compaction for all entity kinds.
    ///
    /// Leaves the same state as [`Self::remove_stereo_bonds`]; unchanged families retain their counts.
    ///
    /// # Panics
    ///
    /// Panics when a supplied id is outside the current relation table.
    pub fn tracked_remove_stereo_bonds(&mut self, ids: &[StereoBondId]) -> MoleculeCompaction {
        let raw: Vec<RelationId> = ids.iter().map(|&i| i.into()).collect();
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::identity(self.atom_count()),
                Compaction::identity(self.bond_count()),
            ),
            Compaction::identity(self.dative_bond_count()),
            Compaction::identity(self.aromatic_system_count()),
            Compaction::identity(self.multicenter_bond_count()),
            Compaction::identity(self.noncovalent_bond_count()),
            Compaction::identity(self.stereo_atom_count()),
            Compaction::new(self.stereo_bond_count(), ids.to_vec())
                .expect("removed entities belong to the source table"),
        );
        self.stereo_bonds.remove_relations(&raw);
        self.constraints.compact(&compaction);
        compaction
    }

    // -- Topological removal --------------------------------------------------

    /// Remove topology directly from the editor, cascading dependent relations and constraints.
    ///
    /// This is the low-level dense topology-removal primitive. It removes the
    /// requested atoms and bonds, cascades relations whose participants were
    /// removed, and compacts molecule-level constraints. It does not build rollback
    /// data; checked transactions capture the removed payloads before calling
    /// this method.
    pub fn remove(&mut self, atoms: &[AtomId], bonds: &[BondId]) {
        self.tracked_remove(atoms, bonds);
    }

    /// Remove topology and return the source-to-result compaction for all eight entity kinds.
    ///
    /// Leaves the same state as [`Self::remove`], including cascading relation and constraint
    /// removal. Every component retains the source count from before removal.
    pub fn tracked_remove(&mut self, atoms: &[AtomId], bonds: &[BondId]) -> MoleculeCompaction {
        let nodes: Vec<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        let edges: Vec<EdgeId> = bonds.iter().map(|&b| EdgeId::from(b)).collect();
        let compaction = self.graph.tracked_remove_cascading(&nodes, &edges);

        let new_atoms = compaction.nodes().compact_vec(&self.atoms);
        let new_bonds = compaction.edges().compact_vec(&self.bonds);
        self.atoms = Arc::new(new_atoms);
        self.bonds = Arc::new(new_bonds);

        // Each entity set reports the relation ids its own compaction consumed, so the drop is
        // discovered once rather than traversed separately. A stereo element whose site or any
        // ligand atom or bond was removed drops out the same way (cascade), and the reported ids
        // feed `MoleculeCompaction` so rollback (`restore_topology`) can reinsert them.
        let dative = mem::replace(
            &mut self.dative_bonds,
            FixedVarSetStorage::Shared(Arc::new(FixedVarBirelationSet::default())),
        );
        let (dative, removed_dative) = dative.compact(&compaction);
        self.dative_bonds = dative;

        let aromatic = mem::replace(
            &mut self.aromatic_systems,
            VarSetStorage::Shared(Arc::new(VarRelationSet::default())),
        );
        let (aromatic, removed_aromatic) = aromatic.compact(&compaction);
        self.aromatic_systems = aromatic;

        let multicenter = mem::replace(
            &mut self.multicenter_bonds,
            VarSetStorage::Shared(Arc::new(VarRelationSet::default())),
        );
        let (multicenter, removed_multicenter) = multicenter.compact(&compaction);
        self.multicenter_bonds = multicenter;

        let noncovalent = mem::replace(
            &mut self.noncovalent_bonds,
            FixedSetStorage::Shared(Arc::new(FixedRelationSet::default())),
        );
        let (noncovalent, removed_noncovalent) = noncovalent.compact(&compaction);
        self.noncovalent_bonds = noncovalent;

        let stereo_atoms = mem::replace(
            &mut self.stereo_atoms,
            FixedVarSetStorage::Shared(Arc::new(FixedVarBirelationSet::default())),
        );
        let (stereo_atoms, removed_stereo_atoms) = stereo_atoms.compact(&compaction);
        self.stereo_atoms = stereo_atoms;

        let stereo_bonds = mem::replace(
            &mut self.stereo_bonds,
            FixedVarSetStorage::Shared(Arc::new(FixedVarBirelationSet::default())),
        );
        let (stereo_bonds, removed_stereo_bonds) = stereo_bonds.compact(&compaction);
        self.stereo_bonds = stereo_bonds;

        let id_compaction = MoleculeCompaction::new(
            compaction,
            Compaction::new(
                removed_dative.source_count(),
                removed_dative
                    .removed()
                    .iter()
                    .map(|&id| DativeBondId::from(id))
                    .collect(),
            )
            .expect("relation compaction preserves its source count"),
            Compaction::new(
                removed_aromatic.source_count(),
                removed_aromatic
                    .removed()
                    .iter()
                    .map(|&id| AromaticSystemId::from(id))
                    .collect(),
            )
            .expect("relation compaction preserves its source count"),
            Compaction::new(
                removed_multicenter.source_count(),
                removed_multicenter
                    .removed()
                    .iter()
                    .map(|&id| MulticenterBondId::from(id))
                    .collect(),
            )
            .expect("relation compaction preserves its source count"),
            Compaction::new(
                removed_noncovalent.source_count(),
                removed_noncovalent
                    .removed()
                    .iter()
                    .map(|&id| NoncovalentBondId::from(id))
                    .collect(),
            )
            .expect("relation compaction preserves its source count"),
            Compaction::new(
                removed_stereo_atoms.source_count(),
                removed_stereo_atoms
                    .removed()
                    .iter()
                    .map(|&id| StereoAtomId::from(id))
                    .collect(),
            )
            .expect("relation compaction preserves its source count"),
            Compaction::new(
                removed_stereo_bonds.source_count(),
                removed_stereo_bonds
                    .removed()
                    .iter()
                    .map(|&id| StereoBondId::from(id))
                    .collect(),
            )
            .expect("relation compaction preserves its source count"),
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
    ) {
        self.restore_atoms(atoms, undo_compaction);
        self.restore_bonds(bonds, undo_compaction);
        self.restore_dative_bonds(overlays.dative_bonds, undo_compaction);
        self.restore_aromatic_systems(overlays.aromatic_systems, undo_compaction);
        self.restore_multicenter_bonds(overlays.multicenter_bonds, undo_compaction);
        self.restore_noncovalent_bonds(overlays.noncovalent_bonds, undo_compaction);
        self.restore_stereo_atoms(overlays.stereo_atoms, undo_compaction);
        self.restore_stereo_bonds(overlays.stereo_bonds, undo_compaction);
    }

    // -- Restore primitives ---------------------------------------------------

    fn restore_atoms(&mut self, removed: Vec<RemovedAtom>, undo_compaction: &UndoCompaction) {
        let mut next = vec![None; self.atoms.len() + removed.len()];
        for removed in removed {
            next[removed.id.index()] = Some(removed.attributes);
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
        let mut old_bonds: Vec<Option<BondForm>> = vec![None; self.bonds.len() + removed.len()];
        for removed in removed {
            old_endpoints[removed.id.index()] = Some(removed.endpoints);
            old_bonds[removed.id.index()] = Some(removed.attributes);
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
                removed.attributes,
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
                removed.attributes,
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
                removed.attributes,
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
                removed.attributes,
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
            next[removed.id.index()] = Some((
                [NodeId::from(removed.site)],
                removed.ligands,
                removed.attributes,
            ));
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
            next[removed.id.index()] = Some((
                [EdgeId::from(removed.site)],
                removed.ligands,
                removed.attributes,
            ));
        }
        self.stereo_bonds =
            FixedVarSetStorage::Mutable(next.into_iter().map(Option::unwrap).collect());
    }

    /// Materialize the editor's current state without consuming it, after checking molecule
    /// integrity.
    ///
    /// Subsequent editor changes are independent of the returned immutable snapshot.
    ///
    /// # Errors
    ///
    /// Returns [`MoleculeIntegrityError`] when the transient editor state cannot be published as a
    /// molecule.
    pub fn snapshot(&self) -> Result<Molecule, MoleculeIntegrityError> {
        self.clone().try_build()
    }

    /// Publish the editor's current state after checking molecule integrity.
    pub fn try_build(self) -> Result<Molecule, MoleculeIntegrityError> {
        Molecule::try_from_arcs(
            self.graph,
            self.atoms,
            self.bonds,
            DativeBonds::from_arc(self.dative_bonds.into_arc()),
            AromaticSystems::from_arc(self.aromatic_systems.into_arc()),
            MulticenterBonds::from_arc(self.multicenter_bonds.into_arc()),
            NoncovalentBonds::from_arc(self.noncovalent_bonds.into_arc()),
            StereoAtoms::from_arc(self.stereo_atoms.into_arc()),
            StereoBonds::from_arc(self.stereo_bonds.into_arc()),
            self.constraints,
        )
    }

    /// Publish editor state whose molecule integrity is established by the producer.
    ///
    /// # Panics
    ///
    /// Panics when the editor does not contain a representation-integral molecule. Use
    /// [`Self::try_build`] for independently assembled or potentially conflicting edits.
    pub fn build(self) -> Molecule {
        self.try_build()
            .unwrap_or_else(|error| panic!("invalid molecule editor state: {error}"))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

    use rstest::*;
    use umol_chem::element::Element;
    use umol_perm::Permutation;

    use super::*;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::dative::DativeBondForm;
    use crate::ir::ligand::StereoLigandKind;
    use crate::ir::noncovalent::NoncovalentBondKind;
    use crate::ir::stereo::StereoKind;
    use crate::mol_dsl;

    #[derive(Debug)]
    struct CloneCounted {
        count: Arc<AtomicUsize>,
    }

    impl Clone for CloneCounted {
        fn clone(&self) -> Self {
            self.count.fetch_add(1, AtomicOrdering::Relaxed);
            Self {
                count: Arc::clone(&self.count),
            }
        }
    }

    #[rstest]
    #[case::unique(false, 0)]
    #[case::shared(true, 1)]
    fn test_fixed_set_storage_materialize(#[case] shared: bool, #[case] expected_clones: usize) {
        let count = Arc::new(AtomicUsize::new(0));
        let relation_set: Arc<FixedRelationSet<NodeId, CloneCounted, 2>> =
            Arc::new(FixedRelationSet::new(vec![(
                [NodeId(0), NodeId(1)],
                CloneCounted {
                    count: Arc::clone(&count),
                },
            )]));
        let _shared = shared.then(|| Arc::clone(&relation_set));
        let mut storage = FixedSetStorage::Shared(relation_set);

        storage.materialize();

        assert_eq!(count.load(AtomicOrdering::Relaxed), expected_clones);
    }

    #[rstest]
    #[case::unique(false, 0)]
    #[case::shared(true, 1)]
    fn test_var_set_storage_materialize(#[case] shared: bool, #[case] expected_clones: usize) {
        let count = Arc::new(AtomicUsize::new(0));
        let relation_set: Arc<VarRelationSet<NodeId, CloneCounted>> =
            Arc::new(VarRelationSet::new(vec![(
                vec![NodeId(0), NodeId(1)],
                CloneCounted {
                    count: Arc::clone(&count),
                },
            )]));
        let _shared = shared.then(|| Arc::clone(&relation_set));
        let mut storage = VarSetStorage::Shared(relation_set);

        storage.materialize();

        assert_eq!(count.load(AtomicOrdering::Relaxed), expected_clones);
    }

    #[rstest]
    #[case::unique(false, 0)]
    #[case::shared(true, 1)]
    fn test_fixed_var_set_storage_materialize(
        #[case] shared: bool,
        #[case] expected_clones: usize,
    ) {
        let count = Arc::new(AtomicUsize::new(0));
        let relation_set: Arc<FixedVarBirelationSet<NodeId, 1, NodeId, CloneCounted>> =
            Arc::new(FixedVarBirelationSet::new(vec![(
                [NodeId(0)],
                vec![NodeId(1), NodeId(2)],
                CloneCounted {
                    count: Arc::clone(&count),
                },
            )]));
        let _shared = shared.then(|| Arc::clone(&relation_set));
        let mut storage = FixedVarSetStorage::Shared(relation_set);

        storage.materialize();

        assert_eq!(count.load(AtomicOrdering::Relaxed), expected_clones);
    }

    #[fixture]
    fn triatomic() -> MoleculeEditor {
        let mut b = Molecule::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_atom(AtomForm::from_element(Element::N));
        b.add_atom(AtomForm::from_element(Element::O));
        b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
        b.add_bond(AtomId(1), AtomId(2), BondForm::from_order(2));
        b
    }

    /// Aromatic systems, where the alignment is genuinely used: `on_permutation` reindexes the
    /// electron counts and `is_permutation_invariant` is false for a determinate vector.
    ///
    /// Classes rather than generated inputs — the methods are crate-visible and the property target
    /// cannot reach them.
    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame(vec![AtomId(0), AtomId(1), AtomId(2)], vec![10, 20, 30], true)]
    #[case::reordered_frame_carrying_its_counts(vec![AtomId(2), AtomId(0), AtomId(1)], vec![30, 10, 20], true)]
    #[case::reordered_frame_keeping_its_counts(vec![AtomId(2), AtomId(0), AtomId(1)], vec![10, 20, 30], false)]
    #[case::different_counts(vec![AtomId(0), AtomId(1), AtomId(2)], vec![10, 20, 99], false)]
    #[case::multiset_differs(vec![AtomId(0), AtomId(1), AtomId(3)], vec![10, 20, 30], false)]
    #[case::wrong_arity(vec![AtomId(0), AtomId(1)], vec![10, 20], false)]
    fn test_molecule_editor_aromatic_system_equiv(
        #[case] atoms: Vec<AtomId>,
        #[case] electrons: Vec<i64>,
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..4 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        editor.add_aromatic_system(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![10, 20, 30]),
        );
        let offered = AromaticSystemForm::from_electrons(electrons);

        assert_eq!(editor.aromatic_system_equiv(AromaticSystemId(0), &atoms, &offered), expected);
    }

    /// Multicenter bonds mirror aromatic systems: one frame-bearing factor with a position-indexed
    /// electron vector, so the alignment is used and the two readings agree.
    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame(vec![AtomId(0), AtomId(1), AtomId(2)], vec![10, 20, 30], true)]
    #[case::reordered_frame_carrying_its_counts(vec![AtomId(2), AtomId(0), AtomId(1)], vec![30, 10, 20], true)]
    #[case::reordered_frame_keeping_its_counts(vec![AtomId(2), AtomId(0), AtomId(1)], vec![10, 20, 30], false)]
    #[case::multiset_differs(vec![AtomId(0), AtomId(1), AtomId(3)], vec![10, 20, 30], false)]
    fn test_molecule_editor_multicenter_bond_equiv(
        #[case] atoms: Vec<AtomId>,
        #[case] electrons: Vec<i64>,
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..4 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        editor.add_multicenter_bond(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::from_electrons(vec![10, 20, 30]),
        );
        let offered = MulticenterBondForm::from_electrons(electrons);

        assert_eq!(editor.multicenter_bond_equiv(MulticenterBondId(0), &atoms, &offered), expected);
    }

    /// Noncovalent bonds and dative bonds carry frame-invariant payloads, so the alignment cannot
    /// change the answer and the two readings agree on identity of participants alone.
    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame([AtomId(0), AtomId(1)], true)]
    #[case::reversed_frame([AtomId(1), AtomId(0)], true)]
    #[case::different_pair([AtomId(0), AtomId(2)], false)]
    fn test_molecule_editor_noncovalent_bond_equiv(
        #[case] atoms: [AtomId; 2],
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..3 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        editor.add_noncovalent_bond(
            [AtomId(0), AtomId(1)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        let offered = NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond);

        assert_eq!(editor.noncovalent_bond_equiv(NoncovalentBondId(0), atoms, &offered), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame(AtomId(0), vec![AtomId(1), AtomId(2)], true)]
    #[case::reordered_donors(AtomId(0), vec![AtomId(2), AtomId(1)], true)]
    #[case::different_acceptor(AtomId(1), vec![AtomId(1), AtomId(2)], false)]
    #[case::different_donors(AtomId(0), vec![AtomId(1), AtomId(3)], false)]
    fn test_molecule_editor_dative_bond_equiv(
        #[case] acceptor: AtomId,
        #[case] donors: Vec<AtomId>,
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..4 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        editor.add_dative_bond(
            vec![AtomId(1), AtomId(2)],
            AtomId(0),
            DativeBondForm::from_order(1),
        );
        let offered = DativeBondForm::from_order(1);

        assert_eq!(editor.dative_bond_equiv(DativeBondId(0), acceptor, &donors, &offered), expected);
    }

    /// The editor's structural equality must reject a different atom set even when the payload
    /// carries through unread — which is what an undetermined electron vector does.
    ///
    /// Without deriving the participant action, two undetermined forms compare equal even when
    /// the offered and stored atom sets differ.
    #[rustfmt::skip]
    #[rstest]
    #[case::stored_atoms(vec![AtomId(0), AtomId(1), AtomId(2)], true)]
    #[case::reordered_atoms(vec![AtomId(2), AtomId(0), AtomId(1)], true)]
    #[case::different_atoms(vec![AtomId(0), AtomId(1), AtomId(3)], false)]
    #[case::wrong_arity(vec![AtomId(0), AtomId(1)], false)]
    fn test_molecule_editor_aromatic_system_equiv_undetermined_electrons(
        #[case] atoms: Vec<AtomId>,
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..4 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        editor.add_aromatic_system(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::default(),
        );

        assert_eq!(
            editor.aromatic_system_equiv(
                AromaticSystemId(0),
                &atoms,
                &AromaticSystemForm::default(),
            ),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored_atoms(vec![AtomId(0), AtomId(1), AtomId(2)], true)]
    #[case::different_atoms(vec![AtomId(0), AtomId(1), AtomId(3)], false)]
    fn test_molecule_editor_multicenter_bond_equiv_undetermined_electrons(
        #[case] atoms: Vec<AtomId>,
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..4 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        editor.add_multicenter_bond(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::default(),
        );

        assert_eq!(
            editor.multicenter_bond_equiv(
                MulticenterBondId(0),
                &atoms,
                &MulticenterBondForm::default(),
            ),
            expected,
        );
    }

    /// A tetrahedral centre over four distinct ligands, stored in one frame.
    #[fixture]
    fn stereo_editor() -> MoleculeEditor {
        let mut b = Molecule::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        for element in [Element::F, Element::Cl, Element::Br, Element::I] {
            b.add_atom(AtomForm::from_element(element));
        }
        for ligand in 1..=4 {
            b.add_bond(AtomId(0), AtomId(ligand), BondForm::from_order(1));
        }
        b.add_stereo_atom(
            AtomId(0),
            (1..=4)
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        );
        b
    }

    /// A coset is read against its ligand frame, so the offered configuration is restated into the
    /// stored frame before comparison and the same index under a transposed frame denotes the
    /// opposite arrangement.
    ///
    /// Classes rather than generated inputs — the methods are crate-visible and the property target
    /// cannot reach them.
    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame([1, 2, 3, 4], 0, true)]
    #[case::stored_frame_other_coset([1, 2, 3, 4], 1, false)]
    #[case::transposed_frame_same_coset([2, 1, 3, 4], 0, false)]
    #[case::transposed_frame_other_coset([2, 1, 3, 4], 1, true)]
    #[case::multiset_differs([1, 2, 3, 5], 0, false)]
    fn test_molecule_editor_stereo_atom_equiv(
        stereo_editor: MoleculeEditor,
        #[case] ligands: [u32; 4],
        #[case] coset: u32,
        #[case] expected: bool,
    ) {
        let offered: Vec<StereoLigand> = ligands
            .into_iter()
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let attributes = StereoAtomForm::new(StereoKind::Tetrahedral, coset);

        assert_eq!(
            stereo_editor.stereo_atom_equiv(StereoAtomId(0), AtomId(0), &offered, &attributes),
            expected,
        );
    }

    /// The stereo equality check must transport the configuration into the stored ligand frame.
    ///
    /// A coset is read against a frame, so the same index under a swapped frame denotes the
    /// opposite arrangement. Presenting the stored entry's own configuration against a transposed
    /// frame therefore describes a different stereocentre, and the check must say so.
    ///
    #[rstest]
    fn test_molecule_editor_stereo_atom_equiv_reordered_frame(stereo_editor: MoleculeEditor) {
        let stored: Vec<StereoLigand> = (1..=4)
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let configuration = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);

        assert!(
            stereo_editor.stereo_atom_equiv(StereoAtomId(0), AtomId(0), &stored, &configuration),
            "the stored frame with its own configuration is equivalent to itself",
        );

        let transposed = Permutation::from_image(&[1, 0, 2, 3]);
        assert!(
            !stereo_editor.stereo_atom_equiv(
                StereoAtomId(0),
                AtomId(0),
                &transposed.act(&stored),
                &configuration,
            ),
            "coset 0 against a transposed frame is the opposite arrangement, not the stored one",
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored(BondId(0), vec![2, 3, 4, 5], true)]
    #[case::within_endpoint(BondId(0), vec![3, 2, 4, 5], true)]
    #[case::endpoint_block_swap(BondId(0), vec![4, 5, 2, 3], true)]
    #[case::across_endpoints(BondId(0), vec![2, 4, 3, 5], false)]
    #[case::different_ligand(BondId(0), vec![2, 3, 4, 6], false)]
    #[case::different_site(BondId(1), vec![2, 3, 4, 5], false)]
    fn test_molecule_editor_stereo_bond_equiv(
        #[case] site: BondId,
        #[case] ligand_ids: Vec<u32>,
        #[case] expected: bool,
    ) {
        let mut editor = Molecule::default().edit();
        for _ in 0..7 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        for (first, second) in [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)] {
            editor.add_bond(AtomId(first), AtomId(second), BondForm::from_order(1));
        }
        editor.add_stereo_bond(
            BondId(0),
            [2, 3, 4, 5]
                .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
                .to_vec(),
            StereoBondForm::default(),
        );
        let ligands = ligand_ids
            .into_iter()
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect::<Vec<_>>();

        assert_eq!(
            editor.stereo_bond_equiv(
                StereoBondId(0),
                site,
                &ligands,
                &StereoBondForm::default(),
            ),
            expected,
        );
    }

    #[rstest]
    fn test_molecule_editor_restore_topology(mut triatomic: MoleculeEditor) {
        let expected = triatomic.clone().build();
        let removed_atoms = vec![RemovedAtom {
            id: AtomId(1),
            attributes: triatomic.atom(AtomId(1)).attributes.clone(),
        }];
        let removed_bonds = vec![
            RemovedBond {
                id: BondId(0),
                endpoints: triatomic.bond(BondId(0)).atoms,
                attributes: triatomic.bond(BondId(0)).attributes.clone(),
            },
            RemovedBond {
                id: BondId(1),
                endpoints: triatomic.bond(BondId(1)).atoms,
                attributes: triatomic.bond(BondId(1)).attributes.clone(),
            },
        ];

        let compaction = triatomic.tracked_remove(&[AtomId(1)], &[]);
        triatomic.restore_topology(
            removed_atoms,
            removed_bonds,
            RemovedOverlays::default(),
            &compaction.undo_compaction(),
        );

        assert_eq!(triatomic.build(), expected);
    }

    #[rstest]
    fn test_molecule_editor_remove_added_topology(mut triatomic: MoleculeEditor) {
        let expected = triatomic.clone().build();
        let added_atom = AddedAtom {
            id: triatomic.add_atom(AtomForm::from_element(Element::F)),
            attributes: AtomForm::from_element(Element::F),
        };
        let added_bond = AddedBond {
            id: triatomic.add_bond(AtomId(2), added_atom.id, BondForm::from_order(1)),
            endpoints: [AtomId(2), added_atom.id],
            attributes: BondForm::from_order(1),
        };

        triatomic.remove_added_topology(&[added_atom], &[added_bond]);

        assert_eq!(triatomic.build(), expected);
    }

    #[rstest]
    fn test_molecule_editor_restore_dative_bond() {
        let mut b = Molecule::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_atom(AtomForm::from_element(Element::N));
        b.add_dative_bond(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1));
        b.add_dative_bond(vec![AtomId(1)], AtomId(0), DativeBondForm::from_order(2));
        let expected = b.clone().build();
        let view = b.dative_bond(DativeBondId(0));
        let removed = RemovedDativeBond {
            id: DativeBondId(0),
            atoms: view.atom_ids().collect(),
            attributes: view.attributes.clone(),
        };

        b.remove_dative_bonds(&[DativeBondId(0)]);
        let undo = MoleculeCompaction::new(
            GraphCompaction::new(Compaction::identity(2), Compaction::empty()),
            Compaction::new(2, vec![removed.id])
                .expect("removed entities belong to the source table"),
            Compaction::empty(),
            Compaction::empty(),
            Compaction::empty(),
            Compaction::empty(),
            Compaction::empty(),
        )
        .undo_compaction();
        b.restore_dative_bonds(vec![removed], &undo);

        assert_eq!(b.build(), expected);
    }

    #[rstest]
    fn test_molecule_editor_snapshot(mut triatomic: MoleculeEditor) {
        let snapshot = triatomic
            .snapshot()
            .expect("the editor contains an integral molecule");
        triatomic.add_atom(AtomForm::from_element(Element::F));

        assert_eq!(
            snapshot,
            mol_dsl!(r#"{:atoms ["C" "N" "O"] :bonds [[0 1 "1"] [1 2 "2"]]}"#)
        );
        assert_eq!(
            triatomic.build(),
            mol_dsl!(r#"{:atoms ["C" "N" "O" "F"] :bonds [[0 1 "1"] [1 2 "2"]]}"#)
        );
    }

    #[rstest]
    #[case::parallel_bond(MoleculeIntegrityError::BondsParallel {
        atoms: [AtomId(0), AtomId(1)],
    })]
    fn test_molecule_editor_snapshot_error(
        #[from(triatomic)] mut editor: MoleculeEditor,
        #[case] expected: MoleculeIntegrityError,
    ) {
        editor.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));

        assert_eq!(editor.snapshot(), Err(expected));
    }

    #[rstest]
    #[case::parallel_bond(AtomId(0), AtomId(1), MoleculeIntegrityError::BondsParallel {
        atoms: [AtomId(0), AtomId(1)],
    })]
    fn test_molecule_editor_try_build_error(
        #[from(triatomic)] mut editor: MoleculeEditor,
        #[case] first: AtomId,
        #[case] second: AtomId,
        #[case] expected: MoleculeIntegrityError,
    ) {
        editor.add_bond(first, second, BondForm::from_order(1));

        assert_eq!(editor.try_build(), Err(expected));
    }

    // `edit()` → `build()` reproduces the molecule including both stereo overlays.
    #[rstest]
    fn test_molecule_editor_build() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "F" "Cl"]
                :bonds [[0 1 "1"] [1 2 "2"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 3 4 [:h 0]] :attrs "Th1"}]
                :stereo-bonds [{:site 1 :ligands [0 [:h 1] [:h 2] [:lp 2]] :attrs "Ct1"}]}"#
        );
        assert_eq!(molecule.edit().build(), molecule);
    }

    // `remove` forward-compacts stereo-atom node refs: removing a non-participant
    // shifts the surviving site/ligand ids; removing the site drops the element.
    #[rstest]
    #[case::remaps_surviving(vec![AtomId(0)], vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]])]
    #[case::drops_on_site_removal(vec![AtomId(1)], vec![])]
    fn test_molecule_editor_remove_stereo_atom(
        #[case] remove_atoms: Vec<AtomId>,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "F" "Cl" "Br"]
                :bonds [[1 2 "1"] [1 3 "1"] [1 4 "1"]]
                :stereo-atoms [{:site 1 :ligands [2 3 4 [:h 1]] :attrs "Th1"}]}"#
        );
        let mut editor = molecule.edit();
        editor.remove(&remove_atoms, &[]);
        let surviving: Vec<Vec<AtomId>> = editor
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
    fn test_molecule_editor_remove_stereo_bond(
        #[case] remove_bonds: Vec<BondId>,
        #[case] expected: Vec<BondId>,
    ) {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C"]
                :bonds [[4 5 "1"] [1 2 "2"] [0 1 "1"] [2 3 "1"]]
                :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
        );
        let mut editor = molecule.edit();
        editor.remove(&[], &remove_bonds);
        let surviving: Vec<BondId> = editor
            .build()
            .stereo_bonds()
            .iter()
            .map(|view| view.site_id())
            .collect();
        assert_eq!(surviving, expected);
    }
}
