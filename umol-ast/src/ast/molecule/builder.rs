//! Structural editing for `MoleculeAst`. The AST itself only allows attribute
//! mutation; structural change (add atoms/bonds/relations, remove anything)
//! goes through `MoleculeBuilder`.
//!
//! Storage is lazy: each Arc-shared field stays shared until first write,
//! at which point only that field decomposes to a mutable form. `build`
//! re-wraps everything in `Arc`, reusing untouched shared data.

use std::collections::HashSet;
use std::sync::Arc;
use std::{iter, mem};

use umol_graph_core::{
    EdgeId, FixedRelationSet, Graph, NodeId, RelationId, RelationParticipant, Remapping, Unordered,
    VarRelationSet,
};

use super::super::aromatic::AromaticSystemAst;
use super::super::atom::AtomAst;
use super::super::bond::BondAst;
use super::super::constraint::{Constraint, Constraints};
use super::super::dative::DativeBondAst;
use super::super::edit::{
    AddedAromaticSystem, AddedAtom, AddedBond, AddedDativeBond, AddedMulticenterBond,
    AddedNoncovalentBond, ConstraintUpdate, RemovedAromaticSystem, RemovedAtom, RemovedBond,
    RemovedDativeBond, RemovedMulticenterBond, RemovedNoncovalentBond, RemovedOverlays,
};
use super::super::ids::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::super::multicenter::MulticenterBondAst;
use super::super::noncovalent::NoncovalentBondAst;
use super::super::remap::{IdRemapping, UndoRemapping};
use super::super::views::{
    AromaticSystemBuilderView, AromaticSystemBuilderViewMut, AtomBuilderView, AtomBuilderViewMut,
    BondBuilderView, BondBuilderViewMut, DativeBondBuilderView, DativeBondBuilderViewMut,
    MulticenterBondBuilderView, MulticenterBondBuilderViewMut, NoncovalentBondBuilderView,
    NoncovalentBondBuilderViewMut,
};
use super::MoleculeAst;

#[derive(Clone)]
enum FixedSetStorage<D, const N: usize> {
    Shared(Arc<FixedRelationSet<NodeId, Unordered, D, N>>),
    Mutable(Vec<([NodeId; N], D)>),
}

impl<D: Clone, const N: usize> FixedSetStorage<D, N> {
    fn push(&mut self, parts: [NodeId; N], data: D) -> u32 {
        self.materialize();
        let FixedSetStorage::Mutable(vec) = self else {
            unreachable!()
        };
        let idx = vec.len() as u32;
        vec.push((parts, data));
        idx
    }

    fn materialize(&mut self) {
        if let FixedSetStorage::Shared(arc) = self {
            let entries: Vec<([NodeId; N], D)> = (0..arc.relation_count())
                .map(|i| {
                    let rid = RelationId(i as u32);
                    (*arc.participants(rid), arc.data(rid).clone())
                })
                .collect();
            *self = FixedSetStorage::Mutable(entries);
        }
    }

    fn into_arc(self) -> Arc<FixedRelationSet<NodeId, Unordered, D, N>> {
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

    fn participants(&self, i: usize) -> [NodeId; N] {
        match self {
            FixedSetStorage::Shared(arc) => *arc.participants(RelationId(i as u32)),
            FixedSetStorage::Mutable(vec) => vec[i].0,
        }
    }

    fn apply_remapping(self, remap: &Remapping) -> Self {
        match self {
            FixedSetStorage::Shared(arc) => {
                FixedSetStorage::Shared(Arc::new(arc.apply_remapping(remap)))
            }
            FixedSetStorage::Mutable(vec) => {
                let remapped: Vec<([NodeId; N], D)> = vec
                    .into_iter()
                    .filter_map(|(parts, d)| {
                        let mut new_parts = [NodeId(0); N];
                        for (i, p) in parts.iter().enumerate() {
                            new_parts[i] = remap.map_node(*p)?;
                        }
                        Some((new_parts, d))
                    })
                    .collect();
                FixedSetStorage::Mutable(remapped)
            }
        }
    }

    fn remove_indices(&mut self, indices: &[u32]) {
        if indices.is_empty() {
            return;
        }
        self.materialize();
        let FixedSetStorage::Mutable(vec) = self else {
            unreachable!()
        };
        let remove: HashSet<u32> = indices.iter().copied().collect();
        let mut dst = 0usize;
        for src in 0..vec.len() {
            if !remove.contains(&(src as u32)) {
                vec.swap(dst, src);
                dst += 1;
            }
        }
        vec.truncate(dst);
    }

    fn entries(&self) -> Vec<([NodeId; N], D)> {
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
enum VarSetStorage<D> {
    Shared(Arc<VarRelationSet<NodeId, Unordered, D>>),
    Mutable(Vec<(Vec<NodeId>, D)>),
}

impl<D: Clone> VarSetStorage<D> {
    fn push(&mut self, atoms: Vec<NodeId>, data: D) -> u32 {
        self.materialize();
        let VarSetStorage::Mutable(vec) = self else {
            unreachable!()
        };
        let idx = vec.len() as u32;
        vec.push((atoms, data));
        idx
    }

    fn materialize(&mut self) {
        if let VarSetStorage::Shared(arc) = self {
            let entries: Vec<(Vec<NodeId>, D)> = (0..arc.relation_count())
                .map(|i| {
                    let rid = RelationId(i as u32);
                    (arc.participants(rid).to_vec(), arc.data(rid).clone())
                })
                .collect();
            *self = VarSetStorage::Mutable(entries);
        }
    }

    fn into_arc(self) -> Arc<VarRelationSet<NodeId, Unordered, D>> {
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

    fn participants(&self, i: usize) -> Vec<NodeId> {
        match self {
            VarSetStorage::Shared(arc) => arc.participants(RelationId(i as u32)).to_vec(),
            VarSetStorage::Mutable(vec) => vec[i].0.clone(),
        }
    }

    fn apply_remapping(self, remap: &Remapping) -> Self {
        match self {
            VarSetStorage::Shared(arc) => {
                VarSetStorage::Shared(Arc::new(arc.apply_remapping(remap)))
            }
            VarSetStorage::Mutable(vec) => {
                let remapped: Vec<(Vec<NodeId>, D)> = vec
                    .into_iter()
                    .filter_map(|(atoms, d)| {
                        let mapped: Option<Vec<NodeId>> =
                            atoms.into_iter().map(|n| remap.map_node(n)).collect();
                        mapped.map(|a| (a, d))
                    })
                    .collect();
                VarSetStorage::Mutable(remapped)
            }
        }
    }

    fn remove_indices(&mut self, indices: &[u32]) {
        if indices.is_empty() {
            return;
        }
        self.materialize();
        let VarSetStorage::Mutable(vec) = self else {
            unreachable!()
        };
        let remove: HashSet<u32> = indices.iter().copied().collect();
        let mut dst = 0usize;
        for src in 0..vec.len() {
            if !remove.contains(&(src as u32)) {
                vec.swap(dst, src);
                dst += 1;
            }
        }
        vec.truncate(dst);
    }

    fn entries(&self) -> Vec<(Vec<NodeId>, D)> {
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

fn fixed_relation_removed<D: Clone, const N: usize>(
    storage: &FixedSetStorage<D, N>,
    remap: &Remapping,
) -> Vec<u32> {
    let mut removed = Vec::new();
    for i in 0..storage.relation_count() {
        let parts = storage.participants(i);
        if parts.iter().any(|&p| remap.map_node(p).is_none()) {
            removed.push(i as u32);
        }
    }
    removed
}

fn var_relation_removed<D: Clone>(storage: &VarSetStorage<D>, remap: &Remapping) -> Vec<u32> {
    let mut removed = Vec::new();
    for i in 0..storage.relation_count() {
        let parts = storage.participants(i);
        if parts.iter().any(|&p| remap.map_node(p).is_none()) {
            removed.push(i as u32);
        }
    }
    removed
}

fn restore_var_participants(parts: Vec<NodeId>, undo_remapping: &UndoRemapping) -> Vec<NodeId> {
    let remapping = undo_remapping.forward().graph();
    parts.into_iter().map(|p| p.unmap(remapping)).collect()
}

fn restore_fixed_participants<const N: usize>(
    parts: [NodeId; N],
    undo_remapping: &UndoRemapping,
) -> [NodeId; N] {
    let remapping = undo_remapping.forward().graph();
    parts.map(|p| p.unmap(remapping))
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
    dative_bonds: VarSetStorage<DativeBondAst>,
    aromatic_systems: VarSetStorage<AromaticSystemAst>,
    multicenter_bonds: VarSetStorage<MulticenterBondAst>,
    noncovalent_bonds: FixedSetStorage<NoncovalentBondAst, 2>,
    constraints: Constraints,
}

impl MoleculeBuilder {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn from_parts(
        graph: Graph,
        atoms: Arc<Vec<AtomAst>>,
        bonds: Arc<Vec<BondAst>>,
        dative_bonds: Arc<VarRelationSet<NodeId, Unordered, DativeBondAst>>,
        aromatic_systems: Arc<VarRelationSet<NodeId, Unordered, AromaticSystemAst>>,
        multicenter_bonds: Arc<VarRelationSet<NodeId, Unordered, MulticenterBondAst>>,
        noncovalent_bonds: Arc<FixedRelationSet<NodeId, Unordered, NoncovalentBondAst, 2>>,
        constraints: Constraints,
    ) -> Self {
        Self {
            graph,
            atoms,
            bonds,
            dative_bonds: VarSetStorage::Shared(dative_bonds),
            aromatic_systems: VarSetStorage::Shared(aromatic_systems),
            multicenter_bonds: VarSetStorage::Shared(multicenter_bonds),
            noncovalent_bonds: FixedSetStorage::Shared(noncovalent_bonds),
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
    /// assumes `src` and `tgt` are valid atom ids in the current dense layout.
    pub fn add_bond(&mut self, src: AtomId, tgt: AtomId, bond: BondAst) -> BondId {
        let id = self.graph.add_edge(NodeId::from(src), NodeId::from(tgt));
        Arc::make_mut(&mut self.bonds).push(bond);
        BondId::from(id)
    }

    /// Append a dative-bond overlay directly to the builder.
    ///
    /// The participant list is sorted into dense atom-id order and the
    /// acceptor slot is normalized into that sorted representation.
    pub fn add_dative_bond(
        &mut self,
        donors: Vec<AtomId>,
        acceptor: AtomId,
        mut bond: DativeBondAst,
    ) -> DativeBondId {
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
        bond.acceptor_slot = slot as u8;
        let i = self.dative_bonds.push(participants, bond);
        DativeBondId(i)
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

    /// Add a molecule-level constraint (molecule-scope predicate or
    /// combinator). Unconditional per-entity constraints belong inline on the
    /// entity — use `atom_mut(idx).constraints.set(c)` etc.
    pub fn push_constraint(&mut self, c: Constraint) {
        self.constraints.push(c);
    }

    // -- Attribute access -----------------------------------------------------
    //
    // Mutable views edit entity data in place. Structural add/remove stays on
    // the builder itself because dense removal can remap many unrelated ids.

    pub fn atom(&self, idx: AtomId) -> AtomBuilderView<'_> {
        AtomBuilderView {
            id: idx,
            ast: &self.atoms[idx.index()],
        }
    }

    pub fn atom_mut(&mut self, idx: AtomId) -> AtomBuilderViewMut<'_> {
        let ast = &mut Arc::make_mut(&mut self.atoms)[idx.index()];
        AtomBuilderViewMut { id: idx, ast }
    }

    pub fn bond(&self, idx: BondId) -> BondBuilderView<'_> {
        let endpoints = self.graph.edge_endpoints(EdgeId::from(idx));
        let atoms = [AtomId::from(endpoints[0]), AtomId::from(endpoints[1])];
        BondBuilderView {
            id: idx,
            ast: &self.bonds[idx.index()],
            atoms,
        }
    }

    pub fn bond_mut(&mut self, idx: BondId) -> BondBuilderViewMut<'_> {
        let endpoints = self.graph.edge_endpoints(EdgeId::from(idx));
        let atoms = [AtomId::from(endpoints[0]), AtomId::from(endpoints[1])];
        let ast = &mut Arc::make_mut(&mut self.bonds)[idx.index()];
        BondBuilderViewMut {
            id: idx,
            ast,
            atoms,
        }
    }

    pub fn dative_bond(&self, idx: DativeBondId) -> DativeBondBuilderView<'_> {
        match &self.dative_bonds {
            VarSetStorage::Shared(arc) => {
                let rid = RelationId(idx.0);
                let atoms = arc.participants(rid);
                let ast = arc.data(rid);
                let acceptor_id = AtomId::from(atoms[ast.acceptor_slot as usize]);
                DativeBondBuilderView {
                    id: idx,
                    ast,
                    atoms,
                    acceptor_id,
                }
            }
            VarSetStorage::Mutable(vec) => {
                let entry = &vec[idx.index()];
                let acceptor_id = AtomId::from(entry.0[entry.1.acceptor_slot as usize]);
                DativeBondBuilderView {
                    id: idx,
                    ast: &entry.1,
                    atoms: &entry.0,
                    acceptor_id,
                }
            }
        }
    }

    pub fn dative_bond_mut(&mut self, idx: DativeBondId) -> DativeBondBuilderViewMut<'_> {
        self.dative_bonds.materialize();
        let VarSetStorage::Mutable(vec) = &mut self.dative_bonds else {
            unreachable!()
        };
        let entry = &mut vec[idx.index()];
        let acceptor_id = AtomId::from(entry.0[entry.1.acceptor_slot as usize]);
        DativeBondBuilderViewMut {
            id: idx,
            atoms: &entry.0,
            ast: &mut entry.1,
            acceptor_id,
        }
    }

    pub fn aromatic_system(&self, idx: AromaticSystemId) -> AromaticSystemBuilderView<'_> {
        match &self.aromatic_systems {
            VarSetStorage::Shared(arc) => {
                let rid = RelationId(idx.0);
                AromaticSystemBuilderView {
                    id: idx,
                    ast: arc.data(rid),
                    atoms: arc.participants(rid),
                }
            }
            VarSetStorage::Mutable(vec) => {
                let entry = &vec[idx.index()];
                AromaticSystemBuilderView {
                    id: idx,
                    ast: &entry.1,
                    atoms: &entry.0,
                }
            }
        }
    }

    pub fn aromatic_system_mut(
        &mut self,
        idx: AromaticSystemId,
    ) -> AromaticSystemBuilderViewMut<'_> {
        self.aromatic_systems.materialize();
        let VarSetStorage::Mutable(vec) = &mut self.aromatic_systems else {
            unreachable!()
        };
        let entry = &mut vec[idx.index()];
        AromaticSystemBuilderViewMut {
            id: idx,
            atoms: &entry.0,
            ast: &mut entry.1,
        }
    }

    pub fn multicenter_bond(&self, idx: MulticenterBondId) -> MulticenterBondBuilderView<'_> {
        match &self.multicenter_bonds {
            VarSetStorage::Shared(arc) => {
                let rid = RelationId(idx.0);
                MulticenterBondBuilderView {
                    id: idx,
                    ast: arc.data(rid),
                    atoms: arc.participants(rid),
                }
            }
            VarSetStorage::Mutable(vec) => {
                let entry = &vec[idx.index()];
                MulticenterBondBuilderView {
                    id: idx,
                    ast: &entry.1,
                    atoms: &entry.0,
                }
            }
        }
    }

    pub fn multicenter_bond_mut(
        &mut self,
        idx: MulticenterBondId,
    ) -> MulticenterBondBuilderViewMut<'_> {
        self.multicenter_bonds.materialize();
        let VarSetStorage::Mutable(vec) = &mut self.multicenter_bonds else {
            unreachable!()
        };
        let entry = &mut vec[idx.index()];
        MulticenterBondBuilderViewMut {
            id: idx,
            atoms: &entry.0,
            ast: &mut entry.1,
        }
    }

    pub fn noncovalent_bond(&self, idx: NoncovalentBondId) -> NoncovalentBondBuilderView<'_> {
        match &self.noncovalent_bonds {
            FixedSetStorage::Shared(arc) => {
                let rid = RelationId(idx.0);
                let parts = arc.participants(rid);
                NoncovalentBondBuilderView {
                    id: idx,
                    ast: arc.data(rid),
                    atoms: [AtomId::from(parts[0]), AtomId::from(parts[1])],
                }
            }
            FixedSetStorage::Mutable(vec) => {
                let entry = &vec[idx.index()];
                NoncovalentBondBuilderView {
                    id: idx,
                    ast: &entry.1,
                    atoms: [AtomId::from(entry.0[0]), AtomId::from(entry.0[1])],
                }
            }
        }
    }

    pub fn noncovalent_bond_mut(
        &mut self,
        idx: NoncovalentBondId,
    ) -> NoncovalentBondBuilderViewMut<'_> {
        self.noncovalent_bonds.materialize();
        let FixedSetStorage::Mutable(vec) = &mut self.noncovalent_bonds else {
            unreachable!()
        };
        let entry = &mut vec[idx.index()];
        let atoms = [AtomId::from(entry.0[0]), AtomId::from(entry.0[1])];
        NoncovalentBondBuilderViewMut {
            id: idx,
            ast: &mut entry.1,
            atoms,
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

    // -- Relation removal -----------------------------------------------------

    /// Remove dative-bond overlays directly from the builder.
    ///
    /// This is a low-level dense removal primitive. It remaps molecule-level
    /// constraints but does not build rollback data.
    pub fn remove_dative_bonds(&mut self, indices: &[DativeBondId]) {
        let raw: Vec<u32> = indices.iter().map(|i| i.0).collect();
        self.dative_bonds.remove_indices(&raw);
        let idx_remap = IdRemapping::new(
            Remapping::new(Vec::new(), Vec::new()),
            raw,
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        self.constraints.remap(&idx_remap);
    }

    /// Remove aromatic-system overlays directly from the builder.
    ///
    /// This is a low-level dense removal primitive. It remaps molecule-level
    /// constraints but does not build rollback data.
    pub fn remove_aromatic_systems(&mut self, indices: &[AromaticSystemId]) {
        let raw: Vec<u32> = indices.iter().map(|i| i.0).collect();
        self.aromatic_systems.remove_indices(&raw);
        let idx_remap = IdRemapping::new(
            Remapping::new(Vec::new(), Vec::new()),
            Vec::new(),
            raw,
            Vec::new(),
            Vec::new(),
        );
        self.constraints.remap(&idx_remap);
    }

    /// Remove multicenter-bond overlays directly from the builder.
    ///
    /// This is a low-level dense removal primitive. It remaps molecule-level
    /// constraints but does not build rollback data.
    pub fn remove_multicenter_bonds(&mut self, indices: &[MulticenterBondId]) {
        let raw: Vec<u32> = indices.iter().map(|i| i.0).collect();
        self.multicenter_bonds.remove_indices(&raw);
        let idx_remap = IdRemapping::new(
            Remapping::new(Vec::new(), Vec::new()),
            Vec::new(),
            Vec::new(),
            raw,
            Vec::new(),
        );
        self.constraints.remap(&idx_remap);
    }

    /// Remove noncovalent-bond overlays directly from the builder.
    ///
    /// This is a low-level dense removal primitive. It remaps molecule-level
    /// constraints but does not build rollback data.
    pub fn remove_noncovalent_bonds(&mut self, indices: &[NoncovalentBondId]) {
        let raw: Vec<u32> = indices.iter().map(|i| i.0).collect();
        self.noncovalent_bonds.remove_indices(&raw);
        let idx_remap = IdRemapping::new(
            Remapping::new(Vec::new(), Vec::new()),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            raw,
        );
        self.constraints.remap(&idx_remap);
    }

    // -- Topological removal --------------------------------------------------

    /// Remove topology directly from the builder and return the forward remap.
    ///
    /// This is the low-level dense topology-removal primitive. It removes the
    /// requested atoms and bonds, cascades relations whose participants were
    /// removed, remaps molecule-level constraints, and returns the forward
    /// `IdRemapping` for downstream id holders. It does not build rollback
    /// data; checked transactions capture the removed payloads before calling
    /// this method.
    pub fn remove(&mut self, atoms: &[AtomId], bonds: &[BondId]) -> IdRemapping {
        let nodes: Vec<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        let edges: Vec<EdgeId> = bonds.iter().map(|&b| EdgeId::from(b)).collect();
        let remap = self.graph.remove(&nodes, &edges);

        let new_atoms = remap.apply_to_node_vec(&self.atoms);
        let new_bonds = remap.apply_to_edge_vec(&self.bonds);
        self.atoms = Arc::new(new_atoms);
        self.bonds = Arc::new(new_bonds);

        let removed_dative = var_relation_removed(&self.dative_bonds, &remap);
        let removed_aromatic = var_relation_removed(&self.aromatic_systems, &remap);
        let removed_multicenter = var_relation_removed(&self.multicenter_bonds, &remap);
        let removed_noncovalent = fixed_relation_removed(&self.noncovalent_bonds, &remap);

        let dative = mem::replace(
            &mut self.dative_bonds,
            VarSetStorage::Shared(Arc::new(VarRelationSet::default())),
        );
        self.dative_bonds = dative.apply_remapping(&remap);

        let aromatic = mem::replace(
            &mut self.aromatic_systems,
            VarSetStorage::Shared(Arc::new(VarRelationSet::default())),
        );
        self.aromatic_systems = aromatic.apply_remapping(&remap);

        let multicenter = mem::replace(
            &mut self.multicenter_bonds,
            VarSetStorage::Shared(Arc::new(VarRelationSet::default())),
        );
        self.multicenter_bonds = multicenter.apply_remapping(&remap);

        let noncovalent = mem::replace(
            &mut self.noncovalent_bonds,
            FixedSetStorage::Shared(Arc::new(FixedRelationSet::default())),
        );
        self.noncovalent_bonds = noncovalent.apply_remapping(&remap);

        let idx_remap = IdRemapping::new(
            remap,
            removed_dative,
            removed_aromatic,
            removed_multicenter,
            removed_noncovalent,
        );
        self.constraints.remap(&idx_remap);
        idx_remap
    }

    pub(super) fn remove_added_topology(&mut self, atoms: &[AddedAtom], bonds: &[AddedBond]) {
        let atom_ids: Vec<AtomId> = atoms.iter().map(|a| a.id).collect();
        let bond_ids: Vec<BondId> = bonds.iter().map(|b| b.id).collect();
        self.remove(&atom_ids, &bond_ids);
    }

    pub(super) fn restore_topology(
        &mut self,
        atoms: Vec<RemovedAtom>,
        bonds: Vec<RemovedBond>,
        overlays: RemovedOverlays,
        undo_remapping: &UndoRemapping,
        constraint_update: ConstraintUpdate,
    ) {
        self.restore_atoms(atoms, undo_remapping);
        self.restore_bonds(bonds, undo_remapping);
        self.restore_dative_bonds(overlays.dative_bonds, undo_remapping);
        self.restore_aromatic_systems(overlays.aromatic_systems, undo_remapping);
        self.restore_multicenter_bonds(overlays.multicenter_bonds, undo_remapping);
        self.restore_noncovalent_bonds(overlays.noncovalent_bonds, undo_remapping);
        constraint_update.rollback_into(&mut self.constraints);
    }

    pub(super) fn remove_added_dative_bond(&mut self, added: &AddedDativeBond) {
        self.remove_dative_bonds(&[added.id]);
    }

    pub(super) fn restore_dative_bond(
        &mut self,
        removed: RemovedDativeBond,
        undo_remapping: &UndoRemapping,
    ) {
        self.restore_dative_bonds(vec![removed], undo_remapping);
    }

    pub(super) fn remove_added_aromatic_system(&mut self, added: &AddedAromaticSystem) {
        self.remove_aromatic_systems(&[added.id]);
    }

    pub(super) fn restore_aromatic_system(
        &mut self,
        removed: RemovedAromaticSystem,
        undo_remapping: &UndoRemapping,
    ) {
        self.restore_aromatic_systems(vec![removed], undo_remapping);
    }

    pub(super) fn remove_added_multicenter_bond(&mut self, added: &AddedMulticenterBond) {
        self.remove_multicenter_bonds(&[added.id]);
    }

    pub(super) fn restore_multicenter_bond(
        &mut self,
        removed: RemovedMulticenterBond,
        undo_remapping: &UndoRemapping,
    ) {
        self.restore_multicenter_bonds(vec![removed], undo_remapping);
    }

    pub(super) fn remove_added_noncovalent_bond(&mut self, added: &AddedNoncovalentBond) {
        self.remove_noncovalent_bonds(&[added.id]);
    }

    pub(super) fn restore_noncovalent_bond(
        &mut self,
        removed: RemovedNoncovalentBond,
        undo_remapping: &UndoRemapping,
    ) {
        self.restore_noncovalent_bonds(vec![removed], undo_remapping);
    }

    fn restore_atoms(&mut self, removed: Vec<RemovedAtom>, undo_remapping: &UndoRemapping) {
        let mut next = vec![None; self.atoms.len() + removed.len()];
        for removed in removed {
            next[removed.id.index()] = Some(removed.ast);
        }
        for (idx, atom) in self.atoms.iter().cloned().enumerate() {
            let old = undo_remapping.atom(AtomId(idx as u32));
            next[old.index()] = Some(atom);
        }
        self.atoms = Arc::new(next.into_iter().map(Option::unwrap).collect());
    }

    fn restore_bonds(&mut self, removed: Vec<RemovedBond>, undo_remapping: &UndoRemapping) {
        let mut old_endpoints: Vec<Option<[AtomId; 2]>> =
            vec![None; self.bonds.len() + removed.len()];
        let mut old_bonds: Vec<Option<BondAst>> = vec![None; self.bonds.len() + removed.len()];
        for removed in removed {
            old_endpoints[removed.id.index()] = Some(removed.endpoints);
            old_bonds[removed.id.index()] = Some(removed.ast);
        }
        for (idx, bond) in self.bonds.iter().cloned().enumerate() {
            let old_id = undo_remapping.bond(BondId(idx as u32));
            let endpoints = self.graph.edge_endpoints(EdgeId(idx as u32));
            old_endpoints[old_id.index()] = Some([
                undo_remapping.atom(AtomId::from(endpoints[0])),
                undo_remapping.atom(AtomId::from(endpoints[1])),
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

    fn restore_dative_bonds(
        &mut self,
        removed: Vec<RemovedDativeBond>,
        undo_remapping: &UndoRemapping,
    ) {
        let current = self.dative_bonds.entries();
        let mut next = vec![None; current.len() + removed.len()];
        for (idx, (parts, data)) in current.into_iter().enumerate() {
            let old_id = undo_remapping.dative_bond(DativeBondId(idx as u32));
            next[old_id.index()] = Some((restore_var_participants(parts, undo_remapping), data));
        }
        for removed in removed {
            next[removed.id.index()] = Some((
                removed.atoms.into_iter().map(NodeId::from).collect(),
                removed.ast,
            ));
        }
        self.dative_bonds = VarSetStorage::Mutable(next.into_iter().map(Option::unwrap).collect());
    }

    fn restore_aromatic_systems(
        &mut self,
        removed: Vec<RemovedAromaticSystem>,
        undo_remapping: &UndoRemapping,
    ) {
        let current = self.aromatic_systems.entries();
        let mut next = vec![None; current.len() + removed.len()];
        for (idx, (parts, data)) in current.into_iter().enumerate() {
            let old_id = undo_remapping.aromatic_system(AromaticSystemId(idx as u32));
            next[old_id.index()] = Some((restore_var_participants(parts, undo_remapping), data));
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

    fn restore_multicenter_bonds(
        &mut self,
        removed: Vec<RemovedMulticenterBond>,
        undo_remapping: &UndoRemapping,
    ) {
        let current = self.multicenter_bonds.entries();
        let mut next = vec![None; current.len() + removed.len()];
        for (idx, (parts, data)) in current.into_iter().enumerate() {
            let old_id = undo_remapping.multicenter_bond(MulticenterBondId(idx as u32));
            next[old_id.index()] = Some((restore_var_participants(parts, undo_remapping), data));
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

    fn restore_noncovalent_bonds(
        &mut self,
        removed: Vec<RemovedNoncovalentBond>,
        undo_remapping: &UndoRemapping,
    ) {
        let current = self.noncovalent_bonds.entries();
        let mut next = vec![None; current.len() + removed.len()];
        for (idx, (parts, data)) in current.into_iter().enumerate() {
            let old_id = undo_remapping.noncovalent_bond(NoncovalentBondId(idx as u32));
            next[old_id.index()] = Some((restore_fixed_participants(parts, undo_remapping), data));
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

    pub fn build(self) -> MoleculeAst {
        MoleculeAst::from_arcs(
            self.graph,
            self.atoms,
            self.bonds,
            self.dative_bonds.into_arc(),
            self.aromatic_systems.into_arc(),
            self.multicenter_bonds.into_arc(),
            self.noncovalent_bonds.into_arc(),
            self.constraints,
        )
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::AtomConstraint;
    use crate::ast::dative::DativeBondAst;
    use crate::ast::DroppedConstraint;

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

        let remap = triatomic.remove(&[AtomId(1)], &[]);
        triatomic.restore_topology(
            removed_atoms,
            removed_bonds,
            RemovedOverlays::default(),
            &remap.undo_remapping(),
            ConstraintUpdate {
                dropped: vec![DroppedConstraint {
                    position: 0,
                    constraint: dropped_constraint,
                }],
                rewritten: Vec::new(),
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
        let undo = IdRemapping::relations(vec![removed.id.0], Vec::new(), Vec::new(), Vec::new())
            .undo_remapping();
        b.restore_dative_bond(removed, &undo);

        assert_eq!(b.build(), expected);
    }
}
