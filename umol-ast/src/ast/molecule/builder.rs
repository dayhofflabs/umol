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
    EdgeId, FixedRelationSet, Graph, NodeId, Remapping, RelationId, VarRelationSet,
};

use super::super::aromatic::AromaticSystemAst;
use super::super::atom::AtomAst;
use super::super::bond::BondAst;
use super::super::constraint::{Constraint, Constraints};
use super::super::dative::DativeBondAst;
use super::super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::super::multicenter::MulticenterBondAst;
use super::super::noncovalent::NoncovalentBondAst;
use super::super::remap::IdRemapping;
use super::super::views::{
    AromaticSystemBuilderView, AromaticSystemBuilderViewMut, AtomBuilderView, AtomBuilderViewMut,
    BondBuilderView, BondBuilderViewMut, DativeBondBuilderView, DativeBondBuilderViewMut,
    MulticenterBondBuilderView, MulticenterBondBuilderViewMut, NoncovalentBondBuilderView,
    NoncovalentBondBuilderViewMut,
};
use super::MoleculeAst;

#[derive(Clone)]
enum FixedSetStorage<R, const N: usize> {
    Shared(Arc<FixedRelationSet<R, N>>),
    Mutable(Vec<([NodeId; N], R)>),
}

impl<R: Clone, const N: usize> FixedSetStorage<R, N> {
    fn push(&mut self, parts: [NodeId; N], data: R) -> u32 {
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
            let entries: Vec<([NodeId; N], R)> = (0..arc.relation_count())
                .map(|i| {
                    let rid = RelationId(i as u32);
                    (*arc.participants(rid), arc.data(rid).clone())
                })
                .collect();
            *self = FixedSetStorage::Mutable(entries);
        }
    }

    fn into_arc(self) -> Arc<FixedRelationSet<R, N>> {
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
                FixedSetStorage::Shared(Arc::new(remap.apply_to_fixed_relation_set(&arc)))
            }
            FixedSetStorage::Mutable(vec) => {
                let remapped: Vec<([NodeId; N], R)> = vec
                    .into_iter()
                    .filter_map(|(parts, d)| {
                        let mut new_parts = [NodeId(0); N];
                        for (i, p) in parts.iter().enumerate() {
                            new_parts[i] = remap.node(*p)?;
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
}

#[derive(Clone)]
enum VarSetStorage<R> {
    Shared(Arc<VarRelationSet<R>>),
    Mutable(Vec<(Vec<NodeId>, R)>),
}

impl<R: Clone> VarSetStorage<R> {
    fn push(&mut self, atoms: Vec<NodeId>, data: R) -> u32 {
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
            let entries: Vec<(Vec<NodeId>, R)> = (0..arc.relation_count())
                .map(|i| {
                    let rid = RelationId(i as u32);
                    (arc.participants(rid).to_vec(), arc.data(rid).clone())
                })
                .collect();
            *self = VarSetStorage::Mutable(entries);
        }
    }

    fn into_arc(self) -> Arc<VarRelationSet<R>> {
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
                VarSetStorage::Shared(Arc::new(remap.apply_to_var_relation_set(&arc)))
            }
            VarSetStorage::Mutable(vec) => {
                let remapped: Vec<(Vec<NodeId>, R)> = vec
                    .into_iter()
                    .filter_map(|(atoms, d)| {
                        let mapped: Option<Vec<NodeId>> =
                            atoms.into_iter().map(|n| remap.node(n)).collect();
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
}

fn fixed_relation_removed<R: Clone, const N: usize>(
    storage: &FixedSetStorage<R, N>,
    remap: &Remapping,
) -> Vec<u32> {
    let mut removed = Vec::new();
    for i in 0..storage.relation_count() {
        let parts = storage.participants(i);
        if parts.iter().any(|&p| remap.node(p).is_none()) {
            removed.push(i as u32);
        }
    }
    removed
}

fn var_relation_removed<R: Clone>(storage: &VarSetStorage<R>, remap: &Remapping) -> Vec<u32> {
    let mut removed = Vec::new();
    for i in 0..storage.relation_count() {
        let parts = storage.participants(i);
        if parts.iter().any(|&p| remap.node(p).is_none()) {
            removed.push(i as u32);
        }
    }
    removed
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
            dative_bonds: VarSetStorage::Shared(dative_bonds),
            aromatic_systems: VarSetStorage::Shared(aromatic_systems),
            multicenter_bonds: VarSetStorage::Shared(multicenter_bonds),
            noncovalent_bonds: FixedSetStorage::Shared(noncovalent_bonds),
            constraints,
        }
    }

    pub fn add_atom(&mut self, atom: AtomAst) -> AtomId {
        let id = self.graph.add_node();
        Arc::make_mut(&mut self.atoms).push(atom);
        AtomId::from(id)
    }

    pub fn add_bond(&mut self, src: AtomId, tgt: AtomId, bond: BondAst) -> BondId {
        let id = self.graph.add_edge(NodeId::from(src), NodeId::from(tgt));
        Arc::make_mut(&mut self.bonds).push(bond);
        BondId::from(id)
    }

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

    pub fn add_aromatic_system(
        &mut self,
        atoms: Vec<AtomId>,
        data: AromaticSystemAst,
    ) -> AromaticSystemId {
        let nodes: Vec<NodeId> = atoms.into_iter().map(NodeId::from).collect();
        let i = self.aromatic_systems.push(nodes, data);
        AromaticSystemId(i)
    }

    pub fn add_multicenter_bond(
        &mut self,
        atoms: Vec<AtomId>,
        data: MulticenterBondAst,
    ) -> MulticenterBondId {
        let nodes: Vec<NodeId> = atoms.into_iter().map(NodeId::from).collect();
        let i = self.multicenter_bonds.push(nodes, data);
        MulticenterBondId(i)
    }

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

    pub fn atom(&self, idx: AtomId) -> AtomBuilderView<'_> {
        AtomBuilderView { id: idx, ast: &self.atoms[idx.index()] }
    }

    pub fn atom_mut(&mut self, idx: AtomId) -> AtomBuilderViewMut<'_> {
        let ast = &mut Arc::make_mut(&mut self.atoms)[idx.index()];
        AtomBuilderViewMut { id: idx, ast }
    }

    pub fn bond(&self, idx: BondId) -> BondBuilderView<'_> {
        let endpoints = self.graph.edge_endpoints(EdgeId::from(idx));
        let atoms = [AtomId::from(endpoints[0]), AtomId::from(endpoints[1])];
        BondBuilderView { id: idx, ast: &self.bonds[idx.index()], atoms }
    }

    pub fn bond_mut(&mut self, idx: BondId) -> BondBuilderViewMut<'_> {
        let endpoints = self.graph.edge_endpoints(EdgeId::from(idx));
        let atoms = [AtomId::from(endpoints[0]), AtomId::from(endpoints[1])];
        let ast = &mut Arc::make_mut(&mut self.bonds)[idx.index()];
        BondBuilderViewMut { id: idx, ast, atoms }
    }

    pub fn dative_bond(&self, idx: DativeBondId) -> DativeBondBuilderView<'_> {
        match &self.dative_bonds {
            VarSetStorage::Shared(arc) => {
                let rid = RelationId(idx.0);
                let atoms = arc.participants(rid);
                let ast = arc.data(rid);
                let acceptor_id = AtomId::from(atoms[ast.acceptor_slot as usize]);
                DativeBondBuilderView { id: idx, ast, atoms, acceptor_id }
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
                AromaticSystemBuilderView { id: idx, ast: &entry.1, atoms: &entry.0 }
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
                MulticenterBondBuilderView { id: idx, ast: &entry.1, atoms: &entry.0 }
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

    pub fn remove_dative_bonds(&mut self, indices: &[DativeBondId]) {
        let raw: Vec<u32> = indices.iter().map(|i| i.0).collect();
        self.dative_bonds.remove_indices(&raw);
        let idx_remap = IdRemapping::new(
            Remapping {
                removed_nodes: Vec::new(),
                removed_edges: Vec::new(),
            },
            raw,
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        self.constraints.remap(&idx_remap);
    }

    pub fn remove_aromatic_systems(&mut self, indices: &[AromaticSystemId]) {
        let raw: Vec<u32> = indices.iter().map(|i| i.0).collect();
        self.aromatic_systems.remove_indices(&raw);
        let idx_remap = IdRemapping::new(
            Remapping {
                removed_nodes: Vec::new(),
                removed_edges: Vec::new(),
            },
            Vec::new(),
            raw,
            Vec::new(),
            Vec::new(),
        );
        self.constraints.remap(&idx_remap);
    }

    pub fn remove_multicenter_bonds(&mut self, indices: &[MulticenterBondId]) {
        let raw: Vec<u32> = indices.iter().map(|i| i.0).collect();
        self.multicenter_bonds.remove_indices(&raw);
        let idx_remap = IdRemapping::new(
            Remapping {
                removed_nodes: Vec::new(),
                removed_edges: Vec::new(),
            },
            Vec::new(),
            Vec::new(),
            raw,
            Vec::new(),
        );
        self.constraints.remap(&idx_remap);
    }

    pub fn remove_noncovalent_bonds(&mut self, indices: &[NoncovalentBondId]) {
        let raw: Vec<u32> = indices.iter().map(|i| i.0).collect();
        self.noncovalent_bonds.remove_indices(&raw);
        let idx_remap = IdRemapping::new(
            Remapping {
                removed_nodes: Vec::new(),
                removed_edges: Vec::new(),
            },
            Vec::new(),
            Vec::new(),
            Vec::new(),
            raw,
        );
        self.constraints.remap(&idx_remap);
    }

    // -- Topological removal --------------------------------------------------

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
