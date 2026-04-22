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

use umol_graph_core::relation::RelationId;
use umol_graph_core::{EdgeId, FixedRelationSet, Graph, NodeId, Remapping, VarRelationSet};

use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::constraint::{
    AromaticSystemConstraint, AtomConstraint, BondConstraint, Constraint, Constraints,
    DativeBondConstraint, MulticenterBondConstraint, NoncovalentBondConstraint,
};
use crate::ast::dative::{DativeBondAst, DativeDirection};
use crate::ast::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use crate::ast::multicenter::MulticenterBondAst;
use crate::ast::noncovalent::NoncovalentBondAst;
use crate::ast::remap::IdxRemapping;

use super::MoleculeAst;

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

pub struct MoleculeBuilder {
    graph: Graph,
    atoms: Arc<Vec<AtomAst>>,
    bonds: Arc<Vec<BondAst>>,
    dative_bonds: FixedSetStorage<DativeBondAst, 2>,
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
            dative_bonds: FixedSetStorage::Shared(dative_bonds),
            aromatic_systems: VarSetStorage::Shared(aromatic_systems),
            multicenter_bonds: VarSetStorage::Shared(multicenter_bonds),
            noncovalent_bonds: FixedSetStorage::Shared(noncovalent_bonds),
            constraints,
        }
    }

    pub fn add_atom(&mut self, atom: AtomAst) -> AtomIdx {
        let id = self.graph.add_node();
        Arc::make_mut(&mut self.atoms).push(atom);
        AtomIdx::from(id)
    }

    pub fn add_bond(&mut self, src: AtomIdx, tgt: AtomIdx, bond: BondAst) -> BondIdx {
        let id = self.graph.add_edge(NodeId::from(src), NodeId::from(tgt));
        Arc::make_mut(&mut self.bonds).push(bond);
        BondIdx::from(id)
    }

    pub fn add_dative_bond(
        &mut self,
        donor: AtomIdx,
        acceptor: AtomIdx,
        mut bond: DativeBondAst,
    ) -> DativeBondIdx {
        bond.direction = if donor.0 <= acceptor.0 {
            DativeDirection::Forward
        } else {
            DativeDirection::Reverse
        };
        let i = self
            .dative_bonds
            .push([NodeId::from(donor), NodeId::from(acceptor)], bond);
        DativeBondIdx(i)
    }

    pub fn add_aromatic_system(
        &mut self,
        atoms: Vec<AtomIdx>,
        data: AromaticSystemAst,
    ) -> AromaticSystemIdx {
        let nodes: Vec<NodeId> = atoms.into_iter().map(NodeId::from).collect();
        let i = self.aromatic_systems.push(nodes, data);
        AromaticSystemIdx(i)
    }

    pub fn add_multicenter_bond(
        &mut self,
        atoms: Vec<AtomIdx>,
        data: MulticenterBondAst,
    ) -> MulticenterBondIdx {
        let nodes: Vec<NodeId> = atoms.into_iter().map(NodeId::from).collect();
        let i = self.multicenter_bonds.push(nodes, data);
        MulticenterBondIdx(i)
    }

    pub fn add_noncovalent_bond(
        &mut self,
        ends: [AtomIdx; 2],
        bond: NoncovalentBondAst,
    ) -> NoncovalentBondIdx {
        let i = self
            .noncovalent_bonds
            .push([NodeId::from(ends[0]), NodeId::from(ends[1])], bond);
        NoncovalentBondIdx(i)
    }

    pub fn push_atom_constraint(&mut self, idx: AtomIdx, c: AtomConstraint) {
        self.constraints.push_atom(idx, c);
    }

    pub fn push_bond_constraint(&mut self, idx: BondIdx, c: BondConstraint) {
        self.constraints.push_bond(idx, c);
    }

    pub fn push_dative_bond_constraint(&mut self, idx: DativeBondIdx, c: DativeBondConstraint) {
        self.constraints.push_dative_bond(idx, c);
    }

    pub fn push_aromatic_system_constraint(
        &mut self,
        idx: AromaticSystemIdx,
        c: AromaticSystemConstraint,
    ) {
        self.constraints.push_aromatic_system(idx, c);
    }

    pub fn push_multicenter_bond_constraint(
        &mut self,
        idx: MulticenterBondIdx,
        c: MulticenterBondConstraint,
    ) {
        self.constraints.push_multicenter_bond(idx, c);
    }

    pub fn push_noncovalent_bond_constraint(
        &mut self,
        idx: NoncovalentBondIdx,
        c: NoncovalentBondConstraint,
    ) {
        self.constraints.push_noncovalent_bond(idx, c);
    }

    pub fn push_molecule_constraint(&mut self, c: Constraint) {
        self.constraints.push_molecule(c);
    }

    // -- Attribute mutation ---------------------------------------------------

    pub fn atom_mut(&mut self, idx: AtomIdx) -> &mut AtomAst {
        &mut Arc::make_mut(&mut self.atoms)[idx.index()]
    }

    pub fn bond_mut(&mut self, idx: BondIdx) -> &mut BondAst {
        &mut Arc::make_mut(&mut self.bonds)[idx.index()]
    }

    pub fn dative_bond_mut(&mut self, idx: DativeBondIdx) -> &mut DativeBondAst {
        self.dative_bonds.materialize();
        let FixedSetStorage::Mutable(vec) = &mut self.dative_bonds else {
            unreachable!()
        };
        &mut vec[idx.index()].1
    }

    pub fn aromatic_system_mut(&mut self, idx: AromaticSystemIdx) -> &mut AromaticSystemAst {
        self.aromatic_systems.materialize();
        let VarSetStorage::Mutable(vec) = &mut self.aromatic_systems else {
            unreachable!()
        };
        &mut vec[idx.index()].1
    }

    pub fn multicenter_bond_mut(&mut self, idx: MulticenterBondIdx) -> &mut MulticenterBondAst {
        self.multicenter_bonds.materialize();
        let VarSetStorage::Mutable(vec) = &mut self.multicenter_bonds else {
            unreachable!()
        };
        &mut vec[idx.index()].1
    }

    pub fn noncovalent_bond_mut(&mut self, idx: NoncovalentBondIdx) -> &mut NoncovalentBondAst {
        self.noncovalent_bonds.materialize();
        let FixedSetStorage::Mutable(vec) = &mut self.noncovalent_bonds else {
            unreachable!()
        };
        &mut vec[idx.index()].1
    }

    pub fn constraints_mut(&mut self) -> &mut Constraints {
        &mut self.constraints
    }

    // -- Relation removal -----------------------------------------------------

    pub fn remove_dative_bonds(&mut self, indices: &[DativeBondIdx]) {
        let raw: Vec<u32> = indices.iter().map(|i| i.0).collect();
        self.dative_bonds.remove_indices(&raw);
        let idx_remap = IdxRemapping::new(
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

    pub fn remove_aromatic_systems(&mut self, indices: &[AromaticSystemIdx]) {
        let raw: Vec<u32> = indices.iter().map(|i| i.0).collect();
        self.aromatic_systems.remove_indices(&raw);
        let idx_remap = IdxRemapping::new(
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

    pub fn remove_multicenter_bonds(&mut self, indices: &[MulticenterBondIdx]) {
        let raw: Vec<u32> = indices.iter().map(|i| i.0).collect();
        self.multicenter_bonds.remove_indices(&raw);
        let idx_remap = IdxRemapping::new(
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

    pub fn remove_noncovalent_bonds(&mut self, indices: &[NoncovalentBondIdx]) {
        let raw: Vec<u32> = indices.iter().map(|i| i.0).collect();
        self.noncovalent_bonds.remove_indices(&raw);
        let idx_remap = IdxRemapping::new(
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

    pub fn remove(&mut self, atoms: &[AtomIdx], bonds: &[BondIdx]) -> IdxRemapping {
        let nodes: Vec<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        let edges: Vec<EdgeId> = bonds.iter().map(|&b| EdgeId::from(b)).collect();
        let remap = self.graph.remove(&nodes, &edges);

        let new_atoms = remap.apply_to_node_vec(&self.atoms);
        let new_bonds = remap.apply_to_edge_vec(&self.bonds);
        self.atoms = Arc::new(new_atoms);
        self.bonds = Arc::new(new_bonds);

        let removed_dative = fixed_relation_removed(&self.dative_bonds, &remap);
        let removed_aromatic = var_relation_removed(&self.aromatic_systems, &remap);
        let removed_multicenter = var_relation_removed(&self.multicenter_bonds, &remap);
        let removed_noncovalent = fixed_relation_removed(&self.noncovalent_bonds, &remap);

        let dative = mem::replace(
            &mut self.dative_bonds,
            FixedSetStorage::Shared(Arc::new(FixedRelationSet::default())),
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

        let idx_remap = IdxRemapping::new(
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
