//! Structural editing for `MoleculeAst`. The AST itself only allows attribute mutation;
//! structural change (add atoms/bonds/relations, remove anything) goes through `MoleculeBuilder`.
//!
//! Storage is lazy: each Arc-shared field stays shared until first write, at which point
//! only that field decomposes to a mutable form. `build` re-wraps everything in `Arc`,
//! reusing untouched shared data.

use std::mem;
use std::sync::Arc;

use umol_graph_core::relation::RelationId;
use umol_graph_core::{EdgeId, FixedRelationSet, Graph, NodeId, Remapping, VarRelationSet};

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::constraint::{MoleculeConstraint, MoleculeConstraints};
use super::molecule::MoleculeAst;
use super::multicenter::MulticenterBondAst;
use super::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};

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
}

impl<R: Clone, const N: usize> FixedSetStorage<R, N> {
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
}

pub struct MoleculeBuilder {
    graph: Graph,
    atoms: Arc<Vec<AtomAst>>,
    bonds: Arc<Vec<BondAst>>,
    dative_bonds: FixedSetStorage<BondAst, 2>,
    noncovalent_bonds: FixedSetStorage<BondAst, 2>,
    aromatic_systems: VarSetStorage<AromaticSystemAst>,
    multicenter_bonds: VarSetStorage<MulticenterBondAst>,
    constraints: MoleculeConstraints,
}

impl MoleculeBuilder {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn from_parts(
        graph: Graph,
        atoms: Arc<Vec<AtomAst>>,
        bonds: Arc<Vec<BondAst>>,
        dative_bonds: Arc<FixedRelationSet<BondAst, 2>>,
        noncovalent_bonds: Arc<FixedRelationSet<BondAst, 2>>,
        aromatic_systems: Arc<VarRelationSet<AromaticSystemAst>>,
        multicenter_bonds: Arc<VarRelationSet<MulticenterBondAst>>,
        constraints: MoleculeConstraints,
    ) -> Self {
        Self {
            graph,
            atoms,
            bonds,
            dative_bonds: FixedSetStorage::Shared(dative_bonds),
            noncovalent_bonds: FixedSetStorage::Shared(noncovalent_bonds),
            aromatic_systems: VarSetStorage::Shared(aromatic_systems),
            multicenter_bonds: VarSetStorage::Shared(multicenter_bonds),
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
        bond: BondAst,
    ) -> DativeBondIdx {
        let i = self
            .dative_bonds
            .push([NodeId::from(donor), NodeId::from(acceptor)], bond);
        DativeBondIdx(i)
    }

    pub fn add_noncovalent_bond(
        &mut self,
        ends: [AtomIdx; 2],
        bond: BondAst,
    ) -> NoncovalentBondIdx {
        let i = self
            .noncovalent_bonds
            .push([NodeId::from(ends[0]), NodeId::from(ends[1])], bond);
        NoncovalentBondIdx(i)
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

    pub fn push_constraint(&mut self, c: MoleculeConstraint) {
        self.constraints.insert(c);
    }

    pub fn remove(&mut self, atoms: &[AtomIdx], bonds: &[BondIdx]) -> Remapping {
        let nodes: Vec<NodeId> = atoms.iter().map(|&a| NodeId::from(a)).collect();
        let edges: Vec<EdgeId> = bonds.iter().map(|&b| EdgeId::from(b)).collect();
        let remap = self.graph.remove(&nodes, &edges);

        let new_atoms = remap.apply_to_node_vec(&self.atoms);
        let new_bonds = remap.apply_to_edge_vec(&self.bonds);
        self.atoms = Arc::new(new_atoms);
        self.bonds = Arc::new(new_bonds);

        let dative = mem::replace(
            &mut self.dative_bonds,
            FixedSetStorage::Shared(Arc::new(FixedRelationSet::default())),
        );
        self.dative_bonds = dative.apply_remapping(&remap);

        let noncovalent = mem::replace(
            &mut self.noncovalent_bonds,
            FixedSetStorage::Shared(Arc::new(FixedRelationSet::default())),
        );
        self.noncovalent_bonds = noncovalent.apply_remapping(&remap);

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

        remap
    }

    pub fn build(self) -> MoleculeAst {
        MoleculeAst::from_arcs(
            self.graph,
            self.atoms,
            self.bonds,
            self.dative_bonds.into_arc(),
            self.noncovalent_bonds.into_arc(),
            self.aromatic_systems.into_arc(),
            self.multicenter_bonds.into_arc(),
            self.constraints,
        )
    }
}
