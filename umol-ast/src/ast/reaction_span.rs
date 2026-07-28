//! Reaction span AST: the superimposed `L ∪_K R` graph encoding a reaction's DPO rule span.
//!
//! Materialized superimposed graph carrying, per atom/bond, both its before and after state plus a
//! membership tag. The DPO span `L ←K─ R` is read off the tags — `K = Unchanged ∪ Modified`,
//! `L = K ∪ Removed`, `R = K ∪ Added` — and `rhs()` / `lhs()` project the two sides back to
//! a `MoleculeAst`. `Modified` (a preserved entity relabeled across the reaction) is the
//! relabeling-DPO reading: the entity persists in `K`, its label resolved per side.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::Hash;

use umol_graph_core::{
    Compaction, Correspondence, EdgeId, FixedRelationSet, FixedVarBirelationSet, Graph, NodeId,
    Ordered, RelationId, Unordered, VarRelationSet,
};

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::constraint::{Constraint, Constraints};
use super::correspondence::MoleculeCorrespondence;
use super::dative::DativeBondAst;
use super::delta::{
    apply_aromatic_change, apply_atom_change, apply_bond_change, apply_dative_change,
    apply_multicenter_change, apply_noncovalent_change, apply_stereo_atom_change,
    apply_stereo_bond_change, remap_delta, AromaticSystemDelta, AtomDelta, BondDelta,
    ConstraintDelta, ConstraintSpan, DativeBondDelta, Delta, Deltas, EntityFold, EntitySpan,
    MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta,
};
use super::error::Contradiction;
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigand;
use super::molecule::{MoleculeAst, MoleculeParts};
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;
use super::reaction::ReactionAst;
use super::remap::{IdCompaction, IdRemapping};
use super::stereo::{StereoAtomAst, StereoBondAst};
use super::traits::{Canonicalize, EntityPatch};

/// The superimposed reaction graph — the reaction's DPO rule span, materialized. The union
/// topology is the `lhs` id space (deleted entities kept as nodes/edges) with created entities
/// appended; `atoms` / `bonds` are indexed parallel to the graph's nodes / edges.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionSpanAst {
    graph: Graph,
    atoms: Vec<EntitySpan<AtomAst>>,
    bonds: Vec<EntitySpan<BondAst>>,
    dative_bonds:
        FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, EntitySpan<DativeBondAst>>,
    aromatic_systems: VarRelationSet<NodeId, Unordered, EntitySpan<AromaticSystemAst>>,
    multicenter_bonds: VarRelationSet<NodeId, Unordered, EntitySpan<MulticenterBondAst>>,
    noncovalent_bonds: FixedRelationSet<NodeId, Unordered, EntitySpan<NoncovalentBondAst>, 2>,
    stereo_atoms:
        FixedVarBirelationSet<NodeId, Ordered, 1, StereoLigand, Ordered, EntitySpan<StereoAtomAst>>,
    stereo_bonds:
        FixedVarBirelationSet<EdgeId, Ordered, 1, StereoLigand, Ordered, EntitySpan<StereoBondAst>>,
    constraints: Vec<ConstraintSpan>,
}

/// R id → union id for one entity family: a mated right id reuses its left partner; a right-exposed
/// id is appended after the left family's ids.
fn union_map<Id: Copy + Ord + Hash + From<usize>>(
    correspondence: &Correspondence<Id>,
    left_count: usize,
) -> HashMap<Id, Id> {
    let mut map: HashMap<Id, Id> = correspondence
        .mates()
        .iter()
        .map(|&(l, r)| (r, l))
        .collect();
    for (offset, &r) in correspondence.right_exposed().iter().enumerate() {
        map.insert(r, Id::from(left_count + offset));
    }
    map
}

/// Recover one family's correspondence from its span column (union order): each `Unchanged`/`Modified`
/// pairs the running left/right ids, `Removed` advances only left, `Added` only right.
fn recover_correspondence<'a, Id, T: 'a>(
    column: impl Iterator<Item = &'a EntitySpan<T>>,
) -> Correspondence<Id>
where
    Id: Copy + Ord + From<usize>,
{
    let mut mates = Vec::new();
    let (mut left, mut right) = (0usize, 0usize);
    for span in column {
        match span {
            EntitySpan::Unchanged(_) | EntitySpan::Modified { .. } => {
                mates.push((Id::from(left), Id::from(right)));
                left += 1;
                right += 1;
            }
            EntitySpan::Removed(_) => left += 1,
            EntitySpan::Added(_) => right += 1,
        }
    }
    Correspondence::new(mates, left, right)
}

impl MoleculeAst {
    /// The deltas transforming `self` (`L`) into `rhs` (`R`) under the per-entity correspondence
    /// `correspondence`: superimpose the two sides into a span, then read off its operational
    /// reaction. Applying the result to `self` reconstructs `rhs` (up to appended-atom renumbering).
    pub fn difference_to(
        &self,
        rhs: &MoleculeAst,
        correspondence: &MoleculeCorrespondence,
    ) -> Deltas {
        ReactionSpanAst::superimpose(self, rhs, correspondence)
            .to_reaction()
            .deltas
    }
}

#[allow(clippy::type_complexity)]
impl ReactionSpanAst {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn from_parts(
        graph: Graph,
        atoms: Vec<EntitySpan<AtomAst>>,
        bonds: Vec<EntitySpan<BondAst>>,
        dative_bonds: FixedVarBirelationSet<
            NodeId,
            Ordered,
            1,
            NodeId,
            Unordered,
            EntitySpan<DativeBondAst>,
        >,
        aromatic_systems: VarRelationSet<NodeId, Unordered, EntitySpan<AromaticSystemAst>>,
        multicenter_bonds: VarRelationSet<NodeId, Unordered, EntitySpan<MulticenterBondAst>>,
        noncovalent_bonds: FixedRelationSet<NodeId, Unordered, EntitySpan<NoncovalentBondAst>, 2>,
        stereo_atoms: FixedVarBirelationSet<
            NodeId,
            Ordered,
            1,
            StereoLigand,
            Ordered,
            EntitySpan<StereoAtomAst>,
        >,
        stereo_bonds: FixedVarBirelationSet<
            EdgeId,
            Ordered,
            1,
            StereoLigand,
            Ordered,
            EntitySpan<StereoBondAst>,
        >,
        constraints: Vec<ConstraintSpan>,
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
        }
    }

    /// Superimpose two molecules over their correspondence into the reaction span. Mated entities
    /// become `Unchanged` / `Modified` carrying both molecules' actual values; entities exposed on
    /// the lhs become `Removed`, those exposed on the rhs `Added`. Lhs-anchored: lhs ids kept,
    /// rhs-exposed entities appended, rhs participants and constraints remapped into that union frame.
    pub fn superimpose(
        lhs: &MoleculeAst,
        rhs: &MoleculeAst,
        correspondence: &MoleculeCorrespondence,
    ) -> ReactionSpanAst {
        let atoms_corr = correspondence.atoms();
        let bonds_corr = correspondence.bonds();
        let lhs_atom_count = lhs.atoms().count();
        let lhs_bond_count = lhs.bonds().count();

        // R id → union id per family.
        let atom_union: HashMap<NodeId, NodeId> = union_map(atoms_corr, lhs_atom_count);
        let bond_union: HashMap<BondId, BondId> = union_map(bonds_corr, lhs_bond_count);

        // Atoms
        let mut atoms: Vec<EntitySpan<AtomAst>> = Vec::new();
        for i in 0..lhs_atom_count {
            let lhs_ast = lhs.atom(AtomId(i as u32)).ast.clone();
            let rhs_ast = atoms_corr
                .right_of(NodeId(i as u32))
                .map(|r| rhs.atom(AtomId::from(r.index())).ast.clone());
            atoms.push(EntitySpan::superimpose(Some(lhs_ast), rhs_ast).unwrap());
        }
        for &r in &atoms_corr.right_exposed() {
            atoms.push(EntitySpan::Added(
                rhs.atom(AtomId::from(r.index())).ast.clone(),
            ));
        }

        // Bonds + union edges
        let mut bonds: Vec<EntitySpan<BondAst>> = Vec::new();
        let mut edges: Vec<[u32; 2]> = Vec::new();
        for i in 0..lhs_bond_count {
            let [a, b] = lhs.raw_graph().edge_endpoints(EdgeId(i as u32));
            edges.push([a.0, b.0]);
            let lhs_ast = lhs.bond(BondId(i as u32)).ast.clone();
            let rhs_ast = bonds_corr
                .right_of(BondId(i as u32))
                .map(|r| rhs.bond(r).ast.clone());
            bonds.push(EntitySpan::superimpose(Some(lhs_ast), rhs_ast).unwrap());
        }
        for &r in &bonds_corr.right_exposed() {
            let [a, b] = rhs.raw_graph().edge_endpoints(EdgeId(r.index() as u32));
            edges.push([atom_union[&a].0, atom_union[&b].0]);
            bonds.push(EntitySpan::Added(rhs.bond(r).ast.clone()));
        }

        // Aromatic systems
        let aromatic_corr = correspondence.aromatic_systems();
        let mut aromatic_entries: Vec<(Vec<NodeId>, EntitySpan<AromaticSystemAst>)> = Vec::new();
        for view in lhs.aromatic_systems().iter() {
            let participants: Vec<NodeId> = view.atom_ids().map(NodeId::from).collect();
            let rhs_ast = aromatic_corr
                .right_of(view.id)
                .map(|r| rhs.aromatic_system(r).ast.clone());
            aromatic_entries.push((
                participants,
                EntitySpan::superimpose(Some(view.ast.clone()), rhs_ast).unwrap(),
            ));
        }
        for &r in &aromatic_corr.right_exposed() {
            let view = rhs.aromatic_system(r);
            let participants: Vec<NodeId> = view
                .atom_ids()
                .map(|a| atom_union[&NodeId::from(a)])
                .collect();
            aromatic_entries.push((participants, EntitySpan::Added(view.ast.clone())));
        }
        let aromatic_systems = VarRelationSet::new(aromatic_entries);

        // Multicenter bonds
        let multicenter_corr = correspondence.multicenter_bonds();
        let mut multicenter_entries: Vec<(Vec<NodeId>, EntitySpan<MulticenterBondAst>)> =
            Vec::new();
        for view in lhs.multicenter_bonds().iter() {
            let participants: Vec<NodeId> = view.atom_ids().map(NodeId::from).collect();
            let rhs_ast = multicenter_corr
                .right_of(view.id)
                .map(|r| rhs.multicenter_bond(r).ast.clone());
            multicenter_entries.push((
                participants,
                EntitySpan::superimpose(Some(view.ast.clone()), rhs_ast).unwrap(),
            ));
        }
        for &r in &multicenter_corr.right_exposed() {
            let view = rhs.multicenter_bond(r);
            let participants: Vec<NodeId> = view
                .atom_ids()
                .map(|a| atom_union[&NodeId::from(a)])
                .collect();
            multicenter_entries.push((participants, EntitySpan::Added(view.ast.clone())));
        }
        let multicenter_bonds = VarRelationSet::new(multicenter_entries);

        // Noncovalent bonds
        let noncovalent_corr = correspondence.noncovalent_bonds();
        let mut noncovalent_entries: Vec<([NodeId; 2], EntitySpan<NoncovalentBondAst>)> =
            Vec::new();
        for view in lhs.noncovalent_bonds().iter() {
            let [a, b] = view.atom_ids();
            let rhs_ast = noncovalent_corr
                .right_of(view.id)
                .map(|r| rhs.noncovalent_bond(r).ast.clone());
            noncovalent_entries.push((
                [NodeId::from(a), NodeId::from(b)],
                EntitySpan::superimpose(Some(view.ast.clone()), rhs_ast).unwrap(),
            ));
        }
        for &r in &noncovalent_corr.right_exposed() {
            let view = rhs.noncovalent_bond(r);
            let [a, b] = view.atom_ids();
            let participants = [atom_union[&NodeId::from(a)], atom_union[&NodeId::from(b)]];
            noncovalent_entries.push((participants, EntitySpan::Added(view.ast.clone())));
        }
        let noncovalent_bonds = FixedRelationSet::new(noncovalent_entries);

        // Dative bonds
        let dative_corr = correspondence.dative_bonds();
        let mut dative_entries: Vec<([NodeId; 1], Vec<NodeId>, EntitySpan<DativeBondAst>)> =
            Vec::new();
        for view in lhs.dative_bonds().iter() {
            let acceptor = [NodeId::from(view.acceptor_id())];
            let donors: Vec<NodeId> = view.donor_ids().map(NodeId::from).collect();
            let rhs_ast = dative_corr
                .right_of(view.id)
                .map(|r| rhs.dative_bond(r).ast.clone());
            dative_entries.push((
                acceptor,
                donors,
                EntitySpan::superimpose(Some(view.ast.clone()), rhs_ast).unwrap(),
            ));
        }
        for &r in &dative_corr.right_exposed() {
            let view = rhs.dative_bond(r);
            let acceptor = [atom_union[&NodeId::from(view.acceptor_id())]];
            let donors: Vec<NodeId> = view
                .donor_ids()
                .map(|a| atom_union[&NodeId::from(a)])
                .collect();
            dative_entries.push((acceptor, donors, EntitySpan::Added(view.ast.clone())));
        }
        let dative_bonds = FixedVarBirelationSet::new(dative_entries);

        // Stereo atoms
        let stereo_atom_corr = correspondence.stereo_atoms();
        let mut stereo_atom_entries: Vec<(
            [NodeId; 1],
            Vec<StereoLigand>,
            EntitySpan<StereoAtomAst>,
        )> = Vec::new();
        for view in lhs.stereo_atoms().iter() {
            let rhs_ast = stereo_atom_corr
                .right_of(view.id)
                .map(|r| rhs.stereo_atom(r).ast.clone());
            stereo_atom_entries.push((
                [NodeId::from(view.site_id())],
                view.ligand_frame(),
                EntitySpan::superimpose(Some(view.ast.clone()), rhs_ast).unwrap(),
            ));
        }
        for &r in &stereo_atom_corr.right_exposed() {
            let view = rhs.stereo_atom(r);
            let site = [atom_union[&NodeId::from(view.site_id())]];
            let ligands: Vec<StereoLigand> = view
                .ligand_frame()
                .iter()
                .map(|l| {
                    StereoLigand::new(
                        AtomId::from(atom_union[&NodeId::from(l.atom_id)].index()),
                        l.kind,
                    )
                })
                .collect();
            stereo_atom_entries.push((site, ligands, EntitySpan::Added(view.ast.clone())));
        }
        let stereo_atoms = FixedVarBirelationSet::new(stereo_atom_entries);

        // Stereo bonds
        let stereo_bond_corr = correspondence.stereo_bonds();
        let mut stereo_bond_entries: Vec<(
            [EdgeId; 1],
            Vec<StereoLigand>,
            EntitySpan<StereoBondAst>,
        )> = Vec::new();
        for view in lhs.stereo_bonds().iter() {
            let rhs_ast = stereo_bond_corr
                .right_of(view.id)
                .map(|r| rhs.stereo_bond(r).ast.clone());
            stereo_bond_entries.push((
                [EdgeId::from(view.site_id())],
                view.ligand_frame(),
                EntitySpan::superimpose(Some(view.ast.clone()), rhs_ast).unwrap(),
            ));
        }
        for &r in &stereo_bond_corr.right_exposed() {
            let view = rhs.stereo_bond(r);
            let site = [EdgeId::from(bond_union[&view.site_id()])];
            let ligands: Vec<StereoLigand> = view
                .ligand_frame()
                .iter()
                .map(|l| {
                    StereoLigand::new(
                        AtomId::from(atom_union[&NodeId::from(l.atom_id)].index()),
                        l.kind,
                    )
                })
                .collect();
            stereo_bond_entries.push((site, ligands, EntitySpan::Added(view.ast.clone())));
        }
        let stereo_bonds = FixedVarBirelationSet::new(stereo_bond_entries);

        // Constraints: R's remapped into the union frame, then set-diffed against L's.
        let remapping = IdRemapping::new(
            atom_union
                .iter()
                .map(|(&r, &u)| (AtomId::from(r.index()), AtomId::from(u.index())))
                .collect(),
            bond_union,
            union_map(dative_corr, lhs.dative_bonds().count()),
            union_map(aromatic_corr, lhs.aromatic_systems().count()),
            union_map(multicenter_corr, lhs.multicenter_bonds().count()),
            union_map(noncovalent_corr, lhs.noncovalent_bonds().count()),
            union_map(stereo_atom_corr, lhs.stereo_atoms().count()),
            union_map(stereo_bond_corr, lhs.stereo_bonds().count()),
        );
        let rhs_constraints: Vec<Constraint> = rhs
            .constraints()
            .iter()
            .cloned()
            .map(|c| c.remap(&remapping))
            .collect();
        let mut constraints: Vec<ConstraintSpan> = Vec::new();
        for c in lhs.constraints().iter() {
            if rhs_constraints.contains(c) {
                constraints.push(ConstraintSpan::Unchanged(c.clone()));
            } else {
                constraints.push(ConstraintSpan::Removed(c.clone()));
            }
        }
        for c in &rhs_constraints {
            if !lhs.constraints().iter().any(|l| l == c) {
                constraints.push(ConstraintSpan::Added(c.clone()));
            }
        }

        let graph = Graph::new(atoms.len(), &edges);
        ReactionSpanAst::from_parts(
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
        )
    }

    /// Recover the per-family correspondence between the two sides, forgetting the values — the
    /// inverse of `superimpose`: `superimpose(self.lhs(), self.rhs(), &self.correspondence())`
    /// reproduces `self`.
    pub fn correspondence(&self) -> MoleculeCorrespondence {
        MoleculeCorrespondence::new(
            recover_correspondence(self.atoms.iter()),
            recover_correspondence(self.bonds.iter()),
            recover_correspondence(
                self.dative_bonds
                    .relation_ids()
                    .map(|r| self.dative_bonds.data(r)),
            ),
            recover_correspondence(
                self.aromatic_systems
                    .relation_ids()
                    .map(|r| self.aromatic_systems.data(r)),
            ),
            recover_correspondence(
                self.multicenter_bonds
                    .relation_ids()
                    .map(|r| self.multicenter_bonds.data(r)),
            ),
            recover_correspondence(
                self.noncovalent_bonds
                    .relation_ids()
                    .map(|r| self.noncovalent_bonds.data(r)),
            ),
            recover_correspondence(
                self.stereo_atoms
                    .relation_ids()
                    .map(|r| self.stereo_atoms.data(r)),
            ),
            recover_correspondence(
                self.stereo_bonds
                    .relation_ids()
                    .map(|r| self.stereo_bonds.data(r)),
            ),
        )
    }

    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    pub fn atoms(&self) -> &[EntitySpan<AtomAst>] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[EntitySpan<BondAst>] {
        &self.bonds
    }

    pub(crate) fn dative_bonds(
        &self,
    ) -> &FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, EntitySpan<DativeBondAst>>
    {
        &self.dative_bonds
    }

    pub(crate) fn aromatic_systems(
        &self,
    ) -> &VarRelationSet<NodeId, Unordered, EntitySpan<AromaticSystemAst>> {
        &self.aromatic_systems
    }

    pub(crate) fn multicenter_bonds(
        &self,
    ) -> &VarRelationSet<NodeId, Unordered, EntitySpan<MulticenterBondAst>> {
        &self.multicenter_bonds
    }

    pub(crate) fn noncovalent_bonds(
        &self,
    ) -> &FixedRelationSet<NodeId, Unordered, EntitySpan<NoncovalentBondAst>, 2> {
        &self.noncovalent_bonds
    }

    pub(crate) fn stereo_atoms(
        &self,
    ) -> &FixedVarBirelationSet<NodeId, Ordered, 1, StereoLigand, Ordered, EntitySpan<StereoAtomAst>>
    {
        &self.stereo_atoms
    }

    pub(crate) fn stereo_bonds(
        &self,
    ) -> &FixedVarBirelationSet<EdgeId, Ordered, 1, StereoLigand, Ordered, EntitySpan<StereoBondAst>>
    {
        &self.stereo_bonds
    }

    pub fn constraints(&self) -> &[ConstraintSpan] {
        &self.constraints
    }

    /// The left-hand molecule: every entity present on the left, in a compacted id space (created
    /// entities dropped).
    pub fn lhs(&self) -> MoleculeAst {
        self.project(Side::Left)
    }

    /// The right-hand molecule: every entity present on the right, in a compacted id space (deleted
    /// entities dropped).
    pub fn rhs(&self) -> MoleculeAst {
        self.project(Side::Right)
    }

    /// Recover the operational `ReactionAst` from the span — the inverse of
    /// `ReactionAst::to_reaction_span`, up to delta normal form. The reaction's `lhs` is `self.lhs()`
    /// (which preserves the original lhs id space); each entity's `EntitySpan` yields its delta, a
    /// `Modified` one via an AST-diff of its left/right values.
    pub fn to_reaction(&self) -> ReactionAst {
        let mut deltas = AtomDelta::deltas_from_states(&self.atoms, |_| ());
        deltas.extend(BondDelta::deltas_from_states(&self.bonds, |edge| {
            let [a, b] = self.graph.edge_endpoints(EdgeId(edge as u32));
            [AtomId::from(a), AtomId::from(b)]
        }));
        let dative_states: Vec<EntitySpan<DativeBondAst>> = (0..self.dative_bonds.count())
            .map(|i| self.dative_bonds.data(RelationId(i as u32)).clone())
            .collect();
        deltas.extend(DativeBondDelta::deltas_from_states(
            &dative_states,
            |index| {
                let rid = RelationId(index as u32);
                (
                    self.dative_bonds
                        .participants_2(rid)
                        .iter()
                        .map(|&n| AtomId::from(n))
                        .collect(),
                    AtomId::from(self.dative_bonds.participants_1(rid)[0]),
                )
            },
        ));
        let aromatic_states: Vec<EntitySpan<AromaticSystemAst>> =
            (0..self.aromatic_systems.count())
                .map(|i| self.aromatic_systems.data(RelationId(i as u32)).clone())
                .collect();
        deltas.extend(AromaticSystemDelta::deltas_from_states(
            &aromatic_states,
            |index| {
                self.aromatic_systems
                    .participants(RelationId(index as u32))
                    .iter()
                    .map(|&n| AtomId::from(n))
                    .collect()
            },
        ));
        let multicenter_states: Vec<EntitySpan<MulticenterBondAst>> =
            (0..self.multicenter_bonds.count())
                .map(|i| self.multicenter_bonds.data(RelationId(i as u32)).clone())
                .collect();
        deltas.extend(MulticenterBondDelta::deltas_from_states(
            &multicenter_states,
            |index| {
                self.multicenter_bonds
                    .participants(RelationId(index as u32))
                    .iter()
                    .map(|&n| AtomId::from(n))
                    .collect()
            },
        ));
        let noncovalent_states: Vec<EntitySpan<NoncovalentBondAst>> =
            (0..self.noncovalent_bonds.count())
                .map(|i| self.noncovalent_bonds.data(RelationId(i as u32)).clone())
                .collect();
        deltas.extend(NoncovalentBondDelta::deltas_from_states(
            &noncovalent_states,
            |index| {
                let [a, b] = *self
                    .noncovalent_bonds
                    .participants(RelationId(index as u32));
                [AtomId::from(a), AtomId::from(b)]
            },
        ));
        // Stereo overlays have no `EntityFold`, so recover their deltas here: `Removed`/`Added` carry
        // the relation's site + ligand frame; `Modified` is the field/constraint diff. Site/ligand
        // ids are the union frame (lhs ids for preserved/removed entities).
        for i in 0..self.stereo_atoms.count() {
            let rid = RelationId(i as u32);
            let id = StereoAtomId::from(rid);
            let site = AtomId::from(self.stereo_atoms.participants_1(rid)[0]);
            let ligands = self.stereo_atoms.participants_2(rid).to_vec();
            match self.stereo_atoms.data(rid) {
                EntitySpan::Unchanged(_) => {}
                EntitySpan::Added(ast) => deltas.push(Delta::StereoAtom(StereoAtomDelta::Add {
                    id,
                    site,
                    ligands,
                    ast: ast.clone(),
                })),
                EntitySpan::Removed(ast) => {
                    deltas.push(Delta::StereoAtom(StereoAtomDelta::Remove {
                        id,
                        site,
                        ligands,
                        ast: ast.clone(),
                    }))
                }
                EntitySpan::Modified {
                    lhs: left,
                    rhs: right,
                } => deltas.extend(
                    StereoAtomDelta::diff(id, left, right)
                        .into_iter()
                        .map(Delta::StereoAtom),
                ),
            }
        }
        for i in 0..self.stereo_bonds.count() {
            let rid = RelationId(i as u32);
            let id = StereoBondId::from(rid);
            let site = BondId::from(self.stereo_bonds.participants_1(rid)[0]);
            let ligands = self.stereo_bonds.participants_2(rid).to_vec();
            match self.stereo_bonds.data(rid) {
                EntitySpan::Unchanged(_) => {}
                EntitySpan::Added(ast) => deltas.push(Delta::StereoBond(StereoBondDelta::Add {
                    id,
                    site,
                    ligands,
                    ast: ast.clone(),
                })),
                EntitySpan::Removed(ast) => {
                    deltas.push(Delta::StereoBond(StereoBondDelta::Remove {
                        id,
                        site,
                        ligands,
                        ast: ast.clone(),
                    }))
                }
                EntitySpan::Modified {
                    lhs: left,
                    rhs: right,
                } => deltas.extend(
                    StereoBondDelta::diff(id, left, right)
                        .into_iter()
                        .map(Delta::StereoBond),
                ),
            }
        }
        for span in &self.constraints {
            match span {
                ConstraintSpan::Added(c) => {
                    deltas.push(Delta::Constraint(ConstraintDelta::Add(c.clone())))
                }
                ConstraintSpan::Removed(c) => {
                    deltas.push(Delta::Constraint(ConstraintDelta::Remove(c.clone())))
                }
                ConstraintSpan::Unchanged(_) => {}
            }
        }
        ReactionAst::new(self.lhs(), Deltas::from_iter(deltas))
    }

    /// Project one `Side` to a `MoleculeAst`: every entity present on that side, in a compacted id
    /// space (entities absent on the side dropped, the survivors renumbered). Overlays are carried
    /// through (an overlay is dropped if its side value is absent or any participant is dropped),
    /// and the side's constraints are compacted through the same compaction (atom/bond and overlay
    /// refs to a dropped entity are dropped).
    fn project(&self, side: Side) -> MoleculeAst {
        let mut compacted: Vec<Option<AtomId>> = vec![None; self.atoms.len()];
        let mut atoms: Vec<AtomAst> = Vec::new();
        let mut removed_nodes: Vec<u32> = Vec::new();
        for (node, change) in self.atoms.iter().enumerate() {
            match entity_side(change, side) {
                Some(ast) => {
                    compacted[node] = Some(AtomId(atoms.len() as u32));
                    atoms.push(ast.clone());
                }
                None => removed_nodes.push(node as u32),
            }
        }
        let mut compacted_bonds: Vec<Option<BondId>> = vec![None; self.bonds.len()];
        let mut bonds: Vec<(AtomId, AtomId, BondAst)> = Vec::new();
        let mut removed_edges: Vec<u32> = Vec::new();
        for (edge, change) in self.bonds.iter().enumerate() {
            let [a, b] = self.graph.edge_endpoints(EdgeId(edge as u32));
            match (
                entity_side(change, side),
                compacted[a.index()],
                compacted[b.index()],
            ) {
                (Some(ast), Some(a), Some(b)) => {
                    compacted_bonds[edge] = Some(BondId(bonds.len() as u32));
                    bonds.push((a, b, ast.clone()));
                }
                _ => removed_edges.push(edge as u32),
            }
        }

        let atom = |n: NodeId| compacted[n.index()];
        let mut dative: Vec<(Vec<AtomId>, AtomId, DativeBondAst)> = Vec::new();
        let mut removed_dative: Vec<RelationId> = Vec::new();
        for i in 0..self.dative_bonds.count() {
            let rid = RelationId(i as u32);
            let acceptor = atom(self.dative_bonds.participants_1(rid)[0]);
            let donors: Option<Vec<AtomId>> = self
                .dative_bonds
                .participants_2(rid)
                .iter()
                .map(|&n| atom(n))
                .collect();
            match (
                entity_side(self.dative_bonds.data(rid), side),
                acceptor,
                donors,
            ) {
                (Some(ast), Some(acceptor), Some(donors)) => {
                    dative.push((donors, acceptor, ast.clone()))
                }
                _ => removed_dative.push(rid),
            }
        }

        let mut aromatic: Vec<(Vec<AtomId>, AromaticSystemAst)> = Vec::new();
        let mut removed_aromatic: Vec<RelationId> = Vec::new();
        for i in 0..self.aromatic_systems.count() {
            let rid = RelationId(i as u32);
            let members: Option<Vec<AtomId>> = self
                .aromatic_systems
                .participants(rid)
                .iter()
                .map(|&n| atom(n))
                .collect();
            match (entity_side(self.aromatic_systems.data(rid), side), members) {
                (Some(ast), Some(members)) => aromatic.push((members, ast.clone())),
                _ => removed_aromatic.push(rid),
            }
        }

        let mut multicenter: Vec<(Vec<AtomId>, MulticenterBondAst)> = Vec::new();
        let mut removed_multicenter: Vec<RelationId> = Vec::new();
        for i in 0..self.multicenter_bonds.count() {
            let rid = RelationId(i as u32);
            let members: Option<Vec<AtomId>> = self
                .multicenter_bonds
                .participants(rid)
                .iter()
                .map(|&n| atom(n))
                .collect();
            match (entity_side(self.multicenter_bonds.data(rid), side), members) {
                (Some(ast), Some(members)) => multicenter.push((members, ast.clone())),
                _ => removed_multicenter.push(rid),
            }
        }

        let mut noncovalent: Vec<(AtomId, AtomId, NoncovalentBondAst)> = Vec::new();
        let mut removed_noncovalent: Vec<RelationId> = Vec::new();
        for i in 0..self.noncovalent_bonds.count() {
            let rid = RelationId(i as u32);
            let [a, b] = *self.noncovalent_bonds.participants(rid);
            match (
                entity_side(self.noncovalent_bonds.data(rid), side),
                atom(a),
                atom(b),
            ) {
                (Some(ast), Some(a), Some(b)) => noncovalent.push((a, b, ast.clone())),
                _ => removed_noncovalent.push(rid),
            }
        }

        let ligands = |ls: &[StereoLigand]| -> Option<Vec<StereoLigand>> {
            ls.iter()
                .map(|l| atom(NodeId::from(l.atom_id)).map(|a| StereoLigand::new(a, l.kind)))
                .collect()
        };
        let mut stereo_atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)> = Vec::new();
        let mut removed_stereo_atoms: Vec<RelationId> = Vec::new();
        for i in 0..self.stereo_atoms.count() {
            let rid = RelationId(i as u32);
            let site = atom(self.stereo_atoms.participants_1(rid)[0]);
            match (
                entity_side(self.stereo_atoms.data(rid), side),
                site,
                ligands(self.stereo_atoms.participants_2(rid)),
            ) {
                (Some(ast), Some(site), Some(ligands)) => {
                    stereo_atoms.push((site, ligands, ast.clone()))
                }
                _ => removed_stereo_atoms.push(rid),
            }
        }

        let mut stereo_bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondAst)> = Vec::new();
        let mut removed_stereo_bonds: Vec<RelationId> = Vec::new();
        for i in 0..self.stereo_bonds.count() {
            let rid = RelationId(i as u32);
            let site = compacted_bonds[self.stereo_bonds.participants_1(rid)[0].index()];
            match (
                entity_side(self.stereo_bonds.data(rid), side),
                site,
                ligands(self.stereo_bonds.participants_2(rid)),
            ) {
                (Some(ast), Some(site), Some(ligands)) => {
                    stereo_bonds.push((site, ligands, ast.clone()))
                }
                _ => removed_stereo_bonds.push(rid),
            }
        }

        let compaction = IdCompaction::new(
            Compaction::new(removed_nodes, removed_edges),
            removed_dative,
            removed_aromatic,
            removed_multicenter,
            removed_noncovalent,
            removed_stereo_atoms,
            removed_stereo_bonds,
        );
        let mut constraints = Constraints::new();
        for span in &self.constraints {
            let value = match side {
                Side::Left => span.lhs(),
                Side::Right => span.rhs(),
            };
            if let Some(c) = value {
                constraints.push(c.clone());
            }
        }
        constraints.compact(&compaction);

        MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            stereo_atoms,
            stereo_bonds,
            constraints,
        })
    }
}

/// Which side of the span a projection reads.
#[derive(Clone, Copy)]
enum Side {
    Left,
    Right,
}

/// The side value of an entity span (`None` if the entity is absent on that side).
fn entity_side<T>(span: &EntitySpan<T>, side: Side) -> Option<&T> {
    match side {
        Side::Left => span.lhs(),
        Side::Right => span.rhs(),
    }
}

impl ReactionAst {
    /// Materialize the superimposed reaction span. Canonicalizes the deltas, then annotates
    /// each `lhs` entity (in its own id space) with its before/after state — `Removed` /
    /// `Added` / `Modified` / `Unchanged` — appending created entities. A `Modified` entity's
    /// right value is its left value with the entity's field/constraint changes applied.
    /// `Err(Contradiction)` if the deltas are inconsistent (or inconsistent with `lhs`).
    pub fn to_reaction_span(&self) -> Result<ReactionSpanAst, Contradiction> {
        let deltas = self.deltas.clone().canonicalize()?;
        let lhs = &self.lhs;
        let atom_count = lhs.atoms().count();
        let bond_count = lhs.bonds().count();

        let mut removed_atoms: HashMap<AtomId, AtomAst> = HashMap::new();
        let mut added_atoms: BTreeMap<AtomId, AtomAst> = BTreeMap::new();
        let mut atom_changes: HashMap<AtomId, Vec<AtomDelta>> = HashMap::new();
        let mut removed_bonds: HashMap<BondId, BondAst> = HashMap::new();
        let mut added_bonds: BTreeMap<BondId, ([AtomId; 2], BondAst)> = BTreeMap::new();
        let mut bond_changes: HashMap<BondId, Vec<BondDelta>> = HashMap::new();
        let mut removed_aromatic: HashMap<AromaticSystemId, AromaticSystemAst> = HashMap::new();
        let mut added_aromatic: BTreeMap<AromaticSystemId, (Vec<AtomId>, AromaticSystemAst)> =
            BTreeMap::new();
        let mut aromatic_changes: HashMap<AromaticSystemId, Vec<AromaticSystemDelta>> =
            HashMap::new();
        let mut removed_dative: HashMap<DativeBondId, DativeBondAst> = HashMap::new();
        let mut added_dative: BTreeMap<DativeBondId, (Vec<AtomId>, AtomId, DativeBondAst)> =
            BTreeMap::new();
        let mut dative_changes: HashMap<DativeBondId, Vec<DativeBondDelta>> = HashMap::new();
        let mut removed_multicenter: HashMap<MulticenterBondId, MulticenterBondAst> =
            HashMap::new();
        let mut added_multicenter: BTreeMap<MulticenterBondId, (Vec<AtomId>, MulticenterBondAst)> =
            BTreeMap::new();
        let mut multicenter_changes: HashMap<MulticenterBondId, Vec<MulticenterBondDelta>> =
            HashMap::new();
        let mut removed_noncovalent: HashMap<NoncovalentBondId, NoncovalentBondAst> =
            HashMap::new();
        let mut added_noncovalent: BTreeMap<NoncovalentBondId, ([AtomId; 2], NoncovalentBondAst)> =
            BTreeMap::new();
        let mut noncovalent_changes: HashMap<NoncovalentBondId, Vec<NoncovalentBondDelta>> =
            HashMap::new();
        let mut removed_stereo_atom: HashMap<StereoAtomId, StereoAtomAst> = HashMap::new();
        let mut added_stereo_atom: BTreeMap<
            StereoAtomId,
            (AtomId, Vec<StereoLigand>, StereoAtomAst),
        > = BTreeMap::new();
        let mut stereo_atom_changes: HashMap<StereoAtomId, Vec<StereoAtomDelta>> = HashMap::new();
        let mut removed_stereo_bond: HashMap<StereoBondId, StereoBondAst> = HashMap::new();
        let mut added_stereo_bond: BTreeMap<
            StereoBondId,
            (BondId, Vec<StereoLigand>, StereoBondAst),
        > = BTreeMap::new();
        let mut stereo_bond_changes: HashMap<StereoBondId, Vec<StereoBondDelta>> = HashMap::new();
        let mut added_constraints: Vec<Constraint> = Vec::new();
        let mut removed_constraints: Vec<Constraint> = Vec::new();

        for delta in deltas.iter() {
            match delta {
                Delta::Atom(atom) => match atom {
                    AtomDelta::Remove { id, ast } => {
                        removed_atoms.insert(*id, ast.clone());
                    }
                    AtomDelta::Add { id, ast } => {
                        added_atoms.insert(*id, ast.clone());
                    }
                    AtomDelta::ModifyField { id, .. } | AtomDelta::ModifyConstraint { id, .. } => {
                        atom_changes.entry(*id).or_default().push(atom.clone());
                    }
                },
                Delta::Bond(bond) => match bond {
                    BondDelta::Remove { id, ast, .. } => {
                        removed_bonds.insert(*id, ast.clone());
                    }
                    BondDelta::Add { id, atoms, ast } => {
                        added_bonds.insert(*id, (*atoms, ast.clone()));
                    }
                    BondDelta::ModifyField { id, .. } | BondDelta::ModifyConstraint { id, .. } => {
                        bond_changes.entry(*id).or_default().push(bond.clone());
                    }
                },
                Delta::AromaticSystem(aromatic) => match aromatic {
                    AromaticSystemDelta::Remove { id, ast, .. } => {
                        removed_aromatic.insert(*id, ast.clone());
                    }
                    AromaticSystemDelta::Add { id, atoms, ast } => {
                        added_aromatic.insert(*id, (atoms.clone(), ast.clone()));
                    }
                    AromaticSystemDelta::ModifyField { id, .. }
                    | AromaticSystemDelta::ModifyConstraint { id, .. } => {
                        aromatic_changes
                            .entry(*id)
                            .or_default()
                            .push(aromatic.clone());
                    }
                },
                Delta::DativeBond(dative) => match dative {
                    DativeBondDelta::Remove { id, ast, .. } => {
                        removed_dative.insert(*id, ast.clone());
                    }
                    DativeBondDelta::Add {
                        id,
                        donors,
                        acceptor,
                        ast,
                    } => {
                        added_dative.insert(*id, (donors.clone(), *acceptor, ast.clone()));
                    }
                    DativeBondDelta::ModifyField { id, .. }
                    | DativeBondDelta::ModifyConstraint { id, .. } => {
                        dative_changes.entry(*id).or_default().push(dative.clone());
                    }
                },
                Delta::MulticenterBond(multicenter) => match multicenter {
                    MulticenterBondDelta::Remove { id, ast, .. } => {
                        removed_multicenter.insert(*id, ast.clone());
                    }
                    MulticenterBondDelta::Add { id, atoms, ast } => {
                        added_multicenter.insert(*id, (atoms.clone(), ast.clone()));
                    }
                    MulticenterBondDelta::ModifyField { id, .. }
                    | MulticenterBondDelta::ModifyConstraint { id, .. } => {
                        multicenter_changes
                            .entry(*id)
                            .or_default()
                            .push(multicenter.clone());
                    }
                },
                Delta::NoncovalentBond(noncovalent) => match noncovalent {
                    NoncovalentBondDelta::Remove { id, ast, .. } => {
                        removed_noncovalent.insert(*id, ast.clone());
                    }
                    NoncovalentBondDelta::Add { id, atoms, ast } => {
                        added_noncovalent.insert(*id, (*atoms, ast.clone()));
                    }
                    NoncovalentBondDelta::ModifyField { id, .. }
                    | NoncovalentBondDelta::ModifyConstraint { id, .. } => {
                        noncovalent_changes
                            .entry(*id)
                            .or_default()
                            .push(noncovalent.clone());
                    }
                },
                Delta::StereoAtom(stereo) => match stereo {
                    StereoAtomDelta::Remove { id, ast, .. } => {
                        removed_stereo_atom.insert(*id, ast.clone());
                    }
                    StereoAtomDelta::Add {
                        id,
                        site,
                        ligands,
                        ast,
                    } => {
                        added_stereo_atom.insert(*id, (*site, ligands.clone(), ast.clone()));
                    }
                    StereoAtomDelta::ModifyField { id, .. }
                    | StereoAtomDelta::ModifyConstraint { id, .. }
                    | StereoAtomDelta::Apply { id, .. }
                    | StereoAtomDelta::Swap { id, .. }
                    | StereoAtomDelta::Mirror { id, .. } => {
                        stereo_atom_changes
                            .entry(*id)
                            .or_default()
                            .push(stereo.clone());
                    }
                },
                Delta::StereoBond(stereo) => match stereo {
                    StereoBondDelta::Remove { id, ast, .. } => {
                        removed_stereo_bond.insert(*id, ast.clone());
                    }
                    StereoBondDelta::Add {
                        id,
                        site,
                        ligands,
                        ast,
                    } => {
                        added_stereo_bond.insert(*id, (*site, ligands.clone(), ast.clone()));
                    }
                    StereoBondDelta::ModifyField { id, .. }
                    | StereoBondDelta::ModifyConstraint { id, .. }
                    | StereoBondDelta::Apply { id, .. }
                    | StereoBondDelta::Swap { id, .. }
                    | StereoBondDelta::Mirror { id, .. } => {
                        stereo_bond_changes
                            .entry(*id)
                            .or_default()
                            .push(stereo.clone());
                    }
                },
                Delta::Constraint(ConstraintDelta::Add(c)) => added_constraints.push(c.clone()),
                Delta::Constraint(ConstraintDelta::Remove(c)) => {
                    removed_constraints.push(c.clone())
                }
            }
        }

        // Union node index keyed by atom id: lhs atoms keep their id, created atoms append.
        let mut atom_index: HashMap<AtomId, usize> =
            HashMap::with_capacity(atom_count + added_atoms.len());
        for node in 0..atom_count {
            atom_index.insert(AtomId(node as u32), node);
        }
        for (offset, &id) in added_atoms.keys().enumerate() {
            atom_index.insert(id, atom_count + offset);
        }
        // Union edge index keyed by bond id (for the stereo-bond site): lhs bonds keep their id,
        // created bonds append — same shape as `atom_index`.
        let mut bond_index: HashMap<BondId, usize> = HashMap::new();
        for edge in 0..bond_count {
            bond_index.insert(BondId(edge as u32), edge);
        }
        for (offset, &id) in added_bonds.keys().enumerate() {
            bond_index.insert(id, bond_count + offset);
        }

        let mut atoms: Vec<EntitySpan<AtomAst>> =
            Vec::with_capacity(atom_count + added_atoms.len());
        for node in 0..atom_count {
            let id = AtomId(node as u32);
            if let Some(ast) = removed_atoms.get(&id) {
                atoms.push(EntitySpan::Removed(ast.clone()));
            } else if let Some(changes) = atom_changes.get(&id) {
                let left = lhs.atom(id).ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_atom_change(&mut right, change)?;
                }
                atoms.push(EntitySpan::Modified {
                    lhs: left,
                    rhs: right,
                });
            } else {
                atoms.push(EntitySpan::Unchanged(lhs.atom(id).ast.clone()));
            }
        }
        for ast in added_atoms.into_values() {
            atoms.push(EntitySpan::Added(ast));
        }

        let mut bonds: Vec<EntitySpan<BondAst>> =
            Vec::with_capacity(bond_count + added_bonds.len());
        let mut edges: Vec<[u32; 2]> = Vec::with_capacity(bond_count + added_bonds.len());
        for edge in 0..bond_count {
            let id = BondId(edge as u32);
            let [a, b] = lhs.raw_graph().edge_endpoints(EdgeId(edge as u32));
            edges.push([a.0, b.0]);
            if let Some(ast) = removed_bonds.get(&id) {
                bonds.push(EntitySpan::Removed(ast.clone()));
            } else if let Some(changes) = bond_changes.get(&id) {
                let left = lhs.bond(id).ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_bond_change(&mut right, change)?;
                }
                bonds.push(EntitySpan::Modified {
                    lhs: left,
                    rhs: right,
                });
            } else {
                bonds.push(EntitySpan::Unchanged(lhs.bond(id).ast.clone()));
            }
        }
        for (atoms, ast) in added_bonds.into_values() {
            edges.push([atom_index[&atoms[0]] as u32, atom_index[&atoms[1]] as u32]);
            bonds.push(EntitySpan::Added(ast));
        }

        // Overlay columns: lhs overlays tagged by their fold (Removed/Modified/Unchanged),
        // created overlays appended; participants mapped to the union id space via `atom_index`.
        let mut aromatic_entries: Vec<(Vec<NodeId>, EntitySpan<AromaticSystemAst>)> = Vec::new();
        for view in lhs.aromatic_systems().iter() {
            let participants: Vec<NodeId> = view
                .atom_ids()
                .map(|a| NodeId(atom_index[&a] as u32))
                .collect();
            if let Some(ast) = removed_aromatic.get(&view.id) {
                aromatic_entries.push((participants, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = aromatic_changes.get(&view.id) {
                let left = view.ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_aromatic_change(&mut right, change)?;
                }
                aromatic_entries.push((
                    participants,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                aromatic_entries.push((participants, EntitySpan::Unchanged(view.ast.clone())));
            }
        }
        for (atoms, ast) in added_aromatic.into_values() {
            let participants: Vec<NodeId> =
                atoms.iter().map(|a| NodeId(atom_index[a] as u32)).collect();
            aromatic_entries.push((participants, EntitySpan::Added(ast)));
        }
        let aromatic_systems = VarRelationSet::new(aromatic_entries);

        let mut multicenter_entries: Vec<(Vec<NodeId>, EntitySpan<MulticenterBondAst>)> =
            Vec::new();
        for view in lhs.multicenter_bonds().iter() {
            let participants: Vec<NodeId> = view
                .atom_ids()
                .map(|a| NodeId(atom_index[&a] as u32))
                .collect();
            if let Some(ast) = removed_multicenter.get(&view.id) {
                multicenter_entries.push((participants, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = multicenter_changes.get(&view.id) {
                let left = view.ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_multicenter_change(&mut right, change)?;
                }
                multicenter_entries.push((
                    participants,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                multicenter_entries.push((participants, EntitySpan::Unchanged(view.ast.clone())));
            }
        }
        for (atoms, ast) in added_multicenter.into_values() {
            let participants: Vec<NodeId> =
                atoms.iter().map(|a| NodeId(atom_index[a] as u32)).collect();
            multicenter_entries.push((participants, EntitySpan::Added(ast)));
        }
        let multicenter_bonds = VarRelationSet::new(multicenter_entries);

        let mut noncovalent_entries: Vec<([NodeId; 2], EntitySpan<NoncovalentBondAst>)> =
            Vec::new();
        for view in lhs.noncovalent_bonds().iter() {
            let [a, b] = view.atom_ids();
            let participants = [NodeId(atom_index[&a] as u32), NodeId(atom_index[&b] as u32)];
            if let Some(ast) = removed_noncovalent.get(&view.id) {
                noncovalent_entries.push((participants, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = noncovalent_changes.get(&view.id) {
                let left = view.ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_noncovalent_change(&mut right, change)?;
                }
                noncovalent_entries.push((
                    participants,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                noncovalent_entries.push((participants, EntitySpan::Unchanged(view.ast.clone())));
            }
        }
        for ([a, b], ast) in added_noncovalent.into_values() {
            let participants = [NodeId(atom_index[&a] as u32), NodeId(atom_index[&b] as u32)];
            noncovalent_entries.push((participants, EntitySpan::Added(ast)));
        }
        let noncovalent_bonds = FixedRelationSet::new(noncovalent_entries);

        let mut dative_entries: Vec<([NodeId; 1], Vec<NodeId>, EntitySpan<DativeBondAst>)> =
            Vec::new();
        for view in lhs.dative_bonds().iter() {
            let acceptor = [NodeId(atom_index[&view.acceptor_id()] as u32)];
            let donors: Vec<NodeId> = view
                .donor_ids()
                .map(|a| NodeId(atom_index[&a] as u32))
                .collect();
            if let Some(ast) = removed_dative.get(&view.id) {
                dative_entries.push((acceptor, donors, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = dative_changes.get(&view.id) {
                let left = view.ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_dative_change(&mut right, change)?;
                }
                dative_entries.push((
                    acceptor,
                    donors,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                dative_entries.push((acceptor, donors, EntitySpan::Unchanged(view.ast.clone())));
            }
        }
        for (donors, acceptor, ast) in added_dative.into_values() {
            let acceptor = [NodeId(atom_index[&acceptor] as u32)];
            let donors: Vec<NodeId> = donors
                .iter()
                .map(|a| NodeId(atom_index[a] as u32))
                .collect();
            dative_entries.push((acceptor, donors, EntitySpan::Added(ast)));
        }
        let dative_bonds = FixedVarBirelationSet::new(dative_entries);

        // Stereo overlays: lhs entities tagged by their fold (Removed/Modified/Unchanged), created
        // ones appended. Site/ligand ids mapped to the union frame via `atom_index`/`bond_index`.
        let mut stereo_atom_entries: Vec<(
            [NodeId; 1],
            Vec<StereoLigand>,
            EntitySpan<StereoAtomAst>,
        )> = Vec::new();
        for view in lhs.stereo_atoms().iter() {
            let site = [NodeId::from(view.site_id())];
            let ligands = view.ligand_frame();
            if let Some(ast) = removed_stereo_atom.get(&view.id) {
                stereo_atom_entries.push((site, ligands, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = stereo_atom_changes.get(&view.id) {
                let left = view.ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_stereo_atom_change(&mut right, change)?;
                }
                stereo_atom_entries.push((
                    site,
                    ligands,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                stereo_atom_entries.push((site, ligands, EntitySpan::Unchanged(view.ast.clone())));
            }
        }
        for (site, ligands, ast) in added_stereo_atom.into_values() {
            let site = [NodeId(atom_index[&site] as u32)];
            let ligands: Vec<StereoLigand> = ligands
                .iter()
                .map(|l| StereoLigand::new(AtomId(atom_index[&l.atom_id] as u32), l.kind))
                .collect();
            stereo_atom_entries.push((site, ligands, EntitySpan::Added(ast)));
        }
        let stereo_atoms = FixedVarBirelationSet::new(stereo_atom_entries);

        let mut stereo_bond_entries: Vec<(
            [EdgeId; 1],
            Vec<StereoLigand>,
            EntitySpan<StereoBondAst>,
        )> = Vec::new();
        for view in lhs.stereo_bonds().iter() {
            let site = [EdgeId::from(view.site_id())];
            let ligands = view.ligand_frame();
            if let Some(ast) = removed_stereo_bond.get(&view.id) {
                stereo_bond_entries.push((site, ligands, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = stereo_bond_changes.get(&view.id) {
                let left = view.ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_stereo_bond_change(&mut right, change)?;
                }
                stereo_bond_entries.push((
                    site,
                    ligands,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                stereo_bond_entries.push((site, ligands, EntitySpan::Unchanged(view.ast.clone())));
            }
        }
        for (site, ligands, ast) in added_stereo_bond.into_values() {
            let site = [EdgeId(bond_index[&site] as u32)];
            let ligands: Vec<StereoLigand> = ligands
                .iter()
                .map(|l| StereoLigand::new(AtomId(atom_index[&l.atom_id] as u32), l.kind))
                .collect();
            stereo_bond_entries.push((site, ligands, EntitySpan::Added(ast)));
        }
        let stereo_bonds = FixedVarBirelationSet::new(stereo_bond_entries);

        let mut constraints: Vec<ConstraintSpan> =
            Vec::with_capacity(lhs.constraints().len() + added_constraints.len());
        for c in lhs.constraints().iter() {
            match removed_constraints.iter().position(|r| r == c) {
                Some(pos) => {
                    removed_constraints.remove(pos);
                    constraints.push(ConstraintSpan::Removed(c.clone()));
                }
                None => constraints.push(ConstraintSpan::Unchanged(c.clone())),
            }
        }
        for c in added_constraints {
            constraints.push(ConstraintSpan::Added(c));
        }

        let graph = Graph::new(atoms.len(), &edges);
        Ok(ReactionSpanAst::from_parts(
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
        ))
    }

    /// The reverse reaction: the product becomes the reactant and every delta is inverted and
    /// re-anchored to the product's (compacted) id space. `reverse().to_reaction_span()` swaps the
    /// sides of `self`'s span. `Err(Contradiction)` if the deltas are inconsistent.
    pub fn reverse(&self) -> Result<ReactionAst, Contradiction> {
        let deltas = self.deltas.clone().canonicalize()?;
        let new_lhs = self.to_reaction_span()?.rhs();
        let atom_count = self.lhs.atoms().count();
        let bond_count = self.lhs.bonds().count();

        let mut removed_atoms: Vec<AtomId> = Vec::new();
        let mut created_atoms: Vec<AtomId> = Vec::new();
        let mut removed_bonds: Vec<BondId> = Vec::new();
        let mut created_bonds: Vec<BondId> = Vec::new();
        let mut removed_dative: Vec<DativeBondId> = Vec::new();
        let mut created_dative: Vec<DativeBondId> = Vec::new();
        let mut removed_aromatic: Vec<AromaticSystemId> = Vec::new();
        let mut created_aromatic: Vec<AromaticSystemId> = Vec::new();
        let mut removed_multicenter: Vec<MulticenterBondId> = Vec::new();
        let mut created_multicenter: Vec<MulticenterBondId> = Vec::new();
        let mut removed_noncovalent: Vec<NoncovalentBondId> = Vec::new();
        let mut created_noncovalent: Vec<NoncovalentBondId> = Vec::new();
        let mut removed_stereo_atom: Vec<StereoAtomId> = Vec::new();
        let mut created_stereo_atom: Vec<StereoAtomId> = Vec::new();
        let mut removed_stereo_bond: Vec<StereoBondId> = Vec::new();
        let mut created_stereo_bond: Vec<StereoBondId> = Vec::new();
        for delta in deltas.iter() {
            match delta {
                Delta::Atom(AtomDelta::Remove { id, .. }) => removed_atoms.push(*id),
                Delta::Atom(AtomDelta::Add { id, .. }) => created_atoms.push(*id),
                Delta::Bond(BondDelta::Remove { id, .. }) => removed_bonds.push(*id),
                Delta::Bond(BondDelta::Add { id, .. }) => created_bonds.push(*id),
                Delta::DativeBond(DativeBondDelta::Remove { id, .. }) => removed_dative.push(*id),
                Delta::DativeBond(DativeBondDelta::Add { id, .. }) => created_dative.push(*id),
                Delta::AromaticSystem(AromaticSystemDelta::Remove { id, .. }) => {
                    removed_aromatic.push(*id)
                }
                Delta::AromaticSystem(AromaticSystemDelta::Add { id, .. }) => {
                    created_aromatic.push(*id)
                }
                Delta::MulticenterBond(MulticenterBondDelta::Remove { id, .. }) => {
                    removed_multicenter.push(*id)
                }
                Delta::MulticenterBond(MulticenterBondDelta::Add { id, .. }) => {
                    created_multicenter.push(*id)
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, .. }) => {
                    removed_noncovalent.push(*id)
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Add { id, .. }) => {
                    created_noncovalent.push(*id)
                }
                Delta::StereoAtom(StereoAtomDelta::Remove { id, .. }) => {
                    removed_stereo_atom.push(*id)
                }
                Delta::StereoAtom(StereoAtomDelta::Add { id, .. }) => created_stereo_atom.push(*id),
                Delta::StereoBond(StereoBondDelta::Remove { id, .. }) => {
                    removed_stereo_bond.push(*id)
                }
                Delta::StereoBond(StereoBondDelta::Add { id, .. }) => created_stereo_bond.push(*id),
                _ => {}
            }
        }

        // Forward → reverse id-space maps, matching `rhs()`'s compaction: survivors take ids in
        // union order (lhs in place, created appended); deleted entities become created in the
        // reverse and take fresh ids after the survivors.
        let remapping = IdRemapping::new(
            reversed_remapping(atom_count, &removed_atoms, &created_atoms),
            reversed_remapping(bond_count, &removed_bonds, &created_bonds),
            reversed_remapping(
                self.lhs.dative_bonds().count(),
                &removed_dative,
                &created_dative,
            ),
            reversed_remapping(
                self.lhs.aromatic_systems().count(),
                &removed_aromatic,
                &created_aromatic,
            ),
            reversed_remapping(
                self.lhs.multicenter_bonds().count(),
                &removed_multicenter,
                &created_multicenter,
            ),
            reversed_remapping(
                self.lhs.noncovalent_bonds().count(),
                &removed_noncovalent,
                &created_noncovalent,
            ),
            reversed_remapping(
                self.lhs.stereo_atoms().count(),
                &removed_stereo_atom,
                &created_stereo_atom,
            ),
            reversed_remapping(
                self.lhs.stereo_bonds().count(),
                &removed_stereo_bond,
                &created_stereo_bond,
            ),
        );

        let reversed: Deltas = deltas
            .iter()
            .map(|delta| remap_delta(delta.clone().inverse(), &remapping))
            .collect();
        Ok(ReactionAst::new(new_lhs, reversed))
    }
}

/// Build a forward → reverse id-space map for one entity kind: surviving lhs ids (those not in
/// `removed`, in id order) then `created` (sorted) take reverse ids `0..k`; `removed` ids (which
/// become created in the reverse) take fresh ids after the survivors.
fn reversed_remapping<Id>(lhs_count: usize, removed: &[Id], created: &[Id]) -> HashMap<Id, Id>
where
    Id: Copy + Eq + Hash + Ord + From<usize>,
{
    let removed_set: HashSet<Id> = removed.iter().copied().collect();
    let mut created: Vec<Id> = created.to_vec();
    created.sort_unstable();
    (0..lhs_count)
        .map(Id::from)
        .filter(|id| !removed_set.contains(id))
        .chain(created)
        .chain(removed.iter().copied())
        .enumerate()
        .map(|(rev, id)| (id, Id::from(rev)))
        .collect()
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::constraint::{Constraint, Constraints, MoleculeConstraint};
    use super::super::delta::Deltas;
    use super::super::edit::{BondFieldChange, NoncovalentBondFieldChange, StereoAtomFieldChange};
    use super::super::ligand::{StereoLigand, StereoLigandKind};
    use super::super::noncovalent::{NoncovalentBondKind, NoncovalentBondKindAst};
    use super::super::stereo::{StereoConfigurationAst, StereoCoset, StereoKind};
    use super::super::value::ValueAst;
    use super::*;

    #[rstest]
    fn test_reaction_ast_to_reaction_span() {
        let reaction = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::C),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new: ValueAst::Lit(2),
                },
            })]),
        );
        let span = reaction.to_reaction_span().unwrap();
        assert_eq!(
            span.atoms(),
            [
                EntitySpan::Unchanged(AtomAst::from_element(Element::C)),
                EntitySpan::Unchanged(AtomAst::from_element(Element::C)),
            ],
        );
        assert_eq!(
            span.bonds(),
            [EntitySpan::Modified {
                lhs: BondAst::from_order(1),
                rhs: BondAst::from_order(2),
            }],
        );
    }

    #[rstest]
    fn test_reaction_span_ast_right(substitution_reaction: ReactionAst) {
        let span = substitution_reaction.to_reaction_span().unwrap();
        assert_eq!(
            span.atoms(),
            [
                EntitySpan::Unchanged(AtomAst::from_element(Element::C)),
                EntitySpan::Removed(AtomAst::from_element(Element::O)),
                EntitySpan::Added(AtomAst::from_element(Element::N)),
            ],
        );
        assert_eq!(
            span.bonds(),
            [
                EntitySpan::Removed(BondAst::from_order(1)),
                EntitySpan::Added(BondAst::from_order(1)),
            ],
        );
        assert_eq!(
            span.rhs(),
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::N),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                ..Default::default()
            }),
        );
    }

    #[rstest]
    fn test_reaction_span_ast_left(substitution_reaction: ReactionAst) {
        let span = substitution_reaction.to_reaction_span().unwrap();
        assert_eq!(
            span.lhs(),
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                ..Default::default()
            }),
        );
    }

    #[rstest]
    #[case::order_change(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::C),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new: ValueAst::Lit(2),
                },
            })]),
        ),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ..Default::default()
        }),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        }),
    )]
    #[case::substitution(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    ast: AtomAst::from_element(Element::N),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(1),
                    atoms: [AtomId(0), AtomId(2)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        ),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        }),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        }),
    )]
    #[case::stereo_atom(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::F),
                    AtomAst::from_element(Element::Cl),
                    AtomAst::from_element(Element::Br),
                    AtomAst::from_element(Element::I),
                ],
                stereo_atoms: vec![(
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                )],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Remove {
                id: StereoAtomId(0),
                site: AtomId(0),
                ligands: vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                ast: StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            })]),
        ),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::F),
                AtomAst::from_element(Element::Cl),
                AtomAst::from_element(Element::Br),
                AtomAst::from_element(Element::I),
            ],
            bonds: vec![],
            ..Default::default()
        }),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::F),
                AtomAst::from_element(Element::Cl),
                AtomAst::from_element(Element::Br),
                AtomAst::from_element(Element::I),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        }),
    )]
    fn test_reaction_ast_reverse(
        #[case] forward: ReactionAst,
        #[case] expected_reactant: MoleculeAst,
        #[case] expected_product: MoleculeAst,
    ) {
        // The reverse reaction's reactant is the forward product; its product is the forward
        // reactant.
        let span = forward.reverse().unwrap().to_reaction_span().unwrap();
        assert_eq!(span.lhs(), expected_reactant);
        assert_eq!(span.rhs(), expected_product);
    }

    #[rstest]
    #[case::identity(
        3,
        vec![],
        vec![],
        HashMap::from([(AtomId(0), AtomId(0)), (AtomId(1), AtomId(1)), (AtomId(2), AtomId(2))])
    )]
    #[case::removed_compacted(
        4,
        vec![AtomId(1), AtomId(3)],
        vec![],
        HashMap::from([
            (AtomId(0), AtomId(0)),
            (AtomId(2), AtomId(1)),
            (AtomId(1), AtomId(2)),
            (AtomId(3), AtomId(3)),
        ])
    )]
    #[case::created_appended_sorted(
        2,
        vec![],
        vec![AtomId(5), AtomId(4)],
        HashMap::from([
            (AtomId(0), AtomId(0)),
            (AtomId(1), AtomId(1)),
            (AtomId(4), AtomId(2)),
            (AtomId(5), AtomId(3)),
        ])
    )]
    #[case::removed_and_created(
        3,
        vec![AtomId(1)],
        vec![AtomId(7)],
        HashMap::from([
            (AtomId(0), AtomId(0)),
            (AtomId(2), AtomId(1)),
            (AtomId(7), AtomId(2)),
            (AtomId(1), AtomId(3)),
        ])
    )]
    fn test_reversed_remapping(
        #[case] lhs_count: usize,
        #[case] removed: Vec<AtomId>,
        #[case] created: Vec<AtomId>,
        #[case] expected: HashMap<AtomId, AtomId>,
    ) {
        assert_eq!(reversed_remapping(lhs_count, &removed, &created), expected);
    }

    #[rstest]
    #[case::unchanged(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                constraints: Constraints::from(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })),
                ..Default::default()
            }),
            Deltas::new(),
        ),
        vec![ConstraintSpan::Unchanged(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }))],
    )]
    #[case::added(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(
                Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
            ))]),
        ),
        vec![ConstraintSpan::Added(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }))],
    )]
    #[case::removed(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                constraints: Constraints::from(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Remove(
                Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
            ))]),
        ),
        vec![ConstraintSpan::Removed(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }))],
    )]
    fn test_reaction_span_ast_constraints(
        #[case] reaction: ReactionAst,
        #[case] expected: Vec<ConstraintSpan>,
    ) {
        assert_eq!(
            reaction.to_reaction_span().unwrap().constraints(),
            expected.as_slice(),
        );
    }

    #[rstest]
    #[case::add(ReactionAst::new(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(
            Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
        ))]),
    ))]
    #[case::remove(ReactionAst::new(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            constraints: Constraints::from(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })),
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Constraint(ConstraintDelta::Remove(
            Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
        ))]),
    ))]
    #[case::dative_add(ReactionAst::new(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::B),
                AtomAst::from_element(Element::N),
            ],
            bonds: vec![],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::DativeBond(DativeBondDelta::Add {
            id: DativeBondId(0),
            donors: vec![AtomId(0), AtomId(2)],
            acceptor: AtomId(1),
            ast: DativeBondAst::from_order(1),
        })]),
    ))]
    #[case::aromatic_add(ReactionAst::new(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)],
            bonds: vec![],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Add {
            id: AromaticSystemId(0),
            atoms: vec![AtomId(0), AtomId(1)],
            ast: AromaticSystemAst::from_electrons(vec![1, 2]),
        })]),
    ))]
    #[case::multicenter_add(ReactionAst::new(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::B),
                AtomAst::from_element(Element::H),
                AtomAst::from_element(Element::B),
            ],
            bonds: vec![],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::MulticenterBond(MulticenterBondDelta::Add {
            id: MulticenterBondId(0),
            atoms: vec![AtomId(0), AtomId(1), AtomId(2)],
            ast: MulticenterBondAst::from_electrons(vec![3, 5, 7]),
        })]),
    ))]
    #[case::noncovalent_add(ReactionAst::new(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
            bonds: vec![],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(0),
            atoms: [AtomId(0), AtomId(1)],
            ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        })]),
    ))]
    #[case::noncovalent_remove(ReactionAst::new(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
            constraints: Constraints::new(),
            ..Default::default()
        }),
        Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
            id: NoncovalentBondId(0),
            atoms: [AtomId(0), AtomId(1)],
            ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        })]),
    ))]
    #[case::noncovalent_modify(ReactionAst::new(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
            constraints: Constraints::new(),
            ..Default::default()
        }),
        Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
            id: NoncovalentBondId(0),
            change: NoncovalentBondFieldChange::Kind {
                old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
                new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
            },
        })]),
    ))]
    #[case::stereo_atom_add(ReactionAst::new(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::F),
                AtomAst::from_element(Element::Cl),
                AtomAst::from_element(Element::Br),
                AtomAst::from_element(Element::I),
            ],
            bonds: vec![],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Add {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            ast: StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        })]),
    ))]
    #[case::stereo_atom_remove(ReactionAst::new(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::F),
                AtomAst::from_element(Element::Cl),
                AtomAst::from_element(Element::Br),
                AtomAst::from_element(Element::I),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        }),
        Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Remove {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            ast: StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        })]),
    ))]
    #[case::stereo_atom_modify(ReactionAst::new(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::F),
                AtomAst::from_element(Element::Cl),
                AtomAst::from_element(Element::Br),
                AtomAst::from_element(Element::I),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        }),
        Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
            id: StereoAtomId(0),
            change: StereoAtomFieldChange::Configuration {
                old: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                new: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            },
        })]),
    ))]
    #[case::stereo_bond_add(ReactionAst::new(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(2)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::StereoBond(StereoBondDelta::Add {
            id: StereoBondId(0),
            site: BondId(1),
            ligands: vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            ast: StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        })]),
    ))]
    fn test_reaction_span_ast_to_reaction(#[case] reaction: ReactionAst) {
        assert_eq!(
            reaction.clone().to_reaction_span().unwrap().to_reaction(),
            reaction,
        );
    }

    // A constraint naming atom 1 survives the left projection (atom 1 present) and is dropped on the
    // right (atom 1 removed) — the projection compaction drops the dangling ref.
    #[rstest]
    fn test_reaction_span_ast_project() {
        let span = ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                constraints: Constraints::from(Constraint::Molecule(
                    MoleculeConstraint::Connected {
                        atoms: Some(vec![AtomId(0), AtomId(1)]),
                    },
                )),
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        )
        .to_reaction_span()
        .unwrap();
        assert_eq!(
            span.lhs().constraints().as_slice(),
            [Constraint::Molecule(MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(0), AtomId(1)]),
            })],
        );
        assert!(span.rhs().constraints().is_empty());
    }

    #[rstest]
    #[case::unchanged(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
                noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::new(),
        ),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
            constraints: Constraints::new(),
            ..Default::default()
        }),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
            constraints: Constraints::new(),
            ..Default::default()
        }),
    )]
    #[case::added(
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
                bonds: vec![],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                id: NoncovalentBondId(0),
                atoms: [AtomId(0), AtomId(1)],
                ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            })]),
        ),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
            bonds: vec![],
            ..Default::default()
        }),
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
            constraints: Constraints::new(),
            ..Default::default()
        }),
    )]
    fn test_reaction_span_ast_project_overlay(
        #[case] reaction: ReactionAst,
        #[case] expected_left: MoleculeAst,
        #[case] expected_right: MoleculeAst,
    ) {
        // An unchanged overlay carries through both projections; an added overlay is absent on the
        // left, present on the right.
        let span = reaction.to_reaction_span().unwrap();
        assert_eq!(span.lhs(), expected_left);
        assert_eq!(span.rhs(), expected_right);
    }

    // C-O with atom 1 (O) and its bond removed, replaced by a new N (atom 2) bonded to C.
    #[fixture]
    fn substitution_reaction() -> ReactionAst {
        ReactionAst::new(
            MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    ast: AtomAst::from_element(Element::N),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(1),
                    atoms: [AtomId(0), AtomId(2)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        )
    }

    #[rstest]
    fn test_reaction_span_ast_superimpose() {
        // atom 0 C stays C (Unchanged); atom 1 C→N (Modified); atom 2 O and its bond are Added.
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
            ],
            ..Default::default()
        });
        let atoms = Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 3);
        let correspondence = MoleculeCorrespondence::induce(&left, &right, atoms);

        let span = ReactionSpanAst::superimpose(&left, &right, &correspondence);

        assert_eq!(span.lhs(), left);
        assert_eq!(span.rhs(), right);
    }

    #[rstest]
    fn test_reaction_span_ast_correspondence() {
        // atom 0 unchanged, 1 modified (C→N), 2 removed (left) with 2 added (right O): all four
        // EntitySpan variants in the atom column.
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 3],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });
        let atoms = Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 3, 3);
        let correspondence = MoleculeCorrespondence::induce(&left, &right, atoms);
        let span = ReactionSpanAst::superimpose(&left, &right, &correspondence);

        // recovers the input correspondence, and inverts `superimpose`.
        assert_eq!(span.correspondence(), correspondence);
        assert_eq!(
            ReactionSpanAst::superimpose(&span.lhs(), &span.rhs(), &span.correspondence()),
            span
        );
    }

    #[rstest]
    fn test_molecule_ast_difference_to() {
        // C-C (order 1) → C-C (order 2), total correspondence: a single bond-order modify.
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(2))],
            ..Default::default()
        });
        let atoms = Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 2, 2);
        let correspondence = MoleculeCorrespondence::induce(&left, &right, atoms);
        assert_eq!(
            left.difference_to(&right, &correspondence),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new: ValueAst::Lit(2),
                },
            })]),
        );
    }

    #[rstest]
    fn test_reaction_span_ast_superimpose_removed() {
        // left has a third atom, a bond to it, and a dative onto it — all left-exposed (Removed).
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 3],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
            ],
            dative: vec![(vec![AtomId(2)], AtomId(1), DativeBondAst::from_order(1))],
            ..Default::default()
        });
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        });
        let atoms = Correspondence::new(vec![(NodeId(0), NodeId(0)), (NodeId(1), NodeId(1))], 3, 2);
        let correspondence = MoleculeCorrespondence::induce(&left, &right, atoms);

        let span = ReactionSpanAst::superimpose(&left, &right, &correspondence);

        assert_eq!(span.lhs(), left);
        assert_eq!(span.rhs(), right);
    }
}
