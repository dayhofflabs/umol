//! Reaction span AST: the superimposed `L ∪_K R` graph encoding a reaction's DPO rule span.
//!
//! Materialized superimposed graph carrying, per atom/bond, both its before and after state plus a
//! membership tag. The DPO span `L ←K─ R` is read off the tags — `K = Unchanged ∪ Modified`,
//! `L = K ∪ Removed`, `R = K ∪ Added` — and `rhs()` / `lhs()` project the two sides back to
//! a `MoleculeAst`. `Modified` (a preserved entity relabeled across the reaction) is the
//! relabeling-DPO reading: the entity persists in `K`, its label resolved per side.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::Hash;

use thiserror::Error;
use umol_graph_core::{
    Correspondence, EdgeId, FixedRelationSet, FixedVarBirelationSet, Graph, NodeId, Ordered,
    RelationId, Remapping, Unordered, VarRelationSet,
};

use super::aromatic::AromaticSystemForm;
use super::atom::AtomForm;
use super::bond::BondForm;
use super::constraint::Constraint;
use super::correspondence::MoleculeCorrespondence;
use super::dative::DativeBondForm;
use super::delta::{
    apply_aromatic_change, apply_atom_change, apply_bond_change, apply_dative_change,
    apply_multicenter_change, apply_noncovalent_change, apply_stereo_atom_change,
    apply_stereo_bond_change, remap_delta, AromaticSystemDelta, AtomDelta, BondDelta,
    ConstraintDelta, ConstraintSpan, DativeBondDelta, Delta, Deltas, EntityFold, EntitySpan,
    MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta,
};
use super::entity::Entity;
use super::error::Contradiction;
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigand;
use super::molecule::{
    validate_constraint_references, validate_entry_references, MoleculeAst, MoleculeEntries,
    MoleculeEntriesError,
};
use super::multicenter::MulticenterBondForm;
use super::noncovalent::NoncovalentBondForm;
use super::reaction::ReactionAst;
use super::remap::IdRemapping;
use super::stereo::{StereoAtomForm, StereoBondForm};
use super::traits::{Canonicalize, EntityPatch};

/// The superimposed reaction graph — the reaction's DPO rule span, materialized. The union
/// topology is the `lhs` id space (deleted entities kept as nodes/edges) with created entities
/// appended; `atoms` / `bonds` are indexed parallel to the graph's nodes / edges.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionSpanAst {
    graph: Graph,
    atoms: Vec<EntitySpan<AtomForm>>,
    bonds: Vec<EntitySpan<BondForm>>,
    dative_bonds:
        FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, EntitySpan<DativeBondForm>>,
    aromatic_systems: VarRelationSet<NodeId, Unordered, EntitySpan<AromaticSystemForm>>,
    multicenter_bonds: VarRelationSet<NodeId, Unordered, EntitySpan<MulticenterBondForm>>,
    noncovalent_bonds: FixedRelationSet<NodeId, Unordered, EntitySpan<NoncovalentBondForm>, 2>,
    stereo_atoms: FixedVarBirelationSet<
        NodeId,
        Ordered,
        1,
        StereoLigand,
        Ordered,
        EntitySpan<StereoAtomForm>,
    >,
    stereo_bonds: FixedVarBirelationSet<
        EdgeId,
        Ordered,
        1,
        StereoLigand,
        Ordered,
        EntitySpan<StereoBondForm>,
    >,
    constraints: Vec<ConstraintSpan>,
}

/// Flat constructor input for [`ReactionSpanAst::from_entries`]. Each [`EntitySpan`] is present on
/// at least one side by construction; a value absent from both sides has no entry representation.
#[derive(Clone, Debug, Default)]
pub struct ReactionSpanEntries {
    pub atoms: Vec<EntitySpan<AtomForm>>,
    pub bonds: Vec<(AtomId, AtomId, EntitySpan<BondForm>)>,
    pub dative: Vec<(Vec<AtomId>, AtomId, EntitySpan<DativeBondForm>)>,
    pub aromatic: Vec<(Vec<AtomId>, EntitySpan<AromaticSystemForm>)>,
    pub multicenter: Vec<(Vec<AtomId>, EntitySpan<MulticenterBondForm>)>,
    pub noncovalent: Vec<(AtomId, AtomId, EntitySpan<NoncovalentBondForm>)>,
    pub stereo_atoms: Vec<(AtomId, Vec<StereoLigand>, EntitySpan<StereoAtomForm>)>,
    pub stereo_bonds: Vec<(BondId, Vec<StereoLigand>, EntitySpan<StereoBondForm>)>,
    pub constraints: Vec<ConstraintSpan>,
}

/// Failure to construct a reaction span from structurally inconsistent entries.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionSpanEntriesError {
    #[error("reaction span entries reference unavailable {entity}")]
    InvalidReference { entity: Entity },
}

impl ReactionSpanAst {
    /// Construct a reaction span from entries whose structural integrity is established by the
    /// caller.
    ///
    /// # Panics
    ///
    /// Panics when a participant, site, ligand, or constraint references an entity absent from the
    /// union frame or from a side on which the referring entry is present.
    pub fn from_entries(entries: ReactionSpanEntries) -> Self {
        Self::try_from_entries(entries)
            .unwrap_or_else(|error| panic!("invalid reaction span entries: {error}"))
    }

    /// Construct a reaction span after checking union-frame and projected-side references.
    ///
    /// # Errors
    ///
    /// Returns [`ReactionSpanEntriesError::InvalidReference`] when a participant, site, ligand, or
    /// constraint references an entity absent from the union frame or from a side on which the
    /// referring entry is present. Chemistry and other semantic properties are not validated.
    pub fn try_from_entries(
        mut entries: ReactionSpanEntries,
    ) -> Result<Self, ReactionSpanEntriesError> {
        normalize_reaction_span_entries(&mut entries);
        validate_reaction_span_entries(&entries)?;

        let ReactionSpanEntries {
            atoms,
            bonds,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            stereo_atoms,
            stereo_bonds,
            constraints,
        } = entries;
        let edges: Vec<[u32; 2]> = bonds
            .iter()
            .map(|(first, second, _)| [first.0, second.0])
            .collect();
        let bond_values = bonds.into_iter().map(|(_, _, span)| span).collect();
        let dative_bonds = FixedVarBirelationSet::new(
            dative
                .into_iter()
                .map(|(donors, acceptor, span)| {
                    (
                        [NodeId::from(acceptor)],
                        donors.into_iter().map(NodeId::from).collect(),
                        span,
                    )
                })
                .collect(),
        );
        let aromatic_systems = VarRelationSet::new(
            aromatic
                .into_iter()
                .map(|(atoms, span)| (atoms.into_iter().map(NodeId::from).collect(), span))
                .collect(),
        );
        let multicenter_bonds = VarRelationSet::new(
            multicenter
                .into_iter()
                .map(|(atoms, span)| (atoms.into_iter().map(NodeId::from).collect(), span))
                .collect(),
        );
        let noncovalent_bonds = FixedRelationSet::new(
            noncovalent
                .into_iter()
                .map(|(first, second, span)| ([NodeId::from(first), NodeId::from(second)], span))
                .collect(),
        );
        let stereo_atoms = FixedVarBirelationSet::new(
            stereo_atoms
                .into_iter()
                .map(|(site, ligands, span)| ([NodeId::from(site)], ligands, span))
                .collect(),
        );
        let stereo_bonds = FixedVarBirelationSet::new(
            stereo_bonds
                .into_iter()
                .map(|(site, ligands, span)| ([EdgeId::from(site)], ligands, span))
                .collect(),
        );

        let span = Self {
            graph: Graph::new(atoms.len(), &edges),
            atoms,
            bonds: bond_values,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            stereo_atoms,
            stereo_bonds,
            constraints,
        };
        for side in [Side::Left, Side::Right] {
            validate_entry_references(&span.project_entries(side)).map_err(
                |error| match error {
                    MoleculeEntriesError::InvalidReference { entity } => {
                        ReactionSpanEntriesError::InvalidReference { entity }
                    }
                },
            )?;
        }
        Ok(span)
    }
}

/// R id → union id for one entity family: a matched right id reuses its left partner; a
/// right-unmatched id is appended after the left family's ids.
fn union_map<Id: Copy + Ord + Hash + From<usize>>(
    correspondence: &Correspondence<Id>,
    left_count: usize,
) -> HashMap<Id, Id> {
    let mut map: HashMap<Id, Id> = correspondence
        .matched_pairs()
        .iter()
        .map(|&(l, r)| (r, l))
        .collect();
    for (offset, &r) in correspondence.right_unmatched().iter().enumerate() {
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
    let mut matched_pairs = Vec::new();
    let (mut left, mut right) = (0usize, 0usize);
    for span in column {
        match span {
            EntitySpan::Unchanged(_) | EntitySpan::Modified { .. } => {
                matched_pairs.push((Id::from(left), Id::from(right)));
                left += 1;
                right += 1;
            }
            EntitySpan::Removed(_) => left += 1,
            EntitySpan::Added(_) => right += 1,
        }
    }
    Correspondence::new(matched_pairs, left, right)
        .unwrap_or_else(|_| panic!("recovered entity spans preserve partial-bijection invariants"))
}

impl MoleculeAst {
    /// The deltas transforming `self` (`L`) into `rhs` (`R`) under the per-entity correspondence
    /// `correspondence`: superimpose the two sides into a span, then read off its operational
    /// reaction. Returns `None` when the correspondence is not compatible with the supplied
    /// molecules.
    ///
    /// # Semantic properties
    ///
    /// Applying a returned delta collection to `self` reconstructs `rhs` in the lhs-anchored
    /// reaction frame: preserved entities retain their lhs ids and rhs-only entities are appended.
    /// The result is equivalent to `rhs` under the induced total correspondence, but need not be
    /// structurally equal when the supplied correspondence changes entity order.
    pub fn difference_to(
        &self,
        rhs: &MoleculeAst,
        correspondence: &MoleculeCorrespondence,
    ) -> Option<Deltas> {
        Some(
            ReactionSpanAst::superimpose(self, rhs, correspondence)?
                .to_reaction()
                .deltas,
        )
    }
}

fn correspondence_is_compatible(
    lhs: &MoleculeAst,
    rhs: &MoleculeAst,
    correspondence: &MoleculeCorrespondence,
) -> bool {
    let counts = [
        (
            correspondence.atoms().left_count(),
            lhs.atoms().count(),
            correspondence.atoms().right_count(),
            rhs.atoms().count(),
        ),
        (
            correspondence.bonds().left_count(),
            lhs.bonds().count(),
            correspondence.bonds().right_count(),
            rhs.bonds().count(),
        ),
        (
            correspondence.dative_bonds().left_count(),
            lhs.dative_bonds().count(),
            correspondence.dative_bonds().right_count(),
            rhs.dative_bonds().count(),
        ),
        (
            correspondence.aromatic_systems().left_count(),
            lhs.aromatic_systems().count(),
            correspondence.aromatic_systems().right_count(),
            rhs.aromatic_systems().count(),
        ),
        (
            correspondence.multicenter_bonds().left_count(),
            lhs.multicenter_bonds().count(),
            correspondence.multicenter_bonds().right_count(),
            rhs.multicenter_bonds().count(),
        ),
        (
            correspondence.noncovalent_bonds().left_count(),
            lhs.noncovalent_bonds().count(),
            correspondence.noncovalent_bonds().right_count(),
            rhs.noncovalent_bonds().count(),
        ),
        (
            correspondence.stereo_atoms().left_count(),
            lhs.stereo_atoms().count(),
            correspondence.stereo_atoms().right_count(),
            rhs.stereo_atoms().count(),
        ),
        (
            correspondence.stereo_bonds().left_count(),
            lhs.stereo_bonds().count(),
            correspondence.stereo_bonds().right_count(),
            rhs.stereo_bonds().count(),
        ),
    ];
    if counts.into_iter().any(
        |(declared_left, actual_left, declared_right, actual_right)| {
            declared_left != actual_left || declared_right != actual_right
        },
    ) {
        return false;
    }

    let atoms = correspondence.atoms();
    let same_atom_set = |left: Vec<AtomId>, mut right: Vec<AtomId>| {
        let Some(mut mapped): Option<Vec<_>> =
            left.into_iter().map(|atom| atoms.right_of(atom)).collect()
        else {
            return false;
        };
        mapped.sort_unstable();
        right.sort_unstable();
        mapped == right
    };
    let same_ligand_set = |left: Vec<StereoLigand>, mut right: Vec<StereoLigand>| {
        let Some(mut mapped): Option<Vec<_>> = left
            .into_iter()
            .map(|ligand| {
                atoms
                    .right_of(ligand.atom_id)
                    .map(|atom| StereoLigand::new(atom, ligand.kind))
            })
            .collect()
        else {
            return false;
        };
        mapped.sort_unstable();
        right.sort_unstable();
        mapped == right
    };

    if !correspondence
        .bonds()
        .matched_pairs()
        .iter()
        .all(|&(left, right)| {
            same_atom_set(
                lhs.bond(left).atom_ids().to_vec(),
                rhs.bond(right).atom_ids().to_vec(),
            )
        })
    {
        return false;
    }
    if !correspondence
        .dative_bonds()
        .matched_pairs()
        .iter()
        .all(|&(left, right)| {
            let lhs = lhs.dative_bond(left);
            let rhs = rhs.dative_bond(right);
            atoms.right_of(lhs.acceptor_id()) == Some(rhs.acceptor_id())
                && same_atom_set(lhs.donor_ids().collect(), rhs.donor_ids().collect())
        })
    {
        return false;
    }
    if !correspondence
        .aromatic_systems()
        .matched_pairs()
        .iter()
        .all(|&(left, right)| {
            same_atom_set(
                lhs.aromatic_system(left).atom_ids().collect(),
                rhs.aromatic_system(right).atom_ids().collect(),
            )
        })
    {
        return false;
    }
    if !correspondence
        .multicenter_bonds()
        .matched_pairs()
        .iter()
        .all(|&(left, right)| {
            same_atom_set(
                lhs.multicenter_bond(left).atom_ids().collect(),
                rhs.multicenter_bond(right).atom_ids().collect(),
            )
        })
    {
        return false;
    }
    if !correspondence
        .noncovalent_bonds()
        .matched_pairs()
        .iter()
        .all(|&(left, right)| {
            same_atom_set(
                lhs.noncovalent_bond(left).atom_ids().to_vec(),
                rhs.noncovalent_bond(right).atom_ids().to_vec(),
            )
        })
    {
        return false;
    }
    if !correspondence
        .stereo_atoms()
        .matched_pairs()
        .iter()
        .all(|&(left, right)| {
            let lhs = lhs.stereo_atom(left);
            let rhs = rhs.stereo_atom(right);
            atoms.right_of(lhs.site_id()) == Some(rhs.site_id())
                && same_ligand_set(lhs.ligand_frame(), rhs.ligand_frame())
        })
    {
        return false;
    }
    correspondence
        .stereo_bonds()
        .matched_pairs()
        .iter()
        .all(|&(left, right)| {
            let lhs = lhs.stereo_bond(left);
            let rhs = rhs.stereo_bond(right);
            correspondence.bonds().right_of(lhs.site_id()) == Some(rhs.site_id())
                && same_ligand_set(lhs.ligand_frame(), rhs.ligand_frame())
        })
}

fn contains_entry(entries: &ReactionSpanEntries, entity: Entity) -> bool {
    match entity {
        Entity::Atom(id) => id.index() < entries.atoms.len(),
        Entity::Bond(id) => id.index() < entries.bonds.len(),
        Entity::DativeBond(id) => id.index() < entries.dative.len(),
        Entity::AromaticSystem(id) => id.index() < entries.aromatic.len(),
        Entity::MulticenterBond(id) => id.index() < entries.multicenter.len(),
        Entity::NoncovalentBond(id) => id.index() < entries.noncovalent.len(),
        Entity::StereoAtom(id) => id.index() < entries.stereo_atoms.len(),
        Entity::StereoBond(id) => id.index() < entries.stereo_bonds.len(),
    }
}

fn validate_reaction_span_entries(
    entries: &ReactionSpanEntries,
) -> Result<(), ReactionSpanEntriesError> {
    let validate = |entity| {
        contains_entry(entries, entity)
            .then_some(())
            .ok_or(ReactionSpanEntriesError::InvalidReference { entity })
    };

    for (first, second, _) in &entries.bonds {
        validate(Entity::Atom(*first))?;
        validate(Entity::Atom(*second))?;
    }
    for (donors, acceptor, _) in &entries.dative {
        validate(Entity::Atom(*acceptor))?;
        for &donor in donors {
            validate(Entity::Atom(donor))?;
        }
    }
    for (atoms, _) in &entries.aromatic {
        for &atom in atoms {
            validate(Entity::Atom(atom))?;
        }
    }
    for (atoms, _) in &entries.multicenter {
        for &atom in atoms {
            validate(Entity::Atom(atom))?;
        }
    }
    for (first, second, _) in &entries.noncovalent {
        validate(Entity::Atom(*first))?;
        validate(Entity::Atom(*second))?;
    }
    for (site, ligands, _) in &entries.stereo_atoms {
        validate(Entity::Atom(*site))?;
        for ligand in ligands {
            validate(Entity::Atom(ligand.atom_id))?;
        }
    }
    for (site, ligands, _) in &entries.stereo_bonds {
        validate(Entity::Bond(*site))?;
        for ligand in ligands {
            validate(Entity::Atom(ligand.atom_id))?;
        }
    }

    let contains = |entity| contains_entry(entries, entity);
    for span in &entries.constraints {
        for constraint in [span.lhs(), span.rhs()].into_iter().flatten() {
            validate_constraint_references(constraint, &contains)
                .map_err(|entity| ReactionSpanEntriesError::InvalidReference { entity })?;
        }
    }
    Ok(())
}

fn normalize_entity_span<T: Canonicalize + Clone>(span: &mut EntitySpan<T>) {
    let unchanged = match span {
        EntitySpan::Modified { lhs, rhs } if lhs.canonical_eq(rhs) => Some(lhs.clone()),
        _ => None,
    };
    if let Some(value) = unchanged {
        *span = EntitySpan::Unchanged(value);
    }
}

fn normalize_reaction_span_entries(entries: &mut ReactionSpanEntries) {
    for span in &mut entries.atoms {
        normalize_entity_span(span);
    }
    for (_, _, span) in &mut entries.bonds {
        normalize_entity_span(span);
    }
    for (_, _, span) in &mut entries.dative {
        normalize_entity_span(span);
    }
    for (_, span) in &mut entries.aromatic {
        normalize_entity_span(span);
    }
    for (_, span) in &mut entries.multicenter {
        normalize_entity_span(span);
    }
    for (_, _, span) in &mut entries.noncovalent {
        normalize_entity_span(span);
    }
    for (_, _, span) in &mut entries.stereo_atoms {
        normalize_entity_span(span);
    }
    for (_, _, span) in &mut entries.stereo_bonds {
        normalize_entity_span(span);
    }
}

/// Map every union id into one projected id space. Present entities occupy the valid dense prefix;
/// absent entities follow it, so a surviving reference to one is retained but rejected by the
/// molecule-entry validator.
fn projected_ids<Id>(presence: impl IntoIterator<Item = bool>) -> HashMap<Id, Id>
where
    Id: Copy + Eq + Hash + From<usize>,
{
    let presence: Vec<bool> = presence.into_iter().collect();
    let mut next_present = 0;
    let mut next_absent = presence.iter().filter(|&&present| present).count();
    presence
        .into_iter()
        .enumerate()
        .map(|(index, present)| {
            let projected = if present {
                let id = next_present;
                next_present += 1;
                id
            } else {
                let id = next_absent;
                next_absent += 1;
                id
            };
            (Id::from(index), Id::from(projected))
        })
        .collect()
}

#[allow(clippy::type_complexity)]
impl ReactionSpanAst {
    /// Superimpose two molecules over their correspondence into the reaction span. Matched entities
    /// become `Unchanged` / `Modified` carrying both molecules' actual values; entities unmatched
    /// on the lhs become `Removed`, those unmatched on the rhs `Added`. Lhs-anchored: lhs ids kept,
    /// right-unmatched entities appended, rhs participants and constraints remapped into that union
    /// frame. Returns `None` when any correspondence family declares counts different from the
    /// supplied molecules or a matched bond, overlay, or stereo entity has incompatible incidence
    /// under the atom correspondence. A compatible correspondence may leave otherwise matchable
    /// entities unmatched; they are represented as removals and additions.
    ///
    /// # Semantic properties
    ///
    /// For a compatible correspondence, the lhs projection is structurally equal to `lhs`. The rhs
    /// projection is `rhs` reindexed into the lhs-anchored reaction frame: preserved entities retain
    /// lhs ids and rhs-only entities are appended. It is equivalent to `rhs` under the induced total
    /// correspondence, but need not be structurally equal when matched pairs cross entity order.
    /// The correspondence recovered from the span relates these normalized projections and
    /// therefore need not equal the source correspondence.
    pub fn superimpose(
        lhs: &MoleculeAst,
        rhs: &MoleculeAst,
        correspondence: &MoleculeCorrespondence,
    ) -> Option<ReactionSpanAst> {
        if !correspondence_is_compatible(lhs, rhs, correspondence) {
            return None;
        }

        let atoms_corr = correspondence.atoms();
        let bonds_corr = correspondence.bonds();
        let lhs_atom_count = lhs.atoms().count();
        let lhs_bond_count = lhs.bonds().count();

        // R id → union id per family.
        let atom_union: HashMap<AtomId, AtomId> = union_map(atoms_corr, lhs_atom_count);
        let bond_union: HashMap<BondId, BondId> = union_map(bonds_corr, lhs_bond_count);
        let participant_remapping = Remapping::new(
            (0..rhs.atoms().count())
                .map(|index| NodeId::from(atom_union[&AtomId::from(index)]))
                .collect(),
            (0..rhs.bonds().count())
                .map(|index| EdgeId::from(bond_union[&BondId::from(index)]))
                .collect(),
        );

        let remapped_rhs_dative: FixedVarBirelationSet<
            NodeId,
            Ordered,
            1,
            NodeId,
            Unordered,
            DativeBondForm,
        > = FixedVarBirelationSet::new(
            rhs.dative_bonds()
                .iter()
                .map(|view| {
                    (
                        [NodeId::from(view.acceptor_id())],
                        view.donor_ids().map(NodeId::from).collect(),
                        view.ast.clone(),
                    )
                })
                .collect(),
        )
        .apply_remapping(&participant_remapping);
        let remapped_rhs_aromatic: VarRelationSet<NodeId, Unordered, AromaticSystemForm> =
            VarRelationSet::new(
                rhs.aromatic_systems()
                    .iter()
                    .map(|view| {
                        (
                            view.atom_ids().map(NodeId::from).collect(),
                            view.ast.clone(),
                        )
                    })
                    .collect(),
            )
            .apply_remapping(&participant_remapping);
        let remapped_rhs_multicenter: VarRelationSet<NodeId, Unordered, MulticenterBondForm> =
            VarRelationSet::new(
                rhs.multicenter_bonds()
                    .iter()
                    .map(|view| {
                        (
                            view.atom_ids().map(NodeId::from).collect(),
                            view.ast.clone(),
                        )
                    })
                    .collect(),
            )
            .apply_remapping(&participant_remapping);
        let remapped_rhs_noncovalent: FixedRelationSet<NodeId, Unordered, NoncovalentBondForm, 2> =
            FixedRelationSet::new(
                rhs.noncovalent_bonds()
                    .iter()
                    .map(|view| {
                        let [first, second] = view.atom_ids();
                        (
                            [NodeId::from(first), NodeId::from(second)],
                            view.ast.clone(),
                        )
                    })
                    .collect(),
            )
            .apply_remapping(&participant_remapping);
        let remapped_rhs_stereo_atoms: FixedVarBirelationSet<
            NodeId,
            Ordered,
            1,
            StereoLigand,
            Ordered,
            StereoAtomForm,
        > = FixedVarBirelationSet::new(
            rhs.stereo_atoms()
                .iter()
                .map(|view| {
                    (
                        [NodeId::from(view.site_id())],
                        view.ligand_frame(),
                        view.ast.clone(),
                    )
                })
                .collect(),
        )
        .apply_remapping(&participant_remapping);
        let remapped_rhs_stereo_bonds: FixedVarBirelationSet<
            EdgeId,
            Ordered,
            1,
            StereoLigand,
            Ordered,
            StereoBondForm,
        > = FixedVarBirelationSet::new(
            rhs.stereo_bonds()
                .iter()
                .map(|view| {
                    (
                        [EdgeId::from(view.site_id())],
                        view.ligand_frame(),
                        view.ast.clone(),
                    )
                })
                .collect(),
        )
        .apply_remapping(&participant_remapping);

        // Atoms
        let mut atoms: Vec<EntitySpan<AtomForm>> = Vec::new();
        for i in 0..lhs_atom_count {
            let lhs_ast = lhs.atom(AtomId(i as u32)).ast.clone();
            let rhs_ast = atoms_corr
                .right_of(AtomId(i as u32))
                .map(|r| rhs.atom(r).ast.clone());
            atoms.push(EntitySpan::superimpose(Some(lhs_ast), rhs_ast).unwrap());
        }
        for &r in &atoms_corr.right_unmatched() {
            atoms.push(EntitySpan::Added(rhs.atom(r).ast.clone()));
        }

        // Bonds
        let mut bonds: Vec<(AtomId, AtomId, EntitySpan<BondForm>)> = Vec::new();
        for i in 0..lhs_bond_count {
            let [a, b] = lhs.raw_graph().edge_endpoints(EdgeId(i as u32));
            let lhs_ast = lhs.bond(BondId(i as u32)).ast.clone();
            let rhs_ast = bonds_corr
                .right_of(BondId(i as u32))
                .map(|r| rhs.bond(r).ast.clone());
            bonds.push((
                AtomId::from(a),
                AtomId::from(b),
                EntitySpan::superimpose(Some(lhs_ast), rhs_ast).unwrap(),
            ));
        }
        for &r in &bonds_corr.right_unmatched() {
            let [a, b] = rhs.raw_graph().edge_endpoints(EdgeId(r.index() as u32));
            bonds.push((
                atom_union[&AtomId::from(a)],
                atom_union[&AtomId::from(b)],
                EntitySpan::Added(rhs.bond(r).ast.clone()),
            ));
        }

        // Aromatic systems
        let aromatic_corr = correspondence.aromatic_systems();
        let mut aromatic: Vec<(Vec<AtomId>, EntitySpan<AromaticSystemForm>)> = Vec::new();
        for view in lhs.aromatic_systems().iter() {
            let participants: Vec<AtomId> = view.atom_ids().collect();
            let rhs_ast = aromatic_corr
                .right_of(view.id)
                .map(|id| remapped_rhs_aromatic.data(id.into()).clone());
            aromatic.push((
                participants,
                EntitySpan::superimpose(Some(view.ast.clone()), rhs_ast).unwrap(),
            ));
        }
        for &r in &aromatic_corr.right_unmatched() {
            let relation_id = RelationId::from(r);
            let participants = remapped_rhs_aromatic
                .participants(relation_id)
                .iter()
                .copied()
                .map(AtomId::from)
                .collect();
            aromatic.push((
                participants,
                EntitySpan::Added(remapped_rhs_aromatic.data(relation_id).clone()),
            ));
        }

        // Multicenter bonds
        let multicenter_corr = correspondence.multicenter_bonds();
        let mut multicenter: Vec<(Vec<AtomId>, EntitySpan<MulticenterBondForm>)> = Vec::new();
        for view in lhs.multicenter_bonds().iter() {
            let participants: Vec<AtomId> = view.atom_ids().collect();
            let rhs_ast = multicenter_corr
                .right_of(view.id)
                .map(|id| remapped_rhs_multicenter.data(id.into()).clone());
            multicenter.push((
                participants,
                EntitySpan::superimpose(Some(view.ast.clone()), rhs_ast).unwrap(),
            ));
        }
        for &r in &multicenter_corr.right_unmatched() {
            let relation_id = RelationId::from(r);
            let participants = remapped_rhs_multicenter
                .participants(relation_id)
                .iter()
                .copied()
                .map(AtomId::from)
                .collect();
            multicenter.push((
                participants,
                EntitySpan::Added(remapped_rhs_multicenter.data(relation_id).clone()),
            ));
        }

        // Noncovalent bonds
        let noncovalent_corr = correspondence.noncovalent_bonds();
        let mut noncovalent: Vec<(AtomId, AtomId, EntitySpan<NoncovalentBondForm>)> = Vec::new();
        for view in lhs.noncovalent_bonds().iter() {
            let [a, b] = view.atom_ids();
            let rhs_ast = noncovalent_corr
                .right_of(view.id)
                .map(|id| remapped_rhs_noncovalent.data(id.into()).clone());
            noncovalent.push((
                a,
                b,
                EntitySpan::superimpose(Some(view.ast.clone()), rhs_ast).unwrap(),
            ));
        }
        for &r in &noncovalent_corr.right_unmatched() {
            let relation_id = RelationId::from(r);
            let &[first, second] = remapped_rhs_noncovalent.participants(relation_id);
            noncovalent.push((
                AtomId::from(first),
                AtomId::from(second),
                EntitySpan::Added(remapped_rhs_noncovalent.data(relation_id).clone()),
            ));
        }

        // Dative bonds
        let dative_corr = correspondence.dative_bonds();
        let mut dative: Vec<(Vec<AtomId>, AtomId, EntitySpan<DativeBondForm>)> = Vec::new();
        for view in lhs.dative_bonds().iter() {
            let acceptor = view.acceptor_id();
            let donors = view.donor_ids().collect();
            let rhs_ast = dative_corr
                .right_of(view.id)
                .map(|id| remapped_rhs_dative.data(id.into()).clone());
            dative.push((
                donors,
                acceptor,
                EntitySpan::superimpose(Some(view.ast.clone()), rhs_ast).unwrap(),
            ));
        }
        for &r in &dative_corr.right_unmatched() {
            let relation_id = RelationId::from(r);
            let acceptor = AtomId::from(remapped_rhs_dative.participants_1(relation_id)[0]);
            let donors = remapped_rhs_dative
                .participants_2(relation_id)
                .iter()
                .copied()
                .map(AtomId::from)
                .collect();
            dative.push((
                donors,
                acceptor,
                EntitySpan::Added(remapped_rhs_dative.data(relation_id).clone()),
            ));
        }

        // Stereo atoms
        let stereo_atom_corr = correspondence.stereo_atoms();
        let mut stereo_atoms: Vec<(AtomId, Vec<StereoLigand>, EntitySpan<StereoAtomForm>)> =
            Vec::new();
        for view in lhs.stereo_atoms().iter() {
            let rhs_ast = stereo_atom_corr
                .right_of(view.id)
                .map(|id| remapped_rhs_stereo_atoms.data(id.into()).clone());
            stereo_atoms.push((
                view.site_id(),
                view.ligand_frame(),
                EntitySpan::superimpose(Some(view.ast.clone()), rhs_ast).unwrap(),
            ));
        }
        for &r in &stereo_atom_corr.right_unmatched() {
            let relation_id = RelationId::from(r);
            let site = AtomId::from(remapped_rhs_stereo_atoms.participants_1(relation_id)[0]);
            let ligands = remapped_rhs_stereo_atoms
                .participants_2(relation_id)
                .to_vec();
            stereo_atoms.push((
                site,
                ligands,
                EntitySpan::Added(remapped_rhs_stereo_atoms.data(relation_id).clone()),
            ));
        }

        // Stereo bonds
        let stereo_bond_corr = correspondence.stereo_bonds();
        let mut stereo_bonds: Vec<(BondId, Vec<StereoLigand>, EntitySpan<StereoBondForm>)> =
            Vec::new();
        for view in lhs.stereo_bonds().iter() {
            let rhs_ast = stereo_bond_corr
                .right_of(view.id)
                .map(|id| remapped_rhs_stereo_bonds.data(id.into()).clone());
            stereo_bonds.push((
                view.site_id(),
                view.ligand_frame(),
                EntitySpan::superimpose(Some(view.ast.clone()), rhs_ast).unwrap(),
            ));
        }
        for &r in &stereo_bond_corr.right_unmatched() {
            let relation_id = RelationId::from(r);
            let site = BondId::from(remapped_rhs_stereo_bonds.participants_1(relation_id)[0]);
            let ligands = remapped_rhs_stereo_bonds
                .participants_2(relation_id)
                .to_vec();
            stereo_bonds.push((
                site,
                ligands,
                EntitySpan::Added(remapped_rhs_stereo_bonds.data(relation_id).clone()),
            ));
        }

        // Constraints: R's remapped into the union frame, then set-diffed against L's.
        let remapping = IdRemapping::new(
            atom_union.iter().map(|(&r, &u)| (r, u)).collect(),
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

        Some(ReactionSpanAst::from_entries(ReactionSpanEntries {
            atoms,
            bonds,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            stereo_atoms,
            stereo_bonds,
            constraints,
        }))
    }

    /// Recover the per-family correspondence between the two normalized side projections,
    /// forgetting the values. For every span,
    /// `superimpose(self.lhs(), self.rhs(), &self.correspondence())` reproduces `Some(self)`.
    /// A source correspondence used to construct the span is not retained when it assigns a
    /// different rhs entity order.
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

    pub fn atoms(&self) -> &[EntitySpan<AtomForm>] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[EntitySpan<BondForm>] {
        &self.bonds
    }

    pub(crate) fn dative_bonds(
        &self,
    ) -> &FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, EntitySpan<DativeBondForm>>
    {
        &self.dative_bonds
    }

    pub(crate) fn aromatic_systems(
        &self,
    ) -> &VarRelationSet<NodeId, Unordered, EntitySpan<AromaticSystemForm>> {
        &self.aromatic_systems
    }

    pub(crate) fn multicenter_bonds(
        &self,
    ) -> &VarRelationSet<NodeId, Unordered, EntitySpan<MulticenterBondForm>> {
        &self.multicenter_bonds
    }

    pub(crate) fn noncovalent_bonds(
        &self,
    ) -> &FixedRelationSet<NodeId, Unordered, EntitySpan<NoncovalentBondForm>, 2> {
        &self.noncovalent_bonds
    }

    pub(crate) fn stereo_atoms(
        &self,
    ) -> &FixedVarBirelationSet<NodeId, Ordered, 1, StereoLigand, Ordered, EntitySpan<StereoAtomForm>>
    {
        &self.stereo_atoms
    }

    pub(crate) fn stereo_bonds(
        &self,
    ) -> &FixedVarBirelationSet<EdgeId, Ordered, 1, StereoLigand, Ordered, EntitySpan<StereoBondForm>>
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
        let mut deltas = Deltas::new();
        AtomDelta::append_deltas_from_states(&self.atoms, |_| (), &mut deltas);
        BondDelta::append_deltas_from_states(
            &self.bonds,
            |edge| {
                let [a, b] = self.graph.edge_endpoints(EdgeId(edge as u32));
                [AtomId::from(a), AtomId::from(b)]
            },
            &mut deltas,
        );
        let dative_states: Vec<EntitySpan<DativeBondForm>> = (0..self.dative_bonds.count())
            .map(|i| self.dative_bonds.data(RelationId(i as u32)).clone())
            .collect();
        DativeBondDelta::append_deltas_from_states(
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
            &mut deltas,
        );
        let aromatic_states: Vec<EntitySpan<AromaticSystemForm>> =
            (0..self.aromatic_systems.count())
                .map(|i| self.aromatic_systems.data(RelationId(i as u32)).clone())
                .collect();
        AromaticSystemDelta::append_deltas_from_states(
            &aromatic_states,
            |index| {
                self.aromatic_systems
                    .participants(RelationId(index as u32))
                    .iter()
                    .map(|&n| AtomId::from(n))
                    .collect()
            },
            &mut deltas,
        );
        let multicenter_states: Vec<EntitySpan<MulticenterBondForm>> =
            (0..self.multicenter_bonds.count())
                .map(|i| self.multicenter_bonds.data(RelationId(i as u32)).clone())
                .collect();
        MulticenterBondDelta::append_deltas_from_states(
            &multicenter_states,
            |index| {
                self.multicenter_bonds
                    .participants(RelationId(index as u32))
                    .iter()
                    .map(|&n| AtomId::from(n))
                    .collect()
            },
            &mut deltas,
        );
        let noncovalent_states: Vec<EntitySpan<NoncovalentBondForm>> =
            (0..self.noncovalent_bonds.count())
                .map(|i| self.noncovalent_bonds.data(RelationId(i as u32)).clone())
                .collect();
        NoncovalentBondDelta::append_deltas_from_states(
            &noncovalent_states,
            |index| {
                let [a, b] = *self
                    .noncovalent_bonds
                    .participants(RelationId(index as u32));
                [AtomId::from(a), AtomId::from(b)]
            },
            &mut deltas,
        );
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
                } => {
                    for delta in StereoAtomDelta::diff(id, left, right) {
                        deltas.push(Delta::StereoAtom(delta));
                    }
                }
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
                } => {
                    for delta in StereoBondDelta::diff(id, left, right) {
                        deltas.push(Delta::StereoBond(delta));
                    }
                }
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
        ReactionAst::new(self.lhs(), deltas)
    }

    /// Project one `Side` into flat molecule entries. Every selected entity is retained. References
    /// to an entity absent from the side map beyond the valid dense prefix so the molecule-entry
    /// validator reports them rather than silently dropping the referring entry.
    fn project_entries(&self, side: Side) -> MoleculeEntries {
        let atom_ids: HashMap<AtomId, AtomId> = projected_ids(
            self.atoms
                .iter()
                .map(|span| entity_side(span, side).is_some()),
        );
        let bond_ids: HashMap<BondId, BondId> = projected_ids(
            self.bonds
                .iter()
                .map(|span| entity_side(span, side).is_some()),
        );
        let dative_ids: HashMap<DativeBondId, DativeBondId> =
            projected_ids((0..self.dative_bonds.count()).map(|i| {
                entity_side(self.dative_bonds.data(RelationId(i as u32)), side).is_some()
            }));
        let aromatic_ids: HashMap<AromaticSystemId, AromaticSystemId> =
            projected_ids((0..self.aromatic_systems.count()).map(|i| {
                entity_side(self.aromatic_systems.data(RelationId(i as u32)), side).is_some()
            }));
        let multicenter_ids: HashMap<MulticenterBondId, MulticenterBondId> =
            projected_ids((0..self.multicenter_bonds.count()).map(|i| {
                entity_side(self.multicenter_bonds.data(RelationId(i as u32)), side).is_some()
            }));
        let noncovalent_ids: HashMap<NoncovalentBondId, NoncovalentBondId> =
            projected_ids((0..self.noncovalent_bonds.count()).map(|i| {
                entity_side(self.noncovalent_bonds.data(RelationId(i as u32)), side).is_some()
            }));
        let stereo_atom_ids: HashMap<StereoAtomId, StereoAtomId> =
            projected_ids((0..self.stereo_atoms.count()).map(|i| {
                entity_side(self.stereo_atoms.data(RelationId(i as u32)), side).is_some()
            }));
        let stereo_bond_ids: HashMap<StereoBondId, StereoBondId> =
            projected_ids((0..self.stereo_bonds.count()).map(|i| {
                entity_side(self.stereo_bonds.data(RelationId(i as u32)), side).is_some()
            }));

        let atoms = self
            .atoms
            .iter()
            .filter_map(|span| entity_side(span, side).cloned())
            .collect();
        let bonds = self
            .bonds
            .iter()
            .enumerate()
            .filter_map(|(index, span)| {
                let ast = entity_side(span, side)?;
                let [first, second] = self.graph.edge_endpoints(EdgeId(index as u32));
                Some((
                    atom_ids[&AtomId::from(first)],
                    atom_ids[&AtomId::from(second)],
                    ast.clone(),
                ))
            })
            .collect();
        let dative = (0..self.dative_bonds.count())
            .filter_map(|index| {
                let id = RelationId(index as u32);
                let ast = entity_side(self.dative_bonds.data(id), side)?;
                let donors = self
                    .dative_bonds
                    .participants_2(id)
                    .iter()
                    .map(|&atom| atom_ids[&AtomId::from(atom)])
                    .collect();
                let acceptor = atom_ids[&AtomId::from(self.dative_bonds.participants_1(id)[0])];
                Some((donors, acceptor, ast.clone()))
            })
            .collect();
        let aromatic = (0..self.aromatic_systems.count())
            .filter_map(|index| {
                let id = RelationId(index as u32);
                let ast = entity_side(self.aromatic_systems.data(id), side)?;
                let atoms = self
                    .aromatic_systems
                    .participants(id)
                    .iter()
                    .map(|&atom| atom_ids[&AtomId::from(atom)])
                    .collect();
                Some((atoms, ast.clone()))
            })
            .collect();
        let multicenter = (0..self.multicenter_bonds.count())
            .filter_map(|index| {
                let id = RelationId(index as u32);
                let ast = entity_side(self.multicenter_bonds.data(id), side)?;
                let atoms = self
                    .multicenter_bonds
                    .participants(id)
                    .iter()
                    .map(|&atom| atom_ids[&AtomId::from(atom)])
                    .collect();
                Some((atoms, ast.clone()))
            })
            .collect();
        let noncovalent = (0..self.noncovalent_bonds.count())
            .filter_map(|index| {
                let id = RelationId(index as u32);
                let ast = entity_side(self.noncovalent_bonds.data(id), side)?;
                let [first, second] = *self.noncovalent_bonds.participants(id);
                Some((
                    atom_ids[&AtomId::from(first)],
                    atom_ids[&AtomId::from(second)],
                    ast.clone(),
                ))
            })
            .collect();
        let remap_ligands = |ligands: &[StereoLigand]| {
            ligands
                .iter()
                .map(|ligand| StereoLigand::new(atom_ids[&ligand.atom_id], ligand.kind))
                .collect()
        };
        let stereo_atoms = (0..self.stereo_atoms.count())
            .filter_map(|index| {
                let id = RelationId(index as u32);
                let ast = entity_side(self.stereo_atoms.data(id), side)?;
                let site = atom_ids[&AtomId::from(self.stereo_atoms.participants_1(id)[0])];
                Some((
                    site,
                    remap_ligands(self.stereo_atoms.participants_2(id)),
                    ast.clone(),
                ))
            })
            .collect();
        let stereo_bonds = (0..self.stereo_bonds.count())
            .filter_map(|index| {
                let id = RelationId(index as u32);
                let ast = entity_side(self.stereo_bonds.data(id), side)?;
                let site = bond_ids[&BondId::from(self.stereo_bonds.participants_1(id)[0])];
                Some((
                    site,
                    remap_ligands(self.stereo_bonds.participants_2(id)),
                    ast.clone(),
                ))
            })
            .collect();

        let remapping = IdRemapping::new(
            atom_ids,
            bond_ids,
            dative_ids,
            aromatic_ids,
            multicenter_ids,
            noncovalent_ids,
            stereo_atom_ids,
            stereo_bond_ids,
        );
        let constraints = self
            .constraints
            .iter()
            .filter_map(|span| match side {
                Side::Left => span.lhs(),
                Side::Right => span.rhs(),
            })
            .cloned()
            .map(|constraint| constraint.remap(&remapping))
            .collect();

        MoleculeEntries {
            atoms,
            bonds,
            dative,
            aromatic,
            multicenter,
            noncovalent,
            stereo_atoms,
            stereo_bonds,
            constraints,
        }
    }

    /// Project one side to the molecule established by the reaction-span construction invariant.
    fn project(&self, side: Side) -> MoleculeAst {
        MoleculeAst::from_entries(self.project_entries(side))
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

        let mut removed_atoms: HashMap<AtomId, AtomForm> = HashMap::new();
        let mut added_atoms: BTreeMap<AtomId, AtomForm> = BTreeMap::new();
        let mut atom_changes: HashMap<AtomId, Vec<AtomDelta>> = HashMap::new();
        let mut removed_bonds: HashMap<BondId, BondForm> = HashMap::new();
        let mut added_bonds: BTreeMap<BondId, ([AtomId; 2], BondForm)> = BTreeMap::new();
        let mut bond_changes: HashMap<BondId, Vec<BondDelta>> = HashMap::new();
        let mut removed_aromatic: HashMap<AromaticSystemId, AromaticSystemForm> = HashMap::new();
        let mut added_aromatic: BTreeMap<AromaticSystemId, (Vec<AtomId>, AromaticSystemForm)> =
            BTreeMap::new();
        let mut aromatic_changes: HashMap<AromaticSystemId, Vec<AromaticSystemDelta>> =
            HashMap::new();
        let mut removed_dative: HashMap<DativeBondId, DativeBondForm> = HashMap::new();
        let mut added_dative: BTreeMap<DativeBondId, (Vec<AtomId>, AtomId, DativeBondForm)> =
            BTreeMap::new();
        let mut dative_changes: HashMap<DativeBondId, Vec<DativeBondDelta>> = HashMap::new();
        let mut removed_multicenter: HashMap<MulticenterBondId, MulticenterBondForm> =
            HashMap::new();
        let mut added_multicenter: BTreeMap<MulticenterBondId, (Vec<AtomId>, MulticenterBondForm)> =
            BTreeMap::new();
        let mut multicenter_changes: HashMap<MulticenterBondId, Vec<MulticenterBondDelta>> =
            HashMap::new();
        let mut removed_noncovalent: HashMap<NoncovalentBondId, NoncovalentBondForm> =
            HashMap::new();
        let mut added_noncovalent: BTreeMap<NoncovalentBondId, ([AtomId; 2], NoncovalentBondForm)> =
            BTreeMap::new();
        let mut noncovalent_changes: HashMap<NoncovalentBondId, Vec<NoncovalentBondDelta>> =
            HashMap::new();
        let mut removed_stereo_atom: HashMap<StereoAtomId, StereoAtomForm> = HashMap::new();
        let mut added_stereo_atom: BTreeMap<
            StereoAtomId,
            (AtomId, Vec<StereoLigand>, StereoAtomForm),
        > = BTreeMap::new();
        let mut stereo_atom_changes: HashMap<StereoAtomId, Vec<StereoAtomDelta>> = HashMap::new();
        let mut removed_stereo_bond: HashMap<StereoBondId, StereoBondForm> = HashMap::new();
        let mut added_stereo_bond: BTreeMap<
            StereoBondId,
            (BondId, Vec<StereoLigand>, StereoBondForm),
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

        let mut atoms: Vec<EntitySpan<AtomForm>> =
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

        let mut bonds: Vec<(AtomId, AtomId, EntitySpan<BondForm>)> =
            Vec::with_capacity(bond_count + added_bonds.len());
        for edge in 0..bond_count {
            let id = BondId(edge as u32);
            let [a, b] = lhs.raw_graph().edge_endpoints(EdgeId(edge as u32));
            let first = AtomId::from(a);
            let second = AtomId::from(b);
            if let Some(ast) = removed_bonds.get(&id) {
                bonds.push((first, second, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = bond_changes.get(&id) {
                let left = lhs.bond(id).ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_bond_change(&mut right, change)?;
                }
                bonds.push((
                    first,
                    second,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                bonds.push((
                    first,
                    second,
                    EntitySpan::Unchanged(lhs.bond(id).ast.clone()),
                ));
            }
        }
        for (atoms, ast) in added_bonds.into_values() {
            bonds.push((
                AtomId(atom_index[&atoms[0]] as u32),
                AtomId(atom_index[&atoms[1]] as u32),
                EntitySpan::Added(ast),
            ));
        }

        // Overlay columns: lhs overlays tagged by their fold (Removed/Modified/Unchanged),
        // created overlays appended; participants mapped to the union id space via `atom_index`.
        let mut aromatic: Vec<(Vec<AtomId>, EntitySpan<AromaticSystemForm>)> = Vec::new();
        for view in lhs.aromatic_systems().iter() {
            let participants: Vec<AtomId> = view
                .atom_ids()
                .map(|a| AtomId(atom_index[&a] as u32))
                .collect();
            if let Some(ast) = removed_aromatic.get(&view.id) {
                aromatic.push((participants, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = aromatic_changes.get(&view.id) {
                let left = view.ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_aromatic_change(&mut right, change)?;
                }
                aromatic.push((
                    participants,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                aromatic.push((participants, EntitySpan::Unchanged(view.ast.clone())));
            }
        }
        for (atoms, ast) in added_aromatic.into_values() {
            let participants = atoms.iter().map(|a| AtomId(atom_index[a] as u32)).collect();
            aromatic.push((participants, EntitySpan::Added(ast)));
        }

        let mut multicenter: Vec<(Vec<AtomId>, EntitySpan<MulticenterBondForm>)> = Vec::new();
        for view in lhs.multicenter_bonds().iter() {
            let participants: Vec<AtomId> = view
                .atom_ids()
                .map(|a| AtomId(atom_index[&a] as u32))
                .collect();
            if let Some(ast) = removed_multicenter.get(&view.id) {
                multicenter.push((participants, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = multicenter_changes.get(&view.id) {
                let left = view.ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_multicenter_change(&mut right, change)?;
                }
                multicenter.push((
                    participants,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                multicenter.push((participants, EntitySpan::Unchanged(view.ast.clone())));
            }
        }
        for (atoms, ast) in added_multicenter.into_values() {
            let participants = atoms.iter().map(|a| AtomId(atom_index[a] as u32)).collect();
            multicenter.push((participants, EntitySpan::Added(ast)));
        }

        let mut noncovalent: Vec<(AtomId, AtomId, EntitySpan<NoncovalentBondForm>)> = Vec::new();
        for view in lhs.noncovalent_bonds().iter() {
            let [a, b] = view.atom_ids();
            let first = AtomId(atom_index[&a] as u32);
            let second = AtomId(atom_index[&b] as u32);
            if let Some(ast) = removed_noncovalent.get(&view.id) {
                noncovalent.push((first, second, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = noncovalent_changes.get(&view.id) {
                let left = view.ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_noncovalent_change(&mut right, change)?;
                }
                noncovalent.push((
                    first,
                    second,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                noncovalent.push((first, second, EntitySpan::Unchanged(view.ast.clone())));
            }
        }
        for ([a, b], ast) in added_noncovalent.into_values() {
            noncovalent.push((
                AtomId(atom_index[&a] as u32),
                AtomId(atom_index[&b] as u32),
                EntitySpan::Added(ast),
            ));
        }

        let mut dative: Vec<(Vec<AtomId>, AtomId, EntitySpan<DativeBondForm>)> = Vec::new();
        for view in lhs.dative_bonds().iter() {
            let acceptor = AtomId(atom_index[&view.acceptor_id()] as u32);
            let donors: Vec<AtomId> = view
                .donor_ids()
                .map(|a| AtomId(atom_index[&a] as u32))
                .collect();
            if let Some(ast) = removed_dative.get(&view.id) {
                dative.push((donors, acceptor, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = dative_changes.get(&view.id) {
                let left = view.ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_dative_change(&mut right, change)?;
                }
                dative.push((
                    donors,
                    acceptor,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                dative.push((donors, acceptor, EntitySpan::Unchanged(view.ast.clone())));
            }
        }
        for (donors, acceptor, ast) in added_dative.into_values() {
            let acceptor = AtomId(atom_index[&acceptor] as u32);
            let donors = donors
                .iter()
                .map(|a| AtomId(atom_index[a] as u32))
                .collect();
            dative.push((donors, acceptor, EntitySpan::Added(ast)));
        }

        // Stereo overlays: lhs entities tagged by their fold (Removed/Modified/Unchanged), created
        // ones appended. Site/ligand ids mapped to the union frame via `atom_index`/`bond_index`.
        let mut stereo_atoms: Vec<(AtomId, Vec<StereoLigand>, EntitySpan<StereoAtomForm>)> =
            Vec::new();
        for view in lhs.stereo_atoms().iter() {
            let site = view.site_id();
            let ligands = view.ligand_frame();
            if let Some(ast) = removed_stereo_atom.get(&view.id) {
                stereo_atoms.push((site, ligands, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = stereo_atom_changes.get(&view.id) {
                let left = view.ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_stereo_atom_change(&mut right, change)?;
                }
                stereo_atoms.push((
                    site,
                    ligands,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                stereo_atoms.push((site, ligands, EntitySpan::Unchanged(view.ast.clone())));
            }
        }
        for (site, ligands, ast) in added_stereo_atom.into_values() {
            let site = AtomId(atom_index[&site] as u32);
            let ligands = ligands
                .iter()
                .map(|l| StereoLigand::new(AtomId(atom_index[&l.atom_id] as u32), l.kind))
                .collect();
            stereo_atoms.push((site, ligands, EntitySpan::Added(ast)));
        }

        let mut stereo_bonds: Vec<(BondId, Vec<StereoLigand>, EntitySpan<StereoBondForm>)> =
            Vec::new();
        for view in lhs.stereo_bonds().iter() {
            let site = view.site_id();
            let ligands = view.ligand_frame();
            if let Some(ast) = removed_stereo_bond.get(&view.id) {
                stereo_bonds.push((site, ligands, EntitySpan::Removed(ast.clone())));
            } else if let Some(changes) = stereo_bond_changes.get(&view.id) {
                let left = view.ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_stereo_bond_change(&mut right, change)?;
                }
                stereo_bonds.push((
                    site,
                    ligands,
                    EntitySpan::Modified {
                        lhs: left,
                        rhs: right,
                    },
                ));
            } else {
                stereo_bonds.push((site, ligands, EntitySpan::Unchanged(view.ast.clone())));
            }
        }
        for (site, ligands, ast) in added_stereo_bond.into_values() {
            let site = BondId(bond_index[&site] as u32);
            let ligands = ligands
                .iter()
                .map(|l| StereoLigand::new(AtomId(atom_index[&l.atom_id] as u32), l.kind))
                .collect();
            stereo_bonds.push((site, ligands, EntitySpan::Added(ast)));
        }

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

        ReactionSpanAst::try_from_entries(ReactionSpanEntries {
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
        .map_err(|_| Contradiction)
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

    use super::super::constraint::{
        AromaticSystemConstraintAst, AtomConstraintAst, BondConstraintAst, Constraint, Constraints,
        DativeBondConstraintAst, MoleculeConstraint, MulticenterBondConstraintAst,
        NoncovalentBondConstraintAst, StereoAtomConstraintAst, StereoBondConstraintAst,
        StereogenicityAst,
    };
    use super::super::delta::Deltas;
    use super::super::edit::{BondFieldChange, NoncovalentBondFieldChange, StereoAtomFieldChange};
    use super::super::ligand::{StereoLigand, StereoLigandKind};
    use super::super::noncovalent::{NoncovalentBondKind, NoncovalentBondKindForm};
    use super::super::stereo::{StereoConfigurationForm, StereoCoset, StereoKind};
    use super::super::value::NumForm;
    use super::*;

    #[rstest]
    fn test_reaction_span_ast_from_entries() {
        let span = ReactionSpanAst::from_entries(ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                EntitySpan::Modified {
                    lhs: AtomForm::from_element(Element::C),
                    rhs: AtomForm::from_element(Element::N),
                },
                EntitySpan::Added(AtomForm::from_element(Element::O)),
                EntitySpan::Removed(AtomForm::from_element(Element::F)),
                EntitySpan::Unchanged(AtomForm::from_element(Element::Cl)),
            ],
            bonds: vec![
                (
                    AtomId(0),
                    AtomId(4),
                    EntitySpan::Unchanged(BondForm::from_order(1)),
                ),
                (
                    AtomId(0),
                    AtomId(1),
                    EntitySpan::Modified {
                        lhs: BondForm::from_order(1),
                        rhs: BondForm::from_order(2),
                    },
                ),
                (
                    AtomId(1),
                    AtomId(2),
                    EntitySpan::Added(BondForm::from_order(1)),
                ),
                (
                    AtomId(1),
                    AtomId(3),
                    EntitySpan::Removed(BondForm::from_order(1)),
                ),
            ],
            dative: vec![(
                vec![AtomId(1)],
                AtomId(0),
                EntitySpan::Unchanged(DativeBondForm::default()),
            )],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(4)],
                EntitySpan::Unchanged(AromaticSystemForm::default()),
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(4)],
                EntitySpan::Unchanged(MulticenterBondForm::default()),
            )],
            noncovalent: vec![(
                AtomId(0),
                AtomId(4),
                EntitySpan::Unchanged(NoncovalentBondForm::default()),
            )],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                EntitySpan::Unchanged(StereoAtomForm::default()),
            )],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                EntitySpan::Unchanged(StereoBondForm::default()),
            )],
            constraints: vec![
                ConstraintSpan::Unchanged(Constraint::Molecule(MoleculeConstraint::Connected {
                    atoms: Some(vec![AtomId(0), AtomId(1)]),
                })),
                ConstraintSpan::Added(Constraint::Molecule(MoleculeConstraint::Connected {
                    atoms: Some(vec![AtomId(1), AtomId(2)]),
                })),
                ConstraintSpan::Removed(Constraint::Molecule(MoleculeConstraint::Connected {
                    atoms: Some(vec![AtomId(1), AtomId(3)]),
                })),
            ],
        });

        assert_eq!(
            span.lhs(),
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::F),
                    AtomForm::from_element(Element::Cl),
                ],
                bonds: vec![
                    (AtomId(0), AtomId(3), BondForm::from_order(1)),
                    (AtomId(0), AtomId(1), BondForm::from_order(1)),
                    (AtomId(1), AtomId(2), BondForm::from_order(1)),
                ],
                dative: vec![(vec![AtomId(1)], AtomId(0), DativeBondForm::default())],
                aromatic: vec![(
                    vec![AtomId(0), AtomId(1), AtomId(3)],
                    AromaticSystemForm::default(),
                )],
                multicenter: vec![(
                    vec![AtomId(0), AtomId(1), AtomId(3)],
                    MulticenterBondForm::default(),
                )],
                noncovalent: vec![(AtomId(0), AtomId(3), NoncovalentBondForm::default(),)],
                stereo_atoms: vec![(
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::default(),
                )],
                stereo_bonds: vec![(
                    BondId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    ],
                    StereoBondForm::default(),
                )],
                constraints: Constraints::from_iter([
                    Constraint::Molecule(MoleculeConstraint::Connected {
                        atoms: Some(vec![AtomId(0), AtomId(1)]),
                    }),
                    Constraint::Molecule(MoleculeConstraint::Connected {
                        atoms: Some(vec![AtomId(1), AtomId(2)]),
                    }),
                ]),
            }),
        );
        assert_eq!(
            span.rhs(),
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::N),
                    AtomForm::from_element(Element::O),
                    AtomForm::from_element(Element::Cl),
                ],
                bonds: vec![
                    (AtomId(0), AtomId(3), BondForm::from_order(1)),
                    (AtomId(0), AtomId(1), BondForm::from_order(2)),
                    (AtomId(1), AtomId(2), BondForm::from_order(1)),
                ],
                dative: vec![(vec![AtomId(1)], AtomId(0), DativeBondForm::default())],
                aromatic: vec![(
                    vec![AtomId(0), AtomId(1), AtomId(3)],
                    AromaticSystemForm::default(),
                )],
                multicenter: vec![(
                    vec![AtomId(0), AtomId(1), AtomId(3)],
                    MulticenterBondForm::default(),
                )],
                noncovalent: vec![(AtomId(0), AtomId(3), NoncovalentBondForm::default(),)],
                stereo_atoms: vec![(
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::default(),
                )],
                stereo_bonds: vec![(
                    BondId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    ],
                    StereoBondForm::default(),
                )],
                constraints: Constraints::from_iter([
                    Constraint::Molecule(MoleculeConstraint::Connected {
                        atoms: Some(vec![AtomId(0), AtomId(1)]),
                    }),
                    Constraint::Molecule(MoleculeConstraint::Connected {
                        atoms: Some(vec![AtomId(1), AtomId(2)]),
                    }),
                ]),
            }),
        );
    }

    #[rstest]
    fn test_reaction_span_ast_from_entries_normalization() {
        let lhs_atom = AtomForm::from_element(Element::C).with_charge(NumForm::Lit(1));
        let rhs_atom = AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([1_i64]));
        assert_ne!(lhs_atom, rhs_atom);
        let span = ReactionSpanAst::from_entries(ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Modified {
                    lhs: lhs_atom.clone(),
                    rhs: rhs_atom,
                },
                EntitySpan::Unchanged(AtomForm::from_element(Element::O)),
            ],
            bonds: vec![(
                AtomId(0),
                AtomId(1),
                EntitySpan::Modified {
                    lhs: BondForm::default(),
                    rhs: BondForm::default(),
                },
            )],
            dative: vec![(
                vec![AtomId(1)],
                AtomId(0),
                EntitySpan::Modified {
                    lhs: DativeBondForm::default(),
                    rhs: DativeBondForm::default(),
                },
            )],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1)],
                EntitySpan::Modified {
                    lhs: AromaticSystemForm::default(),
                    rhs: AromaticSystemForm::default(),
                },
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1)],
                EntitySpan::Modified {
                    lhs: MulticenterBondForm::default(),
                    rhs: MulticenterBondForm::default(),
                },
            )],
            noncovalent: vec![(
                AtomId(0),
                AtomId(1),
                EntitySpan::Modified {
                    lhs: NoncovalentBondForm::default(),
                    rhs: NoncovalentBondForm::default(),
                },
            )],
            stereo_atoms: vec![(
                AtomId(0),
                vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
                EntitySpan::Modified {
                    lhs: StereoAtomForm::default(),
                    rhs: StereoAtomForm::default(),
                },
            )],
            stereo_bonds: vec![(
                BondId(0),
                vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
                EntitySpan::Modified {
                    lhs: StereoBondForm::default(),
                    rhs: StereoBondForm::default(),
                },
            )],
            constraints: Vec::new(),
        });

        assert_eq!(span.atoms()[0], EntitySpan::Unchanged(lhs_atom));
        assert_eq!(span.bonds()[0], EntitySpan::Unchanged(BondForm::default()));
        assert_eq!(
            span.dative_bonds().data(RelationId(0)),
            &EntitySpan::Unchanged(DativeBondForm::default())
        );
        assert_eq!(
            span.aromatic_systems().data(RelationId(0)),
            &EntitySpan::Unchanged(AromaticSystemForm::default())
        );
        assert_eq!(
            span.multicenter_bonds().data(RelationId(0)),
            &EntitySpan::Unchanged(MulticenterBondForm::default())
        );
        assert_eq!(
            span.noncovalent_bonds().data(RelationId(0)),
            &EntitySpan::Unchanged(NoncovalentBondForm::default())
        );
        assert_eq!(
            span.stereo_atoms().data(RelationId(0)),
            &EntitySpan::Unchanged(StereoAtomForm::default())
        );
        assert_eq!(
            span.stereo_bonds().data(RelationId(0)),
            &EntitySpan::Unchanged(StereoBondForm::default())
        );
    }

    #[rstest]
    #[case::bond(
        ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Removed(AtomForm::from_element(Element::C)),
                EntitySpan::Unchanged(AtomForm::from_element(Element::O)),
            ],
            bonds: vec![(
                AtomId(0),
                AtomId(1),
                EntitySpan::Removed(BondForm::from_order(1)),
            )],
            ..Default::default()
        },
        Graph::new(2, &[[0, 1]]),
        vec![
            EntitySpan::Removed(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::O)),
        ],
        vec![EntitySpan::Removed(BondForm::from_order(1))],
        Vec::new(),
    )]
    #[case::constraint(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Removed(AtomForm::from_element(Element::C))],
            constraints: vec![ConstraintSpan::Removed(Constraint::Atom(
                AtomId(0),
                AtomConstraintAst::valence(NumForm::Lit(4)),
            ))],
            ..Default::default()
        },
        Graph::new(1, &[]),
        vec![EntitySpan::Removed(AtomForm::from_element(Element::C))],
        Vec::new(),
        vec![ConstraintSpan::Removed(Constraint::Atom(
            AtomId(0),
            AtomConstraintAst::valence(NumForm::Lit(4)),
        ))],
    )]
    fn test_reaction_span_ast_try_from_entries(
        #[case] entries: ReactionSpanEntries,
        #[case] expected_graph: Graph,
        #[case] expected_atoms: Vec<EntitySpan<AtomForm>>,
        #[case] expected_bonds: Vec<EntitySpan<BondForm>>,
        #[case] expected_constraints: Vec<ConstraintSpan>,
    ) {
        let span = ReactionSpanAst::try_from_entries(entries).unwrap();

        assert_eq!(span.graph(), &expected_graph);
        assert_eq!(span.atoms(), expected_atoms);
        assert_eq!(span.bonds(), expected_bonds);
        assert_eq!(span.constraints(), expected_constraints);
    }

    #[rstest]
    #[should_panic(
        expected = "invalid reaction span entries: reaction span entries reference unavailable atom 1"
    )]
    fn test_reaction_span_ast_from_entries_error() {
        ReactionSpanAst::from_entries(ReactionSpanEntries {
            atoms: vec![EntitySpan::Unchanged(AtomForm::default())],
            bonds: vec![(
                AtomId(0),
                AtomId(1),
                EntitySpan::Unchanged(BondForm::default()),
            )],
            ..Default::default()
        });
    }

    #[rstest]
    #[case::bond_union(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Unchanged(AtomForm::default())],
            bonds: vec![(AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::default()))],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::dative_union(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Unchanged(AtomForm::default())],
            dative: vec![(
                vec![AtomId(1)],
                AtomId(0),
                EntitySpan::Unchanged(DativeBondForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::aromatic_union(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Unchanged(AtomForm::default())],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1)],
                EntitySpan::Unchanged(AromaticSystemForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::multicenter_union(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Unchanged(AtomForm::default())],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1)],
                EntitySpan::Unchanged(MulticenterBondForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::noncovalent_union(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Unchanged(AtomForm::default())],
            noncovalent: vec![(
                AtomId(0),
                AtomId(1),
                EntitySpan::Unchanged(NoncovalentBondForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::stereo_atom_site(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Unchanged(AtomForm::default())],
            stereo_atoms: vec![(
                AtomId(1),
                Vec::new(),
                EntitySpan::Unchanged(StereoAtomForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::stereo_atom_ligand_union(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Unchanged(AtomForm::default())],
            stereo_atoms: vec![(
                AtomId(0),
                vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
                EntitySpan::Unchanged(StereoAtomForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::stereo_bond_site_union(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Unchanged(AtomForm::default())],
            stereo_bonds: vec![(
                BondId(0),
                Vec::new(),
                EntitySpan::Unchanged(StereoBondForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Bond(BondId(0)) },
    )]
    #[case::stereo_bond_ligand_union(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Unchanged(AtomForm::default())],
            bonds: vec![(AtomId(0), AtomId(0), EntitySpan::Unchanged(BondForm::default()))],
            stereo_bonds: vec![(
                BondId(0),
                vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
                EntitySpan::Unchanged(StereoBondForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::constraint_union(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Unchanged(AtomForm::default())],
            constraints: vec![ConstraintSpan::Unchanged(Constraint::Molecule(
                MoleculeConstraint::Connected { atoms: Some(vec![AtomId(1)]) },
            ))],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::bond_lhs(
        ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::default()),
                EntitySpan::Added(AtomForm::default()),
            ],
            bonds: vec![(
                AtomId(0),
                AtomId(1),
                EntitySpan::Unchanged(BondForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::bond_rhs(
        ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::default()),
                EntitySpan::Removed(AtomForm::default()),
            ],
            bonds: vec![(
                AtomId(0),
                AtomId(1),
                EntitySpan::Unchanged(BondForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::dative_rhs(
        ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::default()),
                EntitySpan::Removed(AtomForm::default()),
            ],
            dative: vec![(
                vec![AtomId(0)],
                AtomId(1),
                EntitySpan::Unchanged(DativeBondForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::aromatic_lhs(
        ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::default()),
                EntitySpan::Added(AtomForm::default()),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1)],
                EntitySpan::Unchanged(AromaticSystemForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::multicenter_rhs(
        ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::default()),
                EntitySpan::Unchanged(AtomForm::default()),
                EntitySpan::Removed(AtomForm::default()),
            ],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                EntitySpan::Unchanged(MulticenterBondForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(2)) },
    )]
    #[case::noncovalent_lhs(
        ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::default()),
                EntitySpan::Added(AtomForm::default()),
            ],
            noncovalent: vec![(
                AtomId(0),
                AtomId(1),
                EntitySpan::Unchanged(NoncovalentBondForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::stereo_atom_site_lhs(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Added(AtomForm::default())],
            stereo_atoms: vec![(
                AtomId(0),
                Vec::new(),
                EntitySpan::Unchanged(StereoAtomForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(0)) },
    )]
    #[case::stereo_atom_ligand_rhs(
        ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::default()),
                EntitySpan::Removed(AtomForm::default()),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
                EntitySpan::Unchanged(StereoAtomForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(1)) },
    )]
    #[case::stereo_bond_site_lhs(
        ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::default()),
                EntitySpan::Unchanged(AtomForm::default()),
            ],
            bonds: vec![(
                AtomId(0),
                AtomId(1),
                EntitySpan::Added(BondForm::default()),
            )],
            stereo_bonds: vec![(
                BondId(0),
                Vec::new(),
                EntitySpan::Unchanged(StereoBondForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Bond(BondId(0)) },
    )]
    #[case::stereo_bond_ligand_rhs(
        ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::default()),
                EntitySpan::Unchanged(AtomForm::default()),
                EntitySpan::Removed(AtomForm::default()),
            ],
            bonds: vec![(
                AtomId(0),
                AtomId(1),
                EntitySpan::Unchanged(BondForm::default()),
            )],
            stereo_bonds: vec![(
                BondId(0),
                vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom)],
                EntitySpan::Unchanged(StereoBondForm::default()),
            )],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(2)) },
    )]
    #[case::constraint_lhs(
        ReactionSpanEntries {
            atoms: vec![EntitySpan::Added(AtomForm::default())],
            constraints: vec![ConstraintSpan::Unchanged(Constraint::Atom(
                AtomId(0),
                AtomConstraintAst::valence(NumForm::Lit(4)),
            ))],
            ..Default::default()
        },
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(0)) },
    )]
    fn test_reaction_span_ast_try_from_entries_error(
        #[case] entries: ReactionSpanEntries,
        #[case] expected: ReactionSpanEntriesError,
    ) {
        assert_eq!(ReactionSpanAst::try_from_entries(entries), Err(expected));
    }

    #[rstest]
    #[case::atom(
        Constraint::Atom(AtomId(0), AtomConstraintAst::valence(NumForm::Lit(4))),
        Entity::Atom(AtomId(0))
    )]
    #[case::bond(
        Constraint::Bond(BondId(0), BondConstraintAst::aromatic(false)),
        Entity::Bond(BondId(0))
    )]
    #[case::dative_bond(
        Constraint::DativeBond(DativeBondId(0), DativeBondConstraintAst::aromatic(false)),
        Entity::DativeBond(DativeBondId(0))
    )]
    #[case::aromatic_system(
        Constraint::AromaticSystem(
            AromaticSystemId(0),
            AromaticSystemConstraintAst::electron_count(NumForm::Lit(2)),
        ),
        Entity::AromaticSystem(AromaticSystemId(0))
    )]
    #[case::multicenter_bond(
        Constraint::MulticenterBond(
            MulticenterBondId(0),
            MulticenterBondConstraintAst::electron_count(NumForm::Lit(2)),
        ),
        Entity::MulticenterBond(MulticenterBondId(0))
    )]
    #[case::noncovalent_bond(
        Constraint::NoncovalentBond(
            NoncovalentBondId(0),
            NoncovalentBondConstraintAst::intramolecular(true),
        ),
        Entity::NoncovalentBond(NoncovalentBondId(0))
    )]
    #[case::stereo_atom(
        Constraint::StereoAtom(
            StereoAtomId(0),
            StereoKind::Tetrahedral,
            StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined),
        ),
        Entity::StereoAtom(StereoAtomId(0))
    )]
    #[case::stereo_bond(
        Constraint::StereoBond(
            StereoBondId(0),
            StereoKind::CisTrans,
            StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Undetermined),
        ),
        Entity::StereoBond(StereoBondId(0))
    )]
    fn test_reaction_span_ast_try_from_entries_constraint_error(
        #[case] constraint: Constraint,
        #[case] entity: Entity,
    ) {
        let entries = ReactionSpanEntries {
            constraints: vec![ConstraintSpan::Unchanged(constraint)],
            ..Default::default()
        };

        assert_eq!(
            ReactionSpanAst::try_from_entries(entries),
            Err(ReactionSpanEntriesError::InvalidReference { entity }),
        );
    }

    #[rstest]
    #[case::invalid(
        ReactionSpanEntriesError::InvalidReference { entity: Entity::Atom(AtomId(3)) },
        "reaction span entries reference unavailable atom 3",
    )]
    fn test_reaction_span_entries_error_display(
        #[case] error: ReactionSpanEntriesError,
        #[case] expected: &str,
    ) {
        assert_eq!(error.to_string(), expected);
    }

    #[rstest]
    fn test_reaction_ast_to_reaction_span() {
        assert_eq!(
            ReactionAst::new(
                MoleculeAst::from_entries(MoleculeEntries {
                    atoms: vec![
                        AtomForm::from_element(Element::C),
                        AtomForm::from_element(Element::C),
                        AtomForm::from_element(Element::F),
                        AtomForm::from_element(Element::Cl),
                    ],
                    bonds: vec![
                        (AtomId(0), AtomId(3), BondForm::from_order(1)),
                        (AtomId(0), AtomId(1), BondForm::from_order(1)),
                        (AtomId(1), AtomId(2), BondForm::from_order(1)),
                    ],
                    dative: vec![(vec![AtomId(2)], AtomId(1), DativeBondForm::default(),)],
                    aromatic: vec![(
                        vec![AtomId(0), AtomId(1), AtomId(2)],
                        AromaticSystemForm::default(),
                    )],
                    multicenter: vec![(
                        vec![AtomId(0), AtomId(1), AtomId(2)],
                        MulticenterBondForm::default(),
                    )],
                    noncovalent: vec![(AtomId(0), AtomId(2), NoncovalentBondForm::default(),)],
                    stereo_atoms: vec![(
                        AtomId(2),
                        vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        StereoAtomForm::default(),
                    )],
                    stereo_bonds: vec![(
                        BondId(2),
                        vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        StereoBondForm::default(),
                    )],
                    constraints: Constraints::from_iter([
                        Constraint::Molecule(MoleculeConstraint::Connected {
                            atoms: Some(vec![AtomId(0), AtomId(1)]),
                        }),
                        Constraint::Molecule(MoleculeConstraint::Connected {
                            atoms: Some(vec![AtomId(1), AtomId(2)]),
                        }),
                    ]),
                }),
                Deltas::from_iter([
                    Delta::Atom(AtomDelta::Remove {
                        id: AtomId(2),
                        ast: AtomForm::from_element(Element::F),
                    }),
                    Delta::Atom(AtomDelta::Add {
                        id: AtomId(4),
                        ast: AtomForm::from_element(Element::O),
                    }),
                    Delta::Bond(BondDelta::Remove {
                        id: BondId(2),
                        atoms: [AtomId(1), AtomId(2)],
                        ast: BondForm::from_order(1),
                    }),
                    Delta::Bond(BondDelta::Add {
                        id: BondId(3),
                        atoms: [AtomId(1), AtomId(4)],
                        ast: BondForm::from_order(1),
                    }),
                    Delta::DativeBond(DativeBondDelta::Remove {
                        id: DativeBondId(0),
                        donors: vec![AtomId(2)],
                        acceptor: AtomId(1),
                        ast: DativeBondForm::default(),
                    }),
                    Delta::DativeBond(DativeBondDelta::Add {
                        id: DativeBondId(1),
                        donors: vec![AtomId(4)],
                        acceptor: AtomId(1),
                        ast: DativeBondForm::default(),
                    }),
                    Delta::AromaticSystem(AromaticSystemDelta::Remove {
                        id: AromaticSystemId(0),
                        atoms: vec![AtomId(0), AtomId(1), AtomId(2)],
                        ast: AromaticSystemForm::default(),
                    }),
                    Delta::AromaticSystem(AromaticSystemDelta::Add {
                        id: AromaticSystemId(1),
                        atoms: vec![AtomId(0), AtomId(1), AtomId(4)],
                        ast: AromaticSystemForm::default(),
                    }),
                    Delta::MulticenterBond(MulticenterBondDelta::Remove {
                        id: MulticenterBondId(0),
                        atoms: vec![AtomId(0), AtomId(1), AtomId(2)],
                        ast: MulticenterBondForm::default(),
                    }),
                    Delta::MulticenterBond(MulticenterBondDelta::Add {
                        id: MulticenterBondId(1),
                        atoms: vec![AtomId(0), AtomId(1), AtomId(4)],
                        ast: MulticenterBondForm::default(),
                    }),
                    Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                        id: NoncovalentBondId(0),
                        atoms: [AtomId(0), AtomId(2)],
                        ast: NoncovalentBondForm::default(),
                    }),
                    Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                        id: NoncovalentBondId(1),
                        atoms: [AtomId(0), AtomId(4)],
                        ast: NoncovalentBondForm::default(),
                    }),
                    Delta::StereoAtom(StereoAtomDelta::Remove {
                        id: StereoAtomId(0),
                        site: AtomId(2),
                        ligands: vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        ast: StereoAtomForm::default(),
                    }),
                    Delta::StereoAtom(StereoAtomDelta::Add {
                        id: StereoAtomId(1),
                        site: AtomId(4),
                        ligands: vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        ast: StereoAtomForm::default(),
                    }),
                    Delta::StereoBond(StereoBondDelta::Remove {
                        id: StereoBondId(0),
                        site: BondId(2),
                        ligands: vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        ast: StereoBondForm::default(),
                    }),
                    Delta::StereoBond(StereoBondDelta::Add {
                        id: StereoBondId(1),
                        site: BondId(3),
                        ligands: vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        ast: StereoBondForm::default(),
                    }),
                    Delta::Constraint(ConstraintDelta::Remove(Constraint::Molecule(
                        MoleculeConstraint::Connected {
                            atoms: Some(vec![AtomId(1), AtomId(2)]),
                        },
                    ))),
                    Delta::Constraint(ConstraintDelta::Add(Constraint::Molecule(
                        MoleculeConstraint::Connected {
                            atoms: Some(vec![AtomId(1), AtomId(4)]),
                        },
                    ))),
                ]),
            )
            .to_reaction_span(),
            Ok(ReactionSpanAst::from_entries(ReactionSpanEntries {
                atoms: vec![
                    EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                    EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                    EntitySpan::Removed(AtomForm::from_element(Element::F)),
                    EntitySpan::Unchanged(AtomForm::from_element(Element::Cl)),
                    EntitySpan::Added(AtomForm::from_element(Element::O)),
                ],
                bonds: vec![
                    (
                        AtomId(0),
                        AtomId(3),
                        EntitySpan::Unchanged(BondForm::from_order(1)),
                    ),
                    (
                        AtomId(0),
                        AtomId(1),
                        EntitySpan::Unchanged(BondForm::from_order(1)),
                    ),
                    (
                        AtomId(1),
                        AtomId(2),
                        EntitySpan::Removed(BondForm::from_order(1)),
                    ),
                    (
                        AtomId(1),
                        AtomId(4),
                        EntitySpan::Added(BondForm::from_order(1)),
                    ),
                ],
                dative: vec![
                    (
                        vec![AtomId(2)],
                        AtomId(1),
                        EntitySpan::Removed(DativeBondForm::default()),
                    ),
                    (
                        vec![AtomId(4)],
                        AtomId(1),
                        EntitySpan::Added(DativeBondForm::default()),
                    ),
                ],
                aromatic: vec![
                    (
                        vec![AtomId(0), AtomId(1), AtomId(2)],
                        EntitySpan::Removed(AromaticSystemForm::default()),
                    ),
                    (
                        vec![AtomId(0), AtomId(1), AtomId(4)],
                        EntitySpan::Added(AromaticSystemForm::default()),
                    ),
                ],
                multicenter: vec![
                    (
                        vec![AtomId(0), AtomId(1), AtomId(2)],
                        EntitySpan::Removed(MulticenterBondForm::default()),
                    ),
                    (
                        vec![AtomId(0), AtomId(1), AtomId(4)],
                        EntitySpan::Added(MulticenterBondForm::default()),
                    ),
                ],
                noncovalent: vec![
                    (
                        AtomId(0),
                        AtomId(2),
                        EntitySpan::Removed(NoncovalentBondForm::default()),
                    ),
                    (
                        AtomId(0),
                        AtomId(4),
                        EntitySpan::Added(NoncovalentBondForm::default()),
                    ),
                ],
                stereo_atoms: vec![
                    (
                        AtomId(2),
                        vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        EntitySpan::Removed(StereoAtomForm::default()),
                    ),
                    (
                        AtomId(4),
                        vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        EntitySpan::Added(StereoAtomForm::default()),
                    ),
                ],
                stereo_bonds: vec![
                    (
                        BondId(2),
                        vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        EntitySpan::Removed(StereoBondForm::default()),
                    ),
                    (
                        BondId(3),
                        vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        EntitySpan::Added(StereoBondForm::default()),
                    ),
                ],
                constraints: vec![
                    ConstraintSpan::Unchanged(Constraint::Molecule(
                        MoleculeConstraint::Connected {
                            atoms: Some(vec![AtomId(0), AtomId(1)]),
                        },
                    )),
                    ConstraintSpan::Removed(Constraint::Molecule(MoleculeConstraint::Connected {
                        atoms: Some(vec![AtomId(1), AtomId(2)]),
                    },)),
                    ConstraintSpan::Added(Constraint::Molecule(MoleculeConstraint::Connected {
                        atoms: Some(vec![AtomId(1), AtomId(4)]),
                    },)),
                ],
            })),
        );
    }

    #[rstest]
    #[case::constraint(
        ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
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
                    ast: AtomForm::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondForm::from_order(1),
                }),
                Delta::Constraint(ConstraintDelta::Remove(Constraint::Molecule(
                    MoleculeConstraint::Connected {
                        atoms: Some(vec![AtomId(0), AtomId(1)]),
                    },
                ))),
            ]),
        ),
        ReactionSpanAst::from_entries(ReactionSpanEntries {
            atoms: vec![
                EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                EntitySpan::Removed(AtomForm::from_element(Element::O)),
            ],
            bonds: vec![(
                AtomId(0),
                AtomId(1),
                EntitySpan::Removed(BondForm::from_order(1)),
            )],
            constraints: vec![ConstraintSpan::Removed(Constraint::Molecule(
                MoleculeConstraint::Connected {
                    atoms: Some(vec![AtomId(0), AtomId(1)]),
                },
            ))],
            ..Default::default()
        }),
    )]
    fn test_reaction_ast_to_reaction_span_constraint(
        #[case] reaction: ReactionAst,
        #[case] expected: ReactionSpanAst,
    ) {
        assert_eq!(reaction.to_reaction_span(), Ok(expected));
    }

    #[rstest]
    fn test_reaction_ast_to_reaction_span_constraint_error() {
        let constraint = Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(0), AtomId(1)]),
        });
        let reaction = ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                constraints: Constraints::from(constraint),
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    ast: AtomForm::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondForm::from_order(1),
                }),
            ]),
        );

        assert_eq!(reaction.to_reaction_span(), Err(Contradiction));
    }

    #[rstest]
    fn test_reaction_ast_to_reaction_span_error() {
        assert_eq!(
            ReactionAst::new(
                MoleculeAst::from_entries(MoleculeEntries {
                    atoms: vec![
                        AtomForm::from_element(Element::C),
                        AtomForm::from_element(Element::C),
                    ],
                    bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                    ..Default::default()
                }),
                Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                    id: BondId(0),
                    change: BondFieldChange::Order {
                        old: NumForm::Lit(2),
                        new: NumForm::Lit(3),
                    },
                })]),
            )
            .to_reaction_span(),
            Err(Contradiction),
        );
    }

    #[rstest]
    fn test_reaction_span_ast_right(substitution_reaction: ReactionAst) {
        let span = substitution_reaction.to_reaction_span().unwrap();
        assert_eq!(
            span.atoms(),
            [
                EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                EntitySpan::Removed(AtomForm::from_element(Element::O)),
                EntitySpan::Added(AtomForm::from_element(Element::N)),
            ],
        );
        assert_eq!(
            span.bonds(),
            [
                EntitySpan::Removed(BondForm::from_order(1)),
                EntitySpan::Added(BondForm::from_order(1)),
            ],
        );
        assert_eq!(
            span.rhs(),
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::N),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
        );
    }

    #[rstest]
    fn test_reaction_span_ast_left(substitution_reaction: ReactionAst) {
        let span = substitution_reaction.to_reaction_span().unwrap();
        assert_eq!(
            span.lhs(),
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
        );
    }

    #[rstest]
    #[case::order_change(
        ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::C),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            })]),
        ),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
            ..Default::default()
        }),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
    )]
    #[case::substitution(
        ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    ast: AtomForm::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondForm::from_order(1),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    ast: AtomForm::from_element(Element::N),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(1),
                    atoms: [AtomId(0), AtomId(2)],
                    ast: BondForm::from_order(1),
                }),
            ]),
        ),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
    )]
    #[case::stereo_atom(
        ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::F),
                    AtomForm::from_element(Element::Cl),
                    AtomForm::from_element(Element::Br),
                    AtomForm::from_element(Element::I),
                ],
                stereo_atoms: vec![(
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
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
                ast: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            })]),
        ),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
            ],
            bonds: vec![],
            ..Default::default()
        }),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
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
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                constraints: Constraints::from(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })),
                ..Default::default()
            }),
            Deltas::new(),
        ),
        vec![ConstraintSpan::Unchanged(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }))],
    )]
    #[case::added(
        ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::C)],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
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
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
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
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::C)],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(
            Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
        ))]),
    ))]
    #[case::remove(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            constraints: Constraints::from(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })),
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Constraint(ConstraintDelta::Remove(
            Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
        ))]),
    ))]
    #[case::dative_add(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::B),
                AtomForm::from_element(Element::N),
            ],
            bonds: vec![],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::DativeBond(DativeBondDelta::Add {
            id: DativeBondId(0),
            donors: vec![AtomId(0), AtomId(2)],
            acceptor: AtomId(1),
            ast: DativeBondForm::from_order(1),
        })]),
    ))]
    #[case::aromatic_add(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::C)],
            bonds: vec![],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Add {
            id: AromaticSystemId(0),
            atoms: vec![AtomId(0), AtomId(1)],
            ast: AromaticSystemForm::from_electrons(vec![1, 2]),
        })]),
    ))]
    #[case::multicenter_add(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::B),
                AtomForm::from_element(Element::H),
                AtomForm::from_element(Element::B),
            ],
            bonds: vec![],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::MulticenterBond(MulticenterBondDelta::Add {
            id: MulticenterBondId(0),
            atoms: vec![AtomId(0), AtomId(1), AtomId(2)],
            ast: MulticenterBondForm::from_electrons(vec![3, 5, 7]),
        })]),
    ))]
    #[case::noncovalent_add(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)],
            bonds: vec![],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(0),
            atoms: [AtomId(0), AtomId(1)],
            ast: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        })]),
    ))]
    #[case::noncovalent_remove(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))],
            constraints: Constraints::new(),
            ..Default::default()
        }),
        Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
            id: NoncovalentBondId(0),
            atoms: [AtomId(0), AtomId(1)],
            ast: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        })]),
    ))]
    #[case::noncovalent_modify(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))],
            constraints: Constraints::new(),
            ..Default::default()
        }),
        Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField {
            id: NoncovalentBondId(0),
            change: NoncovalentBondFieldChange::Kind {
                old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
            },
        })]),
    ))]
    #[case::stereo_atom_add(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
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
            ast: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        })]),
    ))]
    #[case::stereo_atom_remove(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
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
            ast: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        })]),
    ))]
    #[case::stereo_atom_modify(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            )],
            constraints: Constraints::new(),
            ..Default::default()
        }),
        Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
            id: StereoAtomId(0),
            change: StereoAtomFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            },
        })]),
    ))]
    #[case::stereo_bond_add(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
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
            ast: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        })]),
    ))]
    fn test_reaction_span_ast_to_reaction(#[case] reaction: ReactionAst) {
        assert_eq!(
            reaction.clone().to_reaction_span().unwrap().to_reaction(),
            reaction,
        );
    }

    #[rstest]
    #[case::unchanged(
        ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)],
                noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::new(),
        ),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))],
            constraints: Constraints::new(),
            ..Default::default()
        }),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))],
            constraints: Constraints::new(),
            ..Default::default()
        }),
    )]
    #[case::added(
        ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)],
                bonds: vec![],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                id: NoncovalentBondId(0),
                atoms: [AtomId(0), AtomId(1)],
                ast: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            })]),
        ),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)],
            bonds: vec![],
            ..Default::default()
        }),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))],
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
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    ast: AtomForm::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondForm::from_order(1),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    ast: AtomForm::from_element(Element::N),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(1),
                    atoms: [AtomId(0), AtomId(2)],
                    ast: BondForm::from_order(1),
                }),
            ]),
        )
    }

    #[rstest]
    fn test_reaction_span_ast_superimpose() {
        let left = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
            ],
            bonds: vec![
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
            ],
            dative: vec![(vec![AtomId(2)], AtomId(1), DativeBondForm::default())],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::default(),
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondForm::default(),
            )],
            noncovalent: vec![(AtomId(0), AtomId(2), NoncovalentBondForm::default())],
            stereo_atoms: vec![(
                AtomId(2),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                ],
                StereoAtomForm::default(),
            )],
            stereo_bonds: vec![(
                BondId(2),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                ],
                StereoBondForm::default(),
            )],
            constraints: Constraints::from_iter([
                Constraint::Molecule(MoleculeConstraint::Connected {
                    atoms: Some(vec![AtomId(0), AtomId(1)]),
                }),
                Constraint::Molecule(MoleculeConstraint::Connected {
                    atoms: Some(vec![AtomId(1), AtomId(2)]),
                }),
            ]),
        });
        let right = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(3), BondForm::from_order(1)),
            ],
            dative: vec![(vec![AtomId(3)], AtomId(1), DativeBondForm::default())],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(3)],
                AromaticSystemForm::default(),
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(3)],
                MulticenterBondForm::default(),
            )],
            noncovalent: vec![(AtomId(0), AtomId(3), NoncovalentBondForm::default())],
            stereo_atoms: vec![(
                AtomId(3),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                ],
                StereoAtomForm::default(),
            )],
            stereo_bonds: vec![(
                BondId(2),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                ],
                StereoBondForm::default(),
            )],
            constraints: Constraints::from_iter([
                Constraint::Molecule(MoleculeConstraint::Connected {
                    atoms: Some(vec![AtomId(0), AtomId(1)]),
                }),
                Constraint::Molecule(MoleculeConstraint::Connected {
                    atoms: Some(vec![AtomId(1), AtomId(3)]),
                }),
            ]),
        });
        let atoms = Correspondence::new(
            vec![
                (AtomId(0), AtomId(0)),
                (AtomId(1), AtomId(1)),
                (AtomId(3), AtomId(2)),
            ],
            4,
            4,
        )
        .expect("correspondence producer preserves partial-bijection invariants");
        let correspondence = MoleculeCorrespondence::induce(&left, &right, atoms)
            .expect("the atom correspondence describes the molecule pair");

        assert_eq!(
            ReactionSpanAst::superimpose(&left, &right, &correspondence),
            Some(ReactionSpanAst::from_entries(ReactionSpanEntries {
                atoms: vec![
                    EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                    EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                    EntitySpan::Removed(AtomForm::from_element(Element::F)),
                    EntitySpan::Unchanged(AtomForm::from_element(Element::Cl)),
                    EntitySpan::Added(AtomForm::from_element(Element::O)),
                ],
                bonds: vec![
                    (
                        AtomId(0),
                        AtomId(3),
                        EntitySpan::Unchanged(BondForm::from_order(1)),
                    ),
                    (
                        AtomId(0),
                        AtomId(1),
                        EntitySpan::Unchanged(BondForm::from_order(1)),
                    ),
                    (
                        AtomId(1),
                        AtomId(2),
                        EntitySpan::Removed(BondForm::from_order(1)),
                    ),
                    (
                        AtomId(1),
                        AtomId(4),
                        EntitySpan::Added(BondForm::from_order(1)),
                    ),
                ],
                dative: vec![
                    (
                        vec![AtomId(2)],
                        AtomId(1),
                        EntitySpan::Removed(DativeBondForm::default()),
                    ),
                    (
                        vec![AtomId(4)],
                        AtomId(1),
                        EntitySpan::Added(DativeBondForm::default()),
                    ),
                ],
                aromatic: vec![
                    (
                        vec![AtomId(0), AtomId(1), AtomId(2)],
                        EntitySpan::Removed(AromaticSystemForm::default()),
                    ),
                    (
                        vec![AtomId(0), AtomId(1), AtomId(4)],
                        EntitySpan::Added(AromaticSystemForm::default()),
                    ),
                ],
                multicenter: vec![
                    (
                        vec![AtomId(0), AtomId(1), AtomId(2)],
                        EntitySpan::Removed(MulticenterBondForm::default()),
                    ),
                    (
                        vec![AtomId(0), AtomId(1), AtomId(4)],
                        EntitySpan::Added(MulticenterBondForm::default()),
                    ),
                ],
                noncovalent: vec![
                    (
                        AtomId(0),
                        AtomId(2),
                        EntitySpan::Removed(NoncovalentBondForm::default()),
                    ),
                    (
                        AtomId(0),
                        AtomId(4),
                        EntitySpan::Added(NoncovalentBondForm::default()),
                    ),
                ],
                stereo_atoms: vec![
                    (
                        AtomId(2),
                        vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        EntitySpan::Removed(StereoAtomForm::default()),
                    ),
                    (
                        AtomId(4),
                        vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        EntitySpan::Added(StereoAtomForm::default()),
                    ),
                ],
                stereo_bonds: vec![
                    (
                        BondId(2),
                        vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        EntitySpan::Removed(StereoBondForm::default()),
                    ),
                    (
                        BondId(3),
                        vec![
                            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        ],
                        EntitySpan::Added(StereoBondForm::default()),
                    ),
                ],
                constraints: vec![
                    ConstraintSpan::Unchanged(Constraint::Molecule(
                        MoleculeConstraint::Connected {
                            atoms: Some(vec![AtomId(0), AtomId(1)]),
                        },
                    )),
                    ConstraintSpan::Removed(Constraint::Molecule(MoleculeConstraint::Connected {
                        atoms: Some(vec![AtomId(1), AtomId(2)]),
                    },)),
                    ConstraintSpan::Added(Constraint::Molecule(MoleculeConstraint::Connected {
                        atoms: Some(vec![AtomId(1), AtomId(4)]),
                    },)),
                ],
            })),
        );
    }

    #[rstest]
    #[case::count(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        }),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        }),
        MoleculeCorrespondence::new(
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
        ),
    )]
    #[case::bond(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(2), BondForm::from_order(1))],
            ..Default::default()
        }),
        MoleculeCorrespondence::new(
            Correspondence::new(
                vec![(AtomId(0), AtomId(0)), (AtomId(1), AtomId(1))],
                2,
                3,
            )
            .unwrap(),
            Correspondence::new(vec![(BondId(0), BondId(0))], 1, 1).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
        ),
    )]
    #[case::aromatic_system(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1)],
                AromaticSystemForm::default(),
            )],
            ..Default::default()
        }),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(2)],
                AromaticSystemForm::default(),
            )],
            ..Default::default()
        }),
        MoleculeCorrespondence::new(
            Correspondence::new(
                vec![(AtomId(0), AtomId(0)), (AtomId(1), AtomId(1))],
                2,
                3,
            )
            .unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(
                vec![(AromaticSystemId(0), AromaticSystemId(0))],
                1,
                1,
            )
            .unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
        ),
    )]
    #[case::stereo_atom(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
                StereoAtomForm::default(),
            )],
            ..Default::default()
        }),
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom)],
                StereoAtomForm::default(),
            )],
            ..Default::default()
        }),
        MoleculeCorrespondence::new(
            Correspondence::new(
                vec![
                    (AtomId(0), AtomId(0)),
                    (AtomId(1), AtomId(1)),
                    (AtomId(2), AtomId(2)),
                ],
                3,
                3,
            )
            .unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![(StereoAtomId(0), StereoAtomId(0))], 1, 1).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
        ),
    )]
    fn test_reaction_span_ast_superimpose_invalid_context(
        #[case] lhs: MoleculeAst,
        #[case] rhs: MoleculeAst,
        #[case] correspondence: MoleculeCorrespondence,
    ) {
        assert_eq!(
            ReactionSpanAst::superimpose(&lhs, &rhs, &correspondence),
            None,
        );
        assert_eq!(lhs.difference_to(&rhs, &correspondence), None);
    }

    #[rstest]
    fn test_reaction_span_ast_superimpose_narrow_correspondence() {
        let lhs = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        });
        let rhs = lhs.clone();
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::new(vec![(AtomId(0), AtomId(0)), (AtomId(1), AtomId(1))], 2, 2)
                .unwrap(),
            Correspondence::new(vec![], 1, 1).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
            Correspondence::new(vec![], 0, 0).unwrap(),
        );

        assert_eq!(
            ReactionSpanAst::superimpose(&lhs, &rhs, &correspondence),
            Some(ReactionSpanAst::from_entries(ReactionSpanEntries {
                atoms: vec![
                    EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                    EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
                ],
                bonds: vec![
                    (
                        AtomId(0),
                        AtomId(1),
                        EntitySpan::Removed(BondForm::from_order(1)),
                    ),
                    (
                        AtomId(0),
                        AtomId(1),
                        EntitySpan::Added(BondForm::from_order(1)),
                    ),
                ],
                ..Default::default()
            })),
        );
    }

    #[rstest]
    fn test_reaction_span_ast_correspondence() {
        // atom 0 unchanged, 1 modified (C→N), 2 removed (left) with 2 added (right O): all four
        // EntitySpan variants in the atom column.
        let left = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        });
        let right = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        });
        let atoms = Correspondence::new(vec![(AtomId(0), AtomId(0)), (AtomId(1), AtomId(1))], 3, 3)
            .expect("correspondence producer preserves partial-bijection invariants");
        let correspondence = MoleculeCorrespondence::induce(&left, &right, atoms)
            .expect("the atom correspondence describes the molecule pair");
        let span = ReactionSpanAst::superimpose(&left, &right, &correspondence).unwrap();

        // recovers the input correspondence, and inverts `superimpose`.
        assert_eq!(span.correspondence(), correspondence);
        assert_eq!(
            ReactionSpanAst::superimpose(&span.lhs(), &span.rhs(), &span.correspondence()),
            Some(span)
        );
    }

    #[rstest]
    fn test_molecule_ast_difference_to() {
        // C-C (order 1) → C-C (order 2), total correspondence: a single bond-order modify.
        let left = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        });
        let right = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
            ..Default::default()
        });
        let atoms = Correspondence::new(vec![(AtomId(0), AtomId(0)), (AtomId(1), AtomId(1))], 2, 2)
            .expect("correspondence producer preserves partial-bijection invariants");
        let correspondence = MoleculeCorrespondence::induce(&left, &right, atoms)
            .expect("the atom correspondence describes the molecule pair");
        assert_eq!(
            left.difference_to(&right, &correspondence),
            Some(Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            })])),
        );
    }
}
