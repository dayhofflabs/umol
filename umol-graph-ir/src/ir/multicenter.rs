//! Multicenter bonds: the molecule's collection and one bond's attribute form.

use std::borrow::Cow;
use std::sync::Arc;

use umol_graph_core::{
    Compaction, GraphCompaction, GraphCorrespondence, GraphRemapping, NodeId, ParticipantPosition,
    RelationId, VarRelationSet,
};
use umol_graph_ir_macros::{Lattice, Normalize};
use umol_perm::DynPermutation;

use super::constraint::{MulticenterBondConstraintForm, MulticenterBondConstraintsForm};
use super::delta::EntitySpan;
use super::electrons::ElectronCountsForm;
use super::error::Contradiction;
use super::frame::MulticenterBondsFrameAction;
use super::id::{AtomId, MulticenterBondId};
use super::num::NumForm;
use super::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
use super::traits::{FrameTransport, Lattice, Normalize, Reframe};

/// A multicenter bond's atom frame and attributes.
///
/// Normalization preserves atom order. Reframing sorts atoms and transports
/// their electron contributions into that frame.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MulticenterBondEntry<'a> {
    pub(crate) atoms: Cow<'a, [AtomId]>,
    pub(crate) attributes: Cow<'a, MulticenterBondForm>,
}

impl Normalize for MulticenterBondEntry<'_> {
    fn normalize(mut self) -> Result<Self, Contradiction> {
        self.attributes = match self.attributes {
            Cow::Borrowed(attributes) => attributes.normalized()?,
            Cow::Owned(attributes) => Cow::Owned(attributes.normalize()?),
        };
        Ok(self)
    }
}

impl FrameTransport for MulticenterBondEntry<'_> {
    type Action = DynPermutation;

    fn reframe_by(mut self, action: &Self::Action) -> Option<Self> {
        self.atoms = Cow::Owned(action.act(&self.atoms)?);
        self.attributes = Cow::Owned(self.attributes.into_owned().reframe_by(action)?);
        Some(self)
    }
}

impl Reframe for MulticenterBondEntry<'_> {
    fn representative_action(&self) -> Self::Action {
        multicenter_bond_representative_action(self.atoms.to_vec())
    }

    fn reframe(self) -> Result<Self, Contradiction> {
        let action = self.representative_action();
        self.normalize()?
            .reframe_by(&action)
            .ok_or(Contradiction)?
            .normalize()
    }

    fn framed_eq(&self, other: &Self) -> bool {
        if self == other {
            return true;
        }
        if DynPermutation::between(&other.atoms, &self.atoms)
            .and_then(|action| other.attributes.as_ref().clone().reframe_by(&action))
            .is_some_and(|attributes| self.attributes.normalized_eq(&attributes))
        {
            return true;
        }
        self.clone().reframe() == other.clone().reframe()
    }
}

/// The molecule's multicenter bonds.
///
/// The atoms bear the participant frame: the per-member electron counts of
/// [`MulticenterBondForm`] are read against it, position by position. Values are issued by checked
/// molecule construction and trusted graph-IR transformations; raw assembly is not a public
/// construction path.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct MulticenterBonds(Arc<VarRelationSet<NodeId, MulticenterBondForm>>);

impl MulticenterBonds {
    pub(crate) fn new(entries: Vec<(Vec<AtomId>, MulticenterBondForm)>) -> Self {
        Self(Arc::new(VarRelationSet::new(
            entries
                .into_iter()
                .map(|(atoms, attributes)| {
                    (atoms.into_iter().map(NodeId::from).collect(), attributes)
                })
                .collect(),
        )))
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: MulticenterBondId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = MulticenterBondId> {
        self.0.ids().map(MulticenterBondId::from)
    }

    /// The atoms of `id`, in their stored frame.
    pub fn atoms(&self, id: MulticenterBondId) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.0
            .participants(RelationId::from(id))
            .iter()
            .map(|&atom| AtomId::from(atom))
    }

    pub fn attributes(&self, id: MulticenterBondId) -> &MulticenterBondForm {
        self.0.data(RelationId::from(id))
    }

    /// The complete bond value, borrowing its attributes and collecting atom ids.
    ///
    /// # Panics
    ///
    /// Panics if `id` is not a bond in this set.
    pub(crate) fn entry(&self, id: MulticenterBondId) -> MulticenterBondEntry<'_> {
        MulticenterBondEntry {
            atoms: Cow::Owned(self.atoms(id).collect()),
            attributes: Cow::Borrowed(self.attributes(id)),
        }
    }

    /// Ids of the multicenter bonds `atom` belongs to. Unlike aromatic systems these may overlap,
    /// so an atom can belong to several; integrity rejects only identical atom sets.
    pub fn incident_ids(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = MulticenterBondId> + '_ {
        self.0
            .incident_to_node(NodeId::from(atom))
            .iter()
            .map(|&id| MulticenterBondId::from(id))
    }

    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.0.has_incident_to_node(NodeId::from(atom))
    }

    /// Whether bond `id` is the one over `atoms` — the known-id sibling of
    /// [`coincident_id`](Self::coincident_id).
    pub fn is_coincident(&self, id: MulticenterBondId, atoms: &[AtomId]) -> bool {
        let query: Vec<NodeId> = atoms.iter().map(|&atom| NodeId::from(atom)).collect();
        self.0.is_coincident(RelationId::from(id), &query)
    }

    /// Id of the entity coinciding with these participants — the one whose participants equal
    /// them as a multiset. The identity question, distinct from lookup.
    pub fn coincident_id(&self, atoms: &[AtomId]) -> Option<MulticenterBondId> {
        // Multicenter bonds anchor on their atoms, so the node index is the one to scan.
        let query: Vec<NodeId> = atoms.iter().map(|&atom| NodeId::from(atom)).collect();
        let anchor = *query.first()?;
        self.0
            .coincident_to_node(anchor, &query)
            .map(MulticenterBondId::from)
    }

    /// The atoms of `id` as graph nodes, for graph-core interop that is not yet typed in graph-IR
    /// ids. The public accessor is [`Self::atoms`].
    pub(crate) fn atom_nodes(&self, id: MulticenterBondId) -> &[NodeId] {
        self.0.participants(RelationId::from(id))
    }

    pub(crate) fn attributes_mut(&mut self, id: MulticenterBondId) -> &mut MulticenterBondForm {
        Arc::make_mut(&mut self.0).data_mut(RelationId::from(id))
    }

    pub(crate) fn attributes_iter_mut(
        &mut self,
    ) -> impl ExactSizeIterator<Item = &mut MulticenterBondForm> {
        Arc::make_mut(&mut self.0)
            .iter_mut()
            .map(|(_, _, attributes)| attributes)
    }

    pub(crate) fn add(
        &mut self,
        atoms: &[AtomId],
        attributes: MulticenterBondForm,
    ) -> MulticenterBondId {
        let nodes: Vec<NodeId> = atoms.iter().copied().map(NodeId::from).collect();
        Arc::make_mut(&mut self.0).add(&nodes, attributes).into()
    }

    /// Append entries in input order and return their ids without retaining a borrow.
    pub(crate) fn extend(
        &mut self,
        entries: Vec<(&[AtomId], MulticenterBondForm)>,
    ) -> impl ExactSizeIterator<Item = MulticenterBondId> + use<> {
        let start = self.count();
        if !entries.is_empty() {
            let nodes: Vec<NodeId> = entries
                .iter()
                .flat_map(|(atoms, _)| atoms.iter().copied().map(NodeId::from))
                .collect();
            let mut remaining = nodes.as_slice();
            let entries = entries
                .into_iter()
                .map(|(atoms, attributes)| {
                    let (atoms, rest) = remaining.split_at(atoms.len());
                    remaining = rest;
                    (atoms, attributes)
                })
                .collect();
            let _ = Arc::make_mut(&mut self.0).extend(entries);
        }
        (start..self.count()).map(MulticenterBondId::from)
    }

    pub(crate) fn remove(&mut self, ids: &[MulticenterBondId]) {
        if ids.is_empty() {
            return;
        }
        let ids: Vec<RelationId> = ids.iter().copied().map(RelationId::from).collect();
        Arc::make_mut(&mut self.0).remove(&ids);
    }

    pub(crate) fn tracked_remove(
        &mut self,
        ids: &[MulticenterBondId],
    ) -> Compaction<MulticenterBondId> {
        if ids.is_empty() {
            return Compaction::identity(self.count());
        }
        let ids: Vec<RelationId> = ids.iter().copied().map(RelationId::from).collect();
        let compaction = Arc::make_mut(&mut self.0).tracked_remove(&ids);
        Compaction::new(
            compaction.source_count(),
            compaction
                .removed()
                .iter()
                .copied()
                .map(MulticenterBondId::from)
                .collect(),
        )
        .expect("relation compaction contains valid multicenter bond ids")
    }

    /// Reinsert saved entries at their original ids after restoring surviving topology ids.
    pub(crate) fn restore(
        &mut self,
        compaction: &Compaction<MulticenterBondId>,
        removed: Vec<(MulticenterBondId, Vec<AtomId>, MulticenterBondForm)>,
    ) {
        if compaction.removed().is_empty() {
            return;
        }
        let relations = Compaction::new(
            compaction.source_count(),
            compaction
                .removed()
                .iter()
                .copied()
                .map(RelationId::from)
                .collect(),
        )
        .expect("multicenter bond compaction contains valid relation ids");
        let removed = removed
            .into_iter()
            .map(|(id, atoms, attributes)| {
                (
                    RelationId::from(id),
                    atoms.into_iter().map(NodeId::from).collect(),
                    attributes,
                )
            })
            .collect();
        Arc::make_mut(&mut self.0).restore(&relations, removed);
    }

    /// Restore the original atom ids in surviving entries, preserving entry ids and attributes.
    pub(crate) fn restore_topology_ids(&mut self, compaction: &GraphCompaction) {
        if compaction.nodes().removed().is_empty() && compaction.edges().removed().is_empty() {
            return;
        }
        Arc::make_mut(&mut self.0).restore_participants(compaction);
    }

    pub(crate) fn replace_atoms(&mut self, id: MulticenterBondId, atoms: &[AtomId]) {
        let nodes: Vec<NodeId> = atoms.iter().copied().map(NodeId::from).collect();
        Arc::make_mut(&mut self.0).replace_participants(id.into(), &nodes);
    }

    pub(crate) fn replace_atom(&mut self, id: MulticenterBondId, position: usize, atom: AtomId) {
        let position = ParticipantPosition(position as u32);
        Arc::make_mut(&mut self.0).replace_participant(id.into(), position, atom.into());
    }

    pub(crate) fn insert_atom(&mut self, id: MulticenterBondId, position: usize, atom: AtomId) {
        let position = ParticipantPosition(position as u32);
        Arc::make_mut(&mut self.0).insert_participant(id.into(), position, atom.into());
    }

    pub(crate) fn remove_atom(&mut self, id: MulticenterBondId, position: usize) {
        let position = ParticipantPosition(position as u32);
        Arc::make_mut(&mut self.0).remove_participant(id.into(), position);
    }

    pub(crate) fn compact(
        &self,
        compaction: &GraphCompaction,
    ) -> (Self, Compaction<MulticenterBondId>) {
        let (set, relations) = self.0.tracked_compact(compaction);
        let relations = Compaction::new(
            relations.source_count(),
            relations
                .removed()
                .iter()
                .copied()
                .map(MulticenterBondId::from)
                .collect(),
        )
        .expect("relation compaction contains valid multicenter bond ids");
        (Self(Arc::new(set)), relations)
    }

    /// Map participant references, preserving entity ids, row order, attributes, and frames.
    ///
    /// # Panics
    /// Panics if a referenced node or edge has no image in `correspondence`.
    pub fn map(&self, correspondence: &GraphCorrespondence) -> Self {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Map participant references, or return `None` if any reference has no image.
    /// Unreferenced nodes and edges need not have images. No entity is dropped.
    pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self> {
        Some(Self(Arc::new(self.0.try_map(correspondence)?)))
    }

    pub(crate) fn remap(&self, remapping: &GraphRemapping) -> Self {
        Self(Arc::new(self.0.remap(remapping)))
    }

    /// Glue `right`, relabelled into this molecule's id space, onto `self`: coinciding bonds meet,
    /// non-coinciding bonds are carried. `None` when a coincident meet is bottom.
    pub(crate) fn glue(&self, right: &Self, correspondence: &GraphCorrespondence) -> Option<Self> {
        self.0
            .pushout(
                &right.map(correspondence).0,
                // Multicenter bonds anchor on their atoms: the node index.
                |set, atoms| {
                    atoms
                        .first()
                        .and_then(|&node| set.coincident_to_node(node, atoms))
                },
                |(left_atoms, left), (right_atoms, right)| {
                    let left_atoms: Vec<AtomId> =
                        left_atoms.iter().map(|&atom| AtomId::from(atom)).collect();
                    let right_atoms: Vec<AtomId> =
                        right_atoms.iter().map(|&atom| AtomId::from(atom)).collect();
                    let action = DynPermutation::between(&right_atoms, &left_atoms)?;
                    right.clone().reframe_by(&action)?.meet(left)
                },
            )
            .map(|object| Self(Arc::new(object)))
    }

    pub(crate) fn into_entries(self) -> Vec<(Vec<AtomId>, MulticenterBondForm)> {
        Arc::try_unwrap(self.0)
            .unwrap_or_else(|shared| (*shared).clone())
            .into_entries()
            .into_iter()
            .map(|(atoms, attributes)| (atoms.into_iter().map(AtomId::from).collect(), attributes))
            .collect()
    }
}

pub(crate) fn reframe_multicenter_bonds_with(
    mut multicenter_bonds: MulticenterBonds,
    mut visit: impl FnMut(MulticenterBondId, &DynPermutation),
) -> Result<MulticenterBonds, Contradiction> {
    let set = Arc::make_mut(&mut multicenter_bonds.0);
    for relation_id in set.ids().collect::<Vec<_>>() {
        let id = MulticenterBondId::from(relation_id);
        let stored = set
            .participants(relation_id)
            .iter()
            .map(|&atom| AtomId::from(atom))
            .collect();
        let action = multicenter_bond_representative_action(stored);
        let attributes = set.data(relation_id).clone().normalize()?;
        *set.data_mut(relation_id) = attributes
            .reframe_by(&action)
            .ok_or(Contradiction)?
            .normalize()?;
        set.permute_participants(relation_id, &participant_order(&action));
        visit(id, &action);
    }
    Ok(multicenter_bonds)
}

impl Normalize for MulticenterBonds {
    fn normalize(mut self) -> Result<Self, Contradiction> {
        for attributes in self.attributes_iter_mut() {
            *attributes = attributes.clone().normalize()?;
        }
        Ok(self)
    }
}

impl FrameTransport for MulticenterBonds {
    type Action = MulticenterBondsFrameAction;

    fn reframe_by(mut self, actions: &Self::Action) -> Option<Self> {
        let set = Arc::make_mut(&mut self.0);
        for relation_id in set.ids().collect::<Vec<_>>() {
            let action = actions.action(MulticenterBondId::from(relation_id))?;
            if action.degree() != set.participants(relation_id).len() {
                return None;
            }
            *set.data_mut(relation_id) = set.data(relation_id).clone().reframe_by(action)?;
            set.permute_participants(relation_id, &participant_order(action));
        }
        Some(self)
    }
}

impl Reframe for MulticenterBonds {
    fn representative_action(&self) -> Self::Action {
        let actions = self
            .ids()
            .map(|id| multicenter_bond_representative_action(self.atoms(id).collect()))
            .collect();
        MulticenterBondsFrameAction::from_vec(actions)
            .expect("every dynamic permutation is a multicenter-bond action")
    }

    fn reframe(self) -> Result<Self, Contradiction> {
        reframe_multicenter_bonds_with(self, |_, _| {})
    }
}

/// The reaction span's multicenter bonds, one [`EntitySpan`] per entity against a single participant frame.
///
/// The `Molecule` peer is [`MulticenterBonds`]. The surface is deliberately duplicated rather than shared
/// through a payload parameter: a type parameter on the molecule-level aggregates would complicate
/// the primary carrier to serve this one. Values are issued by
/// [`ReactionSpan`](super::reaction_span::ReactionSpan).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct MulticenterBondSpans(VarRelationSet<NodeId, EntitySpan<MulticenterBondForm>>);

impl MulticenterBondSpans {
    pub(crate) fn into_entries(self) -> Vec<(Vec<AtomId>, EntitySpan<MulticenterBondForm>)> {
        self.0
            .into_entries()
            .into_iter()
            .map(|(atoms, span)| (atoms.into_iter().map(AtomId::from).collect(), span))
            .collect()
    }

    pub(crate) fn new(entries: Vec<(Vec<AtomId>, EntitySpan<MulticenterBondForm>)>) -> Self {
        Self(VarRelationSet::new(
            entries
                .into_iter()
                .map(|(atoms, span)| (atoms.into_iter().map(NodeId::from).collect(), span))
                .collect(),
        ))
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: MulticenterBondId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = MulticenterBondId> {
        self.0.ids().map(MulticenterBondId::from)
    }

    /// The atoms of `id`, in their stored frame.
    pub fn atoms(&self, id: MulticenterBondId) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.0
            .participants(RelationId::from(id))
            .iter()
            .map(|&atom| AtomId::from(atom))
    }

    pub fn attributes(&self, id: MulticenterBondId) -> &EntitySpan<MulticenterBondForm> {
        self.0.data(RelationId::from(id))
    }

    /// Map participant references, preserving entity ids, row order, attributes, and frames.
    ///
    /// # Panics
    /// Panics if a referenced node or edge has no image in `correspondence`.
    pub fn map(&self, correspondence: &GraphCorrespondence) -> Self {
        self.try_map(correspondence)
            .expect("correspondence must cover every participant reference")
    }

    /// Map participant references, or return `None` if any reference has no image.
    /// Unreferenced nodes and edges need not have images. No entity is dropped.
    pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self> {
        self.0.try_map(correspondence).map(Self)
    }

    pub(crate) fn remap(&self, remapping: &GraphRemapping) -> Self {
        Self(self.0.remap(remapping))
    }
}

impl Normalize for MulticenterBondSpans {
    fn normalize(mut self) -> Result<Self, Contradiction> {
        for id in self.0.ids().collect::<Vec<_>>() {
            *self.0.data_mut(id) = self.0.data(id).clone().normalize()?;
        }
        Ok(self)
    }
}

impl FrameTransport for MulticenterBondSpans {
    type Action = MulticenterBondsFrameAction;

    fn reframe_by(mut self, actions: &Self::Action) -> Option<Self> {
        for relation_id in self.0.ids().collect::<Vec<_>>() {
            let action = actions.action(MulticenterBondId::from(relation_id))?;
            if action.degree() != self.0.participants(relation_id).len() {
                return None;
            }
            *self.0.data_mut(relation_id) = self.0.data(relation_id).clone().reframe_by(action)?;
            self.0
                .permute_participants(relation_id, &participant_order(action));
        }
        Some(self)
    }
}

impl Reframe for MulticenterBondSpans {
    fn representative_action(&self) -> Self::Action {
        let actions = self
            .ids()
            .map(|id| multicenter_bond_representative_action(self.atoms(id).collect()))
            .collect();
        MulticenterBondsFrameAction::from_vec(actions)
            .expect("every dynamic permutation is a multicenter-bond action")
    }

    fn reframe(self) -> Result<Self, Contradiction> {
        reframe_multicenter_bond_spans_with(self, |_, _| {})
    }
}

pub(crate) fn reframe_multicenter_bond_spans_with(
    mut multicenter_bonds: MulticenterBondSpans,
    mut visit: impl FnMut(MulticenterBondId, &DynPermutation),
) -> Result<MulticenterBondSpans, Contradiction> {
    for relation_id in multicenter_bonds.0.ids().collect::<Vec<_>>() {
        let id = MulticenterBondId::from(relation_id);
        let stored = multicenter_bonds
            .0
            .participants(relation_id)
            .iter()
            .map(|&atom| AtomId::from(atom))
            .collect();
        let action = multicenter_bond_representative_action(stored);
        let span = multicenter_bonds.0.data(relation_id).clone().normalize()?;
        *multicenter_bonds.0.data_mut(relation_id) =
            span.reframe_by(&action).ok_or(Contradiction)?.normalize()?;
        multicenter_bonds
            .0
            .permute_participants(relation_id, &participant_order(&action));
        visit(id, &action);
    }
    Ok(multicenter_bonds)
}

pub(crate) fn multicenter_bond_representative_action(frame: Vec<AtomId>) -> DynPermutation {
    let mut image: Vec<usize> = (0..frame.len()).collect();
    image.sort_unstable_by_key(|&position| frame[position]);
    DynPermutation::try_from(image).expect("sorted positions form a permutation")
}

fn participant_order(action: &DynPermutation) -> Vec<ParticipantPosition> {
    action
        .image()
        .iter()
        .map(|&position| ParticipantPosition(position as u32))
        .collect()
}

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Normalize, Lattice)]
pub struct MulticenterBondForm {
    pub electrons: ElectronCountsForm,
    pub charge: NumForm,
    pub unpaired_electrons: UnpairedElectronsForm,
    pub constraints: MulticenterBondConstraintsForm,
}

/// Attribute update for a multicenter bond. Ordinary fields are optional,
/// unpaired-electron components are updated independently, and undetermined constraints remove
/// their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterBondUpdate {
    pub electrons: Option<ElectronCountsForm>,
    pub charge: Option<NumForm>,
    pub unpaired_electrons: UnpairedElectronsUpdate,
    pub constraints: MulticenterBondConstraintsForm,
}

impl From<&str> for MulticenterBondForm {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid multicenter bond string")
    }
}

impl MulticenterBondForm {
    /// Concrete: every inherent field is ground; the constraint channel does
    /// not bear on concreteness.
    pub fn is_concrete(&self) -> bool {
        let MulticenterBondForm {
            electrons,
            charge,
            unpaired_electrons,
            constraints: _,
        } = self;
        electrons.is_ground() && charge.is_ground() && unpaired_electrons.is_ground()
    }
    pub fn new(electrons: ElectronCountsForm) -> Self {
        Self {
            electrons,
            ..Default::default()
        }
    }

    pub fn from_electrons(electrons: Vec<i64>) -> Self {
        Self::new(ElectronCountsForm::Lit(electrons))
    }

    pub fn with_charge(mut self, charge: impl Into<NumForm>) -> Self {
        self.charge = charge.into();
        self
    }

    pub fn with_unpaired_electrons(
        mut self,
        unpaired_electrons: impl Into<UnpairedElectronsForm>,
    ) -> Self {
        self.unpaired_electrons = unpaired_electrons.into();
        self
    }

    /// Add a single constraint, replacing any existing entry of the same
    /// kind (last-wins per `MulticenterBondConstraintsForm::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<MulticenterBondConstraintForm>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `MulticenterBondConstraintsForm::set`).
    /// Does not clear existing constraints; use `bond.constraints.clear()`
    /// or direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<MulticenterBondConstraintForm>,
    {
        for c in constraints {
            self.constraints.set(c.into());
        }
        self
    }

    /// Fill `Undetermined` value-bearing struct fields with zero defaults:
    /// charge → `Lit(0)`, unpaired electrons → closed-shell singlet `(0, 1)`. `electrons`
    /// and `constraints` are preserved. The result is concrete iff `electrons`
    /// is already `Lit`.
    pub fn into_concrete(mut self) -> Self {
        if self.charge.is_undetermined() {
            self.charge = NumForm::Lit(0);
        }
        if self.unpaired_electrons.is_undetermined() {
            self.unpaired_electrons = UnpairedElectronsForm::from((0_u8, 1_u8));
        }
        self
    }

    /// Apply an attribute update, leaving omitted leaves and constraint keys unchanged.
    pub fn update(&self, update: &MulticenterBondUpdate) -> MulticenterBondForm {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        MulticenterBondForm {
            electrons: update
                .electrons
                .clone()
                .unwrap_or_else(|| self.electrons.clone()),
            charge: update.charge.clone().unwrap_or_else(|| self.charge.clone()),
            unpaired_electrons: self.unpaired_electrons.update(&update.unpaired_electrons),
            constraints,
        }
    }

    /// Derive the minimal normalized attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> MulticenterBondUpdate {
        let mut constraints = MulticenterBondConstraintsForm::new();
        for new in other.constraints.iter() {
            if self
                .constraints
                .get(new.key())
                .is_none_or(|old| !old.normalized_eq(new))
            {
                constraints.set(new.clone());
            }
        }
        for old in self.constraints.iter() {
            if other.constraints.get(old.key()).is_none() {
                constraints.set(old.as_undetermined());
            }
        }
        MulticenterBondUpdate {
            electrons: (!self.electrons.normalized_eq(&other.electrons))
                .then(|| other.electrons.clone()),
            charge: (!self.charge.normalized_eq(&other.charge)).then(|| other.charge.clone()),
            unpaired_electrons: self
                .unpaired_electrons
                .difference_to(&other.unpaired_electrons),
            constraints,
        }
    }

    /// Reorder the positional `electrons` by `order`, tracking a participant
    /// reordering; charge / unpaired electrons / constraints are positionless and unchanged.
    pub fn permute(&mut self, order: &[ParticipantPosition]) {
        self.electrons.permute(order);
    }
}

impl FrameTransport for MulticenterBondForm {
    type Action = DynPermutation;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        let Self {
            electrons,
            charge,
            unpaired_electrons,
            constraints,
        } = self;
        Some(Self {
            electrons: electrons.reframe_by(action)?,
            charge,
            unpaired_electrons,
            constraints: constraints.reframe_by(action)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use std::ptr;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::{Correspondence, EdgeId};

    use super::*;
    use crate::ir::error::Contradiction;
    use crate::ir::traits::Normalize;

    #[rstest]
    #[case::singleton(NumForm::lit_set([1]), Ok(NumForm::Lit(1)))]
    #[case::contradiction(NumForm::lit_set([]), Err(Contradiction))]
    fn test_multicenter_bond_entry_normalize(
        #[case] charge: NumForm,
        #[case] expected: Result<NumForm, Contradiction>,
        #[values(false, true)] owned: bool,
    ) {
        let atoms = [AtomId(2), AtomId(0), AtomId(1)];
        let attributes = MulticenterBondForm::from_electrons(vec![3, 1, 2]).with_charge(charge);
        let entry = MulticenterBondEntry {
            atoms: if owned {
                Cow::Owned(atoms.to_vec())
            } else {
                Cow::Borrowed(&atoms)
            },
            attributes: if owned {
                Cow::Owned(attributes.clone())
            } else {
                Cow::Borrowed(&attributes)
            },
        };
        let expected = expected.map(|charge| MulticenterBondEntry {
            atoms: Cow::Borrowed(&atoms),
            attributes: Cow::Owned(
                MulticenterBondForm::from_electrons(vec![3, 1, 2]).with_charge(charge),
            ),
        });
        let normalized = entry.normalize();
        assert_eq!(normalized, expected);
        if let Ok(normalized) = normalized {
            assert_eq!(normalized.clone().normalize(), Ok(normalized));
        }
    }

    #[rstest]
    #[case::swap(vec![1, 0, 2])]
    #[case::cycle(vec![1, 2, 0])]
    #[case::reverse(vec![2, 1, 0])]
    fn test_multicenter_bond_entry_reframe_by(#[case] image: Vec<usize>) {
        let atoms = [AtomId(2), AtomId(0), AtomId(1)];
        let counts = [3, 1, 2];
        let attributes = MulticenterBondForm::from_electrons(counts.to_vec()).with_charge(1);
        let entry = MulticenterBondEntry {
            atoms: Cow::Borrowed(&atoms),
            attributes: Cow::Borrowed(&attributes),
        };
        let action = DynPermutation::try_from(image.clone()).unwrap();
        let expected = MulticenterBondEntry {
            atoms: Cow::Owned(image.iter().map(|&i| atoms[i]).collect()),
            attributes: Cow::Owned(
                MulticenterBondForm::from_electrons(image.iter().map(|&i| counts[i]).collect())
                    .with_charge(1),
            ),
        };
        assert_eq!(
            entry.clone().reframe_by(&DynPermutation::identity(3)),
            Some(entry.clone())
        );
        let transformed = entry.clone().reframe_by(&action).unwrap();
        assert_eq!(transformed, expected);
        assert_eq!(
            transformed.reframe_by(&action.inverse()),
            Some(entry.clone())
        );
        let second = DynPermutation::try_from(vec![1, 0, 2]).unwrap();
        assert_eq!(
            entry
                .clone()
                .reframe_by(&action)
                .unwrap()
                .reframe_by(&second),
            entry.reframe_by(&action.compose(&second).unwrap()),
        );
    }

    #[rstest]
    #[case::short_action(vec![1, 2, 3], 2)]
    #[case::long_action(vec![1, 2, 3], 4)]
    #[case::short_counts(vec![1, 2], 3)]
    #[case::long_counts(vec![1, 2, 3, 4], 3)]
    fn test_multicenter_bond_entry_reframe_by_error(
        #[case] counts: Vec<i64>,
        #[case] degree: usize,
    ) {
        let entry = MulticenterBondEntry {
            atoms: Cow::Borrowed(&[AtomId(0), AtomId(1), AtomId(2)]),
            attributes: Cow::Owned(MulticenterBondForm::from_electrons(counts)),
        };
        assert_eq!(entry.reframe_by(&DynPermutation::identity(degree)), None);
    }

    #[rstest]
    fn test_multicenter_bond_entry_reframe() {
        let attributes =
            MulticenterBondForm::from_electrons(vec![3, 1, 2]).with_charge(NumForm::lit_set([1]));
        let entry = MulticenterBondEntry {
            atoms: Cow::Borrowed(&[AtomId(2), AtomId(0), AtomId(1)]),
            attributes: Cow::Borrowed(&attributes),
        };
        let expected = MulticenterBondEntry {
            atoms: Cow::Owned(vec![AtomId(0), AtomId(1), AtomId(2)]),
            attributes: Cow::Owned(
                MulticenterBondForm::from_electrons(vec![1, 2, 3]).with_charge(1),
            ),
        };
        let action = DynPermutation::try_from(vec![1, 2, 0]).unwrap();
        assert_eq!(entry.representative_action(), action);
        assert_eq!(entry.clone().reframe(), Ok(expected.clone()));
        assert_eq!(entry.tracked_reframe(), Ok((expected.clone(), action)));
        assert_eq!(expected.clone().reframe(), Ok(expected));
    }

    #[rstest]
    #[case::stored(vec![AtomId(0), AtomId(1), AtomId(2)], vec![1, 2, 3], true, true)]
    #[case::reordered(vec![AtomId(2), AtomId(0), AtomId(1)], vec![3, 1, 2], false, true)]
    #[case::misaligned(vec![AtomId(2), AtomId(0), AtomId(1)], vec![1, 2, 3], false, false)]
    #[case::counts(vec![AtomId(0), AtomId(1), AtomId(2)], vec![1, 2, 4], false, false)]
    #[case::membership(vec![AtomId(0), AtomId(1), AtomId(3)], vec![1, 2, 3], false, false)]
    #[case::length(vec![AtomId(0), AtomId(1)], vec![1, 2], false, false)]
    fn test_multicenter_bond_entry_framed_eq(
        #[case] atoms: Vec<AtomId>,
        #[case] counts: Vec<i64>,
        #[case] normalized: bool,
        #[case] framed: bool,
    ) {
        let left = MulticenterBondEntry {
            atoms: Cow::Borrowed(&[AtomId(0), AtomId(1), AtomId(2)]),
            attributes: Cow::Owned(MulticenterBondForm::from_electrons(vec![1, 2, 3])),
        };
        let right = MulticenterBondEntry {
            atoms: Cow::Owned(atoms),
            attributes: Cow::Owned(MulticenterBondForm::from_electrons(counts)),
        };
        assert_eq!(left.normalized_eq(&right), normalized);
        assert_eq!(left.framed_eq(&right), framed);
        assert_eq!(right.framed_eq(&left), framed);
        assert_eq!(
            left.framed_eq(&right),
            left == right || left.reframe() == right.reframe()
        );
    }

    #[rstest]
    #[case::stored(vec![AtomId(0), AtomId(1), AtomId(2)], true)]
    #[case::reordered(vec![AtomId(2), AtomId(0), AtomId(1)], true)]
    #[case::membership(vec![AtomId(0), AtomId(1), AtomId(3)], false)]
    #[case::length(vec![AtomId(0), AtomId(1)], false)]
    fn test_multicenter_bond_entry_framed_eq_undetermined(
        #[case] atoms: Vec<AtomId>,
        #[case] expected: bool,
    ) {
        let attributes = MulticenterBondForm::default();
        let left = MulticenterBondEntry {
            atoms: Cow::Borrowed(&[AtomId(0), AtomId(1), AtomId(2)]),
            attributes: Cow::Borrowed(&attributes),
        };
        let right = MulticenterBondEntry {
            atoms: Cow::Owned(atoms),
            attributes: Cow::Borrowed(&attributes),
        };
        assert_eq!(left.framed_eq(&right), expected);
    }

    #[rstest]
    #[case::reordered(vec![AtomId(2), AtomId(0), AtomId(1)], vec![1, 2, 3], vec![3, 1, 2])]
    #[case::misaligned(vec![AtomId(2), AtomId(0), AtomId(1)], vec![1, 2, 3], vec![1, 2, 3])]
    #[case::membership(vec![AtomId(2), AtomId(0), AtomId(3)], vec![1, 2, 3], vec![3, 1, 2])]
    #[case::one_count_mismatch(vec![AtomId(2), AtomId(0), AtomId(1)], vec![1, 2, 3], vec![3, 1])]
    #[case::both_count_mismatches(vec![AtomId(2), AtomId(0), AtomId(1)], vec![1, 2], vec![3, 1])]
    fn test_multicenter_bond_entry_framed_eq_definition(
        #[case] atoms: Vec<AtomId>,
        #[case] left_counts: Vec<i64>,
        #[case] right_counts: Vec<i64>,
        #[values(NumForm::Lit(0), NumForm::lit_set([]))] left_charge: NumForm,
        #[values(NumForm::lit_set([0]), NumForm::lit_set([]))] right_charge: NumForm,
    ) {
        let left = MulticenterBondEntry {
            atoms: Cow::Borrowed(&[AtomId(0), AtomId(1), AtomId(2)]),
            attributes: Cow::Owned(
                MulticenterBondForm::from_electrons(left_counts).with_charge(left_charge),
            ),
        };
        let right = MulticenterBondEntry {
            atoms: Cow::Owned(atoms),
            attributes: Cow::Owned(
                MulticenterBondForm::from_electrons(right_counts).with_charge(right_charge),
            ),
        };
        let expected = left == right || left.clone().reframe() == right.clone().reframe();
        assert_eq!(left.framed_eq(&right), expected);
        assert_eq!(right.framed_eq(&left), expected);
    }

    #[rstest]
    fn test_multicenter_bonds_entry() {
        let bonds = MulticenterBonds::new(vec![(
            vec![AtomId(2), AtomId(1)],
            MulticenterBondForm::from_electrons(vec![2, 1]),
        )]);
        let before = bonds.clone();
        let entry = bonds.entry(MulticenterBondId(0));
        assert_eq!(
            entry,
            MulticenterBondEntry {
                atoms: Cow::Borrowed(&[AtomId(2), AtomId(1)]),
                attributes: Cow::Owned(MulticenterBondForm::from_electrons(vec![2, 1])),
            }
        );
        assert!(matches!(entry.atoms, Cow::Owned(_)));
        assert!(
            matches!(entry.attributes, Cow::Borrowed(attributes) if ptr::eq(attributes, bonds.attributes(MulticenterBondId(0))))
        );
        let reframed = entry.reframe().unwrap();
        assert_eq!(
            reframed,
            MulticenterBondEntry {
                atoms: Cow::Borrowed(&[AtomId(1), AtomId(2)]),
                attributes: Cow::Owned(MulticenterBondForm::from_electrons(vec![1, 2])),
            }
        );
        assert_eq!(bonds, before);
    }

    #[rstest]
    #[should_panic]
    fn test_multicenter_bonds_entry_error() {
        MulticenterBonds::default().entry(MulticenterBondId(0));
    }

    #[rstest]
    fn test_multicenter_bonds_add() {
        let mut bonds = MulticenterBonds::default();
        let attributes = MulticenterBondForm::from_electrons(vec![1, 2, 3]);
        let id = bonds.add(&[AtomId(4), AtomId(1), AtomId(4)], attributes.clone());

        assert_eq!(id, MulticenterBondId(0));
        assert_eq!(
            bonds.atoms(id).collect::<Vec<_>>(),
            vec![AtomId(4), AtomId(1), AtomId(4)]
        );
        assert_eq!(bonds.attributes(id), &attributes);
        assert_eq!(bonds.incident_ids(AtomId(4)).collect::<Vec<_>>(), vec![id]);
        assert!(!bonds.has_incident(AtomId(2)));
    }

    #[rstest]
    #[case::empty(vec![], vec![MulticenterBondId(0), MulticenterBondId(1), MulticenterBondId(2)])]
    #[case::populated(vec![(vec![AtomId(0), AtomId(1)], MulticenterBondForm::from_electrons(vec![1, 1]))], vec![MulticenterBondId(1), MulticenterBondId(2), MulticenterBondId(3)])]
    fn test_multicenter_bonds_extend(
        #[case] initial: Vec<(Vec<AtomId>, MulticenterBondForm)>,
        #[case] expected_ids: Vec<MulticenterBondId>,
        #[values(false, true)] shared: bool,
    ) {
        let mut bonds = MulticenterBonds::new(initial.clone());
        let original = shared.then(|| bonds.clone());
        let mut expected = initial.clone();
        expected.extend([
            (
                vec![AtomId(4), AtomId(1), AtomId(4)],
                MulticenterBondForm::from_electrons(vec![1, 2, 3]),
            ),
            (vec![], MulticenterBondForm::default()),
            (
                vec![AtomId(1), AtomId(3)],
                MulticenterBondForm::from_electrons(vec![2, 1]),
            ),
        ]);
        let mut ids = {
            let first = [AtomId(4), AtomId(1), AtomId(4)];
            let last = [AtomId(1), AtomId(3)];
            bonds.extend(vec![
                (&first, MulticenterBondForm::from_electrons(vec![1, 2, 3])),
                (&[], MulticenterBondForm::default()),
                (&last, MulticenterBondForm::from_electrons(vec![2, 1])),
            ])
        };

        assert_eq!(bonds, MulticenterBonds::new(expected.clone()));
        for atom in (0..8).map(AtomId) {
            let expected_ids: Vec<_> = expected
                .iter()
                .enumerate()
                .filter(|(_, (atoms, _))| atoms.contains(&atom))
                .map(|(index, _)| MulticenterBondId::from(index))
                .collect();
            assert_eq!(bonds.incident_ids(atom).collect::<Vec<_>>(), expected_ids);
            assert_eq!(bonds.has_incident(atom), !expected_ids.is_empty());
        }

        if let Some(original) = original {
            assert_eq!(original, MulticenterBonds::new(initial));
        }
        drop(bonds);
        assert_eq!(ids.len(), 3);
        assert_eq!(ids.next(), Some(expected_ids[0]));
        assert_eq!(ids.len(), 2);
        assert_eq!(ids.collect::<Vec<_>>(), expected_ids[1..]);
    }

    #[rstest]
    #[case::empty(MulticenterBonds::default())]
    #[case::populated(MulticenterBonds::new(vec![(vec![AtomId(0), AtomId(1)], MulticenterBondForm::from_electrons(vec![1, 1]))]))]
    fn test_multicenter_bonds_extend_identity(#[case] original: MulticenterBonds) {
        let mut bonds = original.clone();
        let mut ids = bonds.extend(vec![]);

        assert_eq!(bonds, original);
        assert!(Arc::ptr_eq(&bonds.0, &original.0));
        assert_eq!(ids.len(), 0);
        assert_eq!(ids.next(), None);
    }

    #[rstest]
    fn test_multicenter_bonds_tracked_remove() {
        let original = MulticenterBonds::new(vec![
            (
                vec![AtomId(0)],
                MulticenterBondForm::from_electrons(vec![1]),
            ),
            (
                vec![AtomId(1)],
                MulticenterBondForm::from_electrons(vec![2]),
            ),
            (
                vec![AtomId(2)],
                MulticenterBondForm::from_electrons(vec![3]),
            ),
        ]);
        let mut bonds = original.clone();
        let compaction = bonds.tracked_remove(&[MulticenterBondId(2), MulticenterBondId(0)]);
        let mut plain = original.clone();
        plain.remove(&[MulticenterBondId(0), MulticenterBondId(2)]);

        assert_eq!(
            compaction,
            Compaction::new(3, vec![MulticenterBondId(0), MulticenterBondId(2)]).unwrap()
        );
        assert_eq!(bonds, plain);
        assert_eq!(
            bonds.atoms(MulticenterBondId(0)).collect::<Vec<_>>(),
            vec![AtomId(1)]
        );
        bonds.restore(
            &compaction,
            vec![
                (
                    MulticenterBondId(2),
                    vec![AtomId(2)],
                    original.attributes(MulticenterBondId(2)).clone(),
                ),
                (
                    MulticenterBondId(0),
                    vec![AtomId(0)],
                    original.attributes(MulticenterBondId(0)).clone(),
                ),
            ],
        );
        assert_eq!(bonds, original);
    }

    #[rstest]
    fn test_multicenter_bonds_replace_atoms() {
        let first = MulticenterBondForm::from_electrons(vec![1, 2]);
        let second = MulticenterBondForm::from_electrons(vec![3]);
        let mut bonds = MulticenterBonds::new(vec![
            (vec![AtomId(1), AtomId(2)], first.clone()),
            (vec![AtomId(3)], second.clone()),
        ]);
        bonds.replace_atoms(MulticenterBondId(0), &[AtomId(4), AtomId(4), AtomId(1)]);

        assert_eq!(
            bonds.atoms(MulticenterBondId(0)).collect::<Vec<_>>(),
            vec![AtomId(4), AtomId(4), AtomId(1)]
        );
        assert_eq!(bonds.attributes(MulticenterBondId(0)), &first);
        assert_eq!(
            bonds.atoms(MulticenterBondId(1)).collect::<Vec<_>>(),
            vec![AtomId(3)]
        );
        assert_eq!(bonds.attributes(MulticenterBondId(1)), &second);
        assert_eq!(
            bonds.incident_ids(AtomId(4)).collect::<Vec<_>>(),
            vec![MulticenterBondId(0)]
        );
        assert!(!bonds.has_incident(AtomId(2)));
    }

    #[rstest]
    fn test_multicenter_bonds_replace_atom() {
        let mut bonds = MulticenterBonds::new(vec![(
            vec![AtomId(1), AtomId(2)],
            MulticenterBondForm::default(),
        )]);
        bonds.replace_atom(MulticenterBondId(0), 1, AtomId(4));
        assert_eq!(
            bonds.atoms(MulticenterBondId(0)).collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4)]
        );
        assert!(!bonds.has_incident(AtomId(2)));
        assert!(bonds.has_incident(AtomId(4)));
    }

    #[rstest]
    fn test_multicenter_bonds_insert_atom() {
        let mut bonds = MulticenterBonds::new(vec![(
            vec![AtomId(1), AtomId(2)],
            MulticenterBondForm::default(),
        )]);
        bonds.insert_atom(MulticenterBondId(0), 1, AtomId(4));
        bonds.insert_atom(MulticenterBondId(0), 3, AtomId(5));
        assert_eq!(
            bonds.atoms(MulticenterBondId(0)).collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4), AtomId(2), AtomId(5)]
        );
        assert_eq!(
            bonds.incident_ids(AtomId(5)).collect::<Vec<_>>(),
            vec![MulticenterBondId(0)]
        );
    }

    #[rstest]
    fn test_multicenter_bonds_remove_atom() {
        let mut bonds = MulticenterBonds::new(vec![(
            vec![AtomId(1), AtomId(2)],
            MulticenterBondForm::default(),
        )]);
        bonds.remove_atom(MulticenterBondId(0), 0);
        bonds.remove_atom(MulticenterBondId(0), 0);
        assert_eq!(
            bonds.atoms(MulticenterBondId(0)).collect::<Vec<_>>(),
            Vec::<AtomId>::new()
        );
        assert!(!bonds.has_incident(AtomId(1)));
        assert!(!bonds.has_incident(AtomId(2)));
    }

    #[rstest]
    fn test_multicenter_bonds_compact() {
        let original = MulticenterBonds::new(vec![
            (
                vec![AtomId(0), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(1), AtomId(4)],
                MulticenterBondForm::from_electrons(vec![3, 4]),
            ),
            (
                vec![AtomId(3), AtomId(5)],
                MulticenterBondForm::from_electrons(vec![5, 6]),
            ),
        ]);
        let graph = GraphCompaction::new(
            Compaction::new(6, vec![NodeId(1)]).unwrap(),
            Compaction::identity(0),
        );
        let (mut compacted, rows) = original.compact(&graph);

        assert_eq!(
            rows,
            Compaction::new(3, vec![MulticenterBondId(1)]).unwrap()
        );
        assert_eq!(
            compacted.atoms(MulticenterBondId(0)).collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1)]
        );
        assert_eq!(
            compacted.atoms(MulticenterBondId(1)).collect::<Vec<_>>(),
            vec![AtomId(2), AtomId(4)]
        );
        assert_eq!(
            compacted.incident_ids(AtomId(1)).collect::<Vec<_>>(),
            vec![MulticenterBondId(0)]
        );

        compacted.restore_topology_ids(&graph);
        compacted.restore(
            &rows,
            vec![(
                MulticenterBondId(1),
                vec![AtomId(1), AtomId(4)],
                original.attributes(MulticenterBondId(1)).clone(),
            )],
        );
        assert_eq!(compacted, original);
    }

    #[rstest]
    #[case::covered(None, None)]
    #[case::missing_atom_0(Some(NodeId(0)), None)]
    #[case::missing_atom_2(Some(NodeId(2)), None)]
    #[case::missing_atom_4(Some(NodeId(4)), None)]
    #[case::missing_atom_6(Some(NodeId(6)), None)]

    fn test_multicenter_bonds_try_map(
        #[case] missing_node: Option<NodeId>,
        #[case] missing_edge: Option<EdgeId>,
    ) {
        let input = MulticenterBonds::new(vec![
            (
                vec![AtomId(4), AtomId(0), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![1, 2, 3]),
            ),
            (
                vec![AtomId(6), AtomId(2), AtomId(0)],
                MulticenterBondForm::from_electrons(vec![2, 3, 4]),
            ),
        ]);
        let correspondence = GraphCorrespondence::new(
            Correspondence::new(
                vec![
                    (NodeId(0), NodeId(5)),
                    (NodeId(2), NodeId(1)),
                    (NodeId(4), NodeId(7)),
                    (NodeId(6), NodeId(3)),
                ]
                .into_iter()
                .filter(|(id, _)| Some(*id) != missing_node)
                .collect(),
                8,
                9,
            )
            .unwrap(),
            Correspondence::new(
                vec![(EdgeId(0), EdgeId(4)), (EdgeId(2), EdgeId(1))]
                    .into_iter()
                    .filter(|(id, _)| Some(*id) != missing_edge)
                    .collect(),
                4,
                6,
            )
            .unwrap(),
        );
        let expected = if missing_node.is_none() && missing_edge.is_none() {
            Some(MulticenterBonds::new(vec![
                (
                    vec![AtomId(7), AtomId(5), AtomId(1)],
                    MulticenterBondForm::from_electrons(vec![1, 2, 3]),
                ),
                (
                    vec![AtomId(3), AtomId(1), AtomId(5)],
                    MulticenterBondForm::from_electrons(vec![2, 3, 4]),
                ),
            ]))
        } else {
            None
        };
        assert_eq!(input.try_map(&correspondence), expected);
        if let Some(expected) = expected {
            assert_eq!(input.map(&correspondence), expected);
        }
    }

    #[rstest]
    #[should_panic(expected = "correspondence must cover every participant reference")]
    fn test_multicenter_bonds_map_error() {
        let input = MulticenterBonds::new(vec![
            (
                vec![AtomId(4), AtomId(0), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![1, 2, 3]),
            ),
            (
                vec![AtomId(6), AtomId(2), AtomId(0)],
                MulticenterBondForm::from_electrons(vec![2, 3, 4]),
            ),
        ]);
        input.map(&GraphCorrespondence::new(
            Correspondence::empty(),
            Correspondence::empty(),
        ));
    }

    #[rstest]
    #[case::covered(None, None)]
    #[case::missing_atom_0(Some(NodeId(0)), None)]
    #[case::missing_atom_2(Some(NodeId(2)), None)]
    #[case::missing_atom_4(Some(NodeId(4)), None)]
    #[case::missing_atom_6(Some(NodeId(6)), None)]

    fn test_multicenter_bond_spans_try_map(
        #[case] missing_node: Option<NodeId>,
        #[case] missing_edge: Option<EdgeId>,
    ) {
        let input = MulticenterBondSpans::new(vec![
            (
                vec![AtomId(4), AtomId(0), AtomId(2)],
                EntitySpan::Modified {
                    lhs: MulticenterBondForm::from_electrons(vec![1, 2, 3]),
                    rhs: MulticenterBondForm::from_electrons(vec![2, 3, 4]),
                },
            ),
            (
                vec![AtomId(6), AtomId(2), AtomId(0)],
                EntitySpan::Added(MulticenterBondForm::from_electrons(vec![2, 3, 4])),
            ),
        ]);
        let correspondence = GraphCorrespondence::new(
            Correspondence::new(
                vec![
                    (NodeId(0), NodeId(5)),
                    (NodeId(2), NodeId(1)),
                    (NodeId(4), NodeId(7)),
                    (NodeId(6), NodeId(3)),
                ]
                .into_iter()
                .filter(|(id, _)| Some(*id) != missing_node)
                .collect(),
                8,
                9,
            )
            .unwrap(),
            Correspondence::new(
                vec![(EdgeId(0), EdgeId(4)), (EdgeId(2), EdgeId(1))]
                    .into_iter()
                    .filter(|(id, _)| Some(*id) != missing_edge)
                    .collect(),
                4,
                6,
            )
            .unwrap(),
        );
        let expected = if missing_node.is_none() && missing_edge.is_none() {
            Some(MulticenterBondSpans::new(vec![
                (
                    vec![AtomId(7), AtomId(5), AtomId(1)],
                    EntitySpan::Modified {
                        lhs: MulticenterBondForm::from_electrons(vec![1, 2, 3]),
                        rhs: MulticenterBondForm::from_electrons(vec![2, 3, 4]),
                    },
                ),
                (
                    vec![AtomId(3), AtomId(1), AtomId(5)],
                    EntitySpan::Added(MulticenterBondForm::from_electrons(vec![2, 3, 4])),
                ),
            ]))
        } else {
            None
        };
        assert_eq!(input.try_map(&correspondence), expected);
        if let Some(expected) = expected {
            assert_eq!(input.map(&correspondence), expected);
        }
    }

    #[rstest]
    #[should_panic(expected = "correspondence must cover every participant reference")]
    fn test_multicenter_bond_spans_map_error() {
        let input = MulticenterBondSpans::new(vec![
            (
                vec![AtomId(4), AtomId(0), AtomId(2)],
                EntitySpan::Modified {
                    lhs: MulticenterBondForm::from_electrons(vec![1, 2, 3]),
                    rhs: MulticenterBondForm::from_electrons(vec![2, 3, 4]),
                },
            ),
            (
                vec![AtomId(6), AtomId(2), AtomId(0)],
                EntitySpan::Added(MulticenterBondForm::from_electrons(vec![2, 3, 4])),
            ),
        ]);
        input.map(&GraphCorrespondence::new(
            Correspondence::empty(),
            Correspondence::empty(),
        ));
    }

    /// Both sides of a `Modified` span are read against one participant list, so one action carries
    /// both. Selection never consults the payload here, so the two sides cannot disagree about it.
    #[rstest]
    fn test_multicenter_bond_spans_reframe() {
        let mut spans = MulticenterBondSpans::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(vec![10, 20, 30]),
                rhs: MulticenterBondForm::from_electrons(vec![11, 21, 31]),
            },
        )]);
        spans.0.permute_participants(
            RelationId(0),
            &[
                ParticipantPosition(2),
                ParticipantPosition(0),
                ParticipantPosition(1),
            ],
        );

        let source = spans.clone();
        let (reframed, actions) = spans.tracked_reframe().expect("the forms are satisfiable");

        assert_eq!(
            reframed.atoms(MulticenterBondId(0)).collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4), AtomId(7)],
        );
        assert_eq!(
            reframed.attributes(MulticenterBondId(0)),
            &EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(vec![20, 30, 10]),
                rhs: MulticenterBondForm::from_electrons(vec![21, 31, 11]),
            },
        );
        assert_eq!(
            actions.action(MulticenterBondId(0)),
            Some(&DynPermutation::try_from(vec![1, 2, 0]).expect("expected action is valid")),
        );
        assert_eq!(source.reframe_by(&actions), Some(reframed));
    }

    #[rstest]
    fn test_multicenter_bond_spans_normalize() {
        let spans = MulticenterBondSpans::new(vec![(
            vec![AtomId(1), AtomId(2)],
            EntitySpan::Modified {
                lhs: MulticenterBondForm::default().with_charge(NumForm::lit_set([0])),
                rhs: MulticenterBondForm::default().with_charge(0),
            },
        )]);

        let normalized = spans.normalize().expect("the forms are satisfiable");

        assert_eq!(
            normalized.attributes(MulticenterBondId(0)),
            &EntitySpan::Unchanged(MulticenterBondForm::default().with_charge(0)),
        );
    }

    #[rstest]
    fn test_multicenter_bond_spans_reframe_identity() {
        let spans = MulticenterBondSpans::new(vec![(
            vec![AtomId(1), AtomId(4)],
            EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(vec![10, 20]),
                rhs: MulticenterBondForm::from_electrons(vec![11, 21]),
            },
        )]);
        let once = spans.reframe().expect("the forms are satisfiable");
        let twice = once.clone().reframe().expect("the forms are satisfiable");
        assert_eq!(twice, once);
    }

    /// A side that declines the frame change takes the whole span with it: one action serves both,
    /// so there is no partial result to keep. Here the rhs electron vector disagrees in length with
    /// the participant frame.
    #[rstest]
    fn test_multicenter_bond_spans_reframe_error() {
        let spans = MulticenterBondSpans::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(vec![10, 20, 30]),
                rhs: MulticenterBondForm::from_electrons(vec![11, 21]),
            },
        )]);
        assert_eq!(spans.reframe(), Err(Contradiction));
    }

    #[rstest]
    fn test_multicenter_bond_spans_framed_eq() {
        let mut unsorted = MulticenterBondSpans::new(vec![(
            vec![AtomId(1), AtomId(4)],
            EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(vec![10, 20]),
                rhs: MulticenterBondForm::from_electrons(vec![11, 21]),
            },
        )]);
        unsorted.0.permute_participants(
            RelationId(0),
            &[ParticipantPosition(1), ParticipantPosition(0)],
        );
        // The permute leaves the payload where it was, so the stored frame now reads
        // atom 4 -> 10, atom 1 -> 20; the selected presentation states the same fact sorted.
        let selected = MulticenterBondSpans::new(vec![(
            vec![AtomId(1), AtomId(4)],
            EntitySpan::Modified {
                lhs: MulticenterBondForm::from_electrons(vec![20, 10]),
                rhs: MulticenterBondForm::from_electrons(vec![21, 11]),
            },
        )]);

        assert!(unsorted.framed_eq(&selected));
        assert!(unsorted != selected);
    }

    /// Storage sorts an `Unordered` factor on construction, so the stored frame is permuted first
    /// to model the frame-preserving storage S5 introduces.
    #[fixture]
    fn unsorted_bond() -> MulticenterBonds {
        let mut systems = MulticenterBonds::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            MulticenterBondForm::from_electrons(vec![10, 20, 30]),
        )]);
        Arc::make_mut(&mut systems.0).permute_participants(
            RelationId(0),
            &[
                ParticipantPosition(2),
                ParticipantPosition(0),
                ParticipantPosition(1),
            ],
        );
        systems
    }

    #[rstest]
    fn test_multicenter_bonds_reframe(unsorted_bond: MulticenterBonds) {
        assert_eq!(
            unsorted_bond
                .atoms(MulticenterBondId(0))
                .collect::<Vec<_>>(),
            vec![AtomId(7), AtomId(1), AtomId(4)],
        );

        let reframed = unsorted_bond.reframe().expect("the form is satisfiable");

        assert_eq!(
            reframed.atoms(MulticenterBondId(0)).collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4), AtomId(7)],
        );
        assert_eq!(
            reframed.attributes(MulticenterBondId(0)),
            &MulticenterBondForm::from_electrons(vec![20, 30, 10]),
        );
    }

    #[rstest]
    fn test_multicenter_bonds_reframe_identity(unsorted_bond: MulticenterBonds) {
        let once = unsorted_bond.reframe().expect("the form is satisfiable");
        let twice = once.clone().reframe().expect("the form is satisfiable");
        assert_eq!(twice, once);
    }

    #[rstest]
    fn test_multicenter_bonds_tracked_reframe(unsorted_bond: MulticenterBonds) {
        let (reframed, actions) = unsorted_bond
            .clone()
            .tracked_reframe()
            .expect("the form is satisfiable");

        let action = actions
            .action(MulticenterBondId(0))
            .expect("the dense action covers the bond");
        assert_eq!(action.image(), [1, 2, 0]);
        assert_eq!(unsorted_bond.reframe_by(&actions), Some(reframed));
    }

    #[rstest]
    fn test_reframe_multicenter_bonds_with(unsorted_bond: MulticenterBonds) {
        let mut visited = None;
        let reframed = reframe_multicenter_bonds_with(unsorted_bond.clone(), |id, action| {
            visited = Some((id, action.clone()));
        })
        .expect("the form is satisfiable");

        assert_eq!(
            visited,
            Some((
                MulticenterBondId(0),
                DynPermutation::try_from(vec![1, 2, 0]).expect("the expected action is valid"),
            )),
        );
        assert_eq!(unsorted_bond.reframe(), Ok(reframed));
    }

    #[rstest]
    fn test_multicenter_bonds_normalize() {
        let bonds = MulticenterBonds::new(vec![(
            vec![AtomId(1), AtomId(2)],
            MulticenterBondForm::default().with_charge(NumForm::lit_set([0])),
        )]);

        let normalized = bonds.normalize().expect("the form is satisfiable");

        assert_eq!(
            normalized.attributes(MulticenterBondId(0)),
            &MulticenterBondForm::default().with_charge(0),
        );
    }

    #[rstest]
    fn test_multicenter_bonds_framed_eq(unsorted_bond: MulticenterBonds) {
        let selected = MulticenterBonds::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            MulticenterBondForm::from_electrons(vec![20, 30, 10]),
        )]);
        assert!(unsorted_bond.framed_eq(&selected));
        assert!(!unsorted_bond.eq(&selected));

        let different = MulticenterBonds::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            MulticenterBondForm::from_electrons(vec![10, 20, 30]),
        )]);
        assert!(!unsorted_bond.framed_eq(&different));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::literal(MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3])),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: MulticenterBondConstraintsForm::new() })]
    fn test_multicenter_bond_form_new(
        #[case] actual: MulticenterBondForm,
        #[case] expected: MulticenterBondForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::literal(MulticenterBondForm::from_electrons(vec![1, 1, 1]),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: MulticenterBondConstraintsForm::new() })]
    fn test_multicenter_bond_form_from_electrons(
        #[case] actual: MulticenterBondForm,
        #[case] expected: MulticenterBondForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_charge(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(-1),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Lit(-1), unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: MulticenterBondConstraintsForm::new() })]
    #[case::with_unpaired_electrons(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((0_u8, 1_u8)),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::closed_shell(),
            constraints: MulticenterBondConstraintsForm::new() })]
    #[case::with_constraint(
        MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintForm::electron_count(2)),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(2)) })]
    #[case::with_constraints_extends(
        MulticenterBondForm::from_electrons(vec![1, 1, 1])
            .with_constraints([MulticenterBondConstraintForm::electron_count(2)]),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(2)) })]
    #[case::with_constraint_replaces_same_kind(
        MulticenterBondForm::from_electrons(vec![1, 1, 1])
            .with_constraint(MulticenterBondConstraintForm::electron_count(2))
            .with_constraint(MulticenterBondConstraintForm::electron_count(4)),
        MulticenterBondForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(4)) })]
    fn test_multicenter_bond_form_with_methods(
        #[case] actual: MulticenterBondForm,
        #[case] expected: MulticenterBondForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_ground_electrons(
        MulticenterBondForm::from_electrons(vec![1; 3]).into_concrete(),
        MulticenterBondForm {
            electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraintsForm::new(),
        },
    )]
    #[case::preserves_set_charge(
        MulticenterBondForm::from_electrons(vec![1; 3]).with_charge(1_i64).into_concrete(),
        MulticenterBondForm {
            electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Lit(1),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraintsForm::new(),
        },
    )]
    #[case::preserves_constraints(
        MulticenterBondForm::from_electrons(vec![1; 3])
            .with_constraint(MulticenterBondConstraintForm::electron_count(3))
            .into_concrete(),
        MulticenterBondForm {
            electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraintsForm::from(
                MulticenterBondConstraintForm::electron_count(3),
            ),
        },
    )]
    fn test_multicenter_bond_form_into_concrete(
        #[case] actual: MulticenterBondForm,
        #[case] expected: MulticenterBondForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electrons(MulticenterBondForm::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])), ..Default::default() }, MulticenterBondForm::from_electrons(vec![2, 2, 2]))]
    #[case::electrons_undetermined(MulticenterBondForm::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { electrons: Some(ElectronCountsForm::Undetermined), ..Default::default() }, MulticenterBondForm::default())]
    #[case::charge(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(0_i64), MulticenterBondUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(-1_i64))]
    #[case::charge_undetermined(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(-1_i64), MulticenterBondUpdate { charge: Some(NumForm::Undetermined), ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]))]
    #[case::unpaired_electrons_count(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 3_u8)), MulticenterBondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Lit(0)), multiplicity: None }, ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((0_u8, 3_u8)))]
    #[case::unpaired_electrons_multiplicity(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 3_u8)), MulticenterBondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 1_u8)))]
    #[case::constraint_set(MulticenterBondForm::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(6_i64)), ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintForm::electron_count(6_i64)))]
    #[case::constraint_replace(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintForm::electron_count(6_i64)), MulticenterBondUpdate { constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(4_i64)), ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintForm::electron_count(4_i64)))]
    #[case::constraint_remove(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_constraint(MulticenterBondConstraintForm::electron_count(6_i64)), MulticenterBondUpdate { constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(NumForm::Undetermined)), ..Default::default() }, MulticenterBondForm::from_electrons(vec![1, 1, 1]))]
    fn test_multicenter_bond_form_update(
        #[case] bond: MulticenterBondForm,
        #[case] update: MulticenterBondUpdate,
        #[case] expected: MulticenterBondForm,
    ) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(-1_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(MulticenterBondConstraintForm::electron_count(6_i64)))]
    fn test_multicenter_bond_form_update_identity(#[case] bond: MulticenterBondForm) {
        assert_eq!(bond.update(&MulticenterBondUpdate::default()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraints(
        MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(0_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(MulticenterBondConstraintForm::electron_count(6_i64)),
        MulticenterBondForm::from_electrons(vec![2, 2, 2]).with_unpaired_electrons((2_u8, 1_u8)),
        MulticenterBondUpdate {
            electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])),
            charge: Some(NumForm::Undetermined),
            unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) },
            constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(NumForm::Undetermined)),
        },
    )]
    fn test_multicenter_bond_form_difference_to(
        #[case] bond: MulticenterBondForm,
        #[case] other: MulticenterBondForm,
        #[case] expected: MulticenterBondUpdate,
    ) {
        assert_eq!(bond.difference_to(&other), expected);
    }

    #[rstest]
    #[case::normalized(
        MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(1_i64),
        MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(NumForm::lit_set([1])),
    )]
    fn test_multicenter_bond_form_difference_to_identity(
        #[case] bond: MulticenterBondForm,
        #[case] other: MulticenterBondForm,
    ) {
        assert_eq!(bond.difference_to(&other), MulticenterBondUpdate::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::positioned(
        MulticenterBondForm::from_electrons(vec![10, 20, 30]).with_charge(-1),
        vec![2, 0, 1],
        Some(MulticenterBondForm::from_electrons(vec![30, 10, 20]).with_charge(-1)),
    )]
    #[case::positioned_degree(
        MulticenterBondForm::from_electrons(vec![10, 20]),
        vec![2, 0, 1],
        None,
    )]
    #[case::dimensionless(
        MulticenterBondForm::default(),
        vec![3, 1, 0, 2],
        Some(MulticenterBondForm::default()),
    )]
    #[case::frame_invariant_constraint(
        MulticenterBondForm::default()
            .with_constraint(MulticenterBondConstraintForm::electron_count(2)),
        vec![1, 0],
        Some(
            MulticenterBondForm::default()
                .with_constraint(MulticenterBondConstraintForm::electron_count(2)),
        ),
    )]
    fn test_multicenter_bond_form_reframe_by(
        #[case] input: MulticenterBondForm,
        #[case] image: Vec<usize>,
        #[case] expected: Option<MulticenterBondForm>,
    ) {
        let action = DynPermutation::try_from(image).expect("case is a permutation");
        assert_eq!(input.reframe_by(&action), expected);
    }

    #[rstest]
    #[case::three_members(
        MulticenterBondForm::from_electrons(vec![10, 20, 30]).with_charge(-1),
        vec![
            ParticipantPosition(2),
            ParticipantPosition(0),
            ParticipantPosition(1),
        ],
        MulticenterBondForm::from_electrons(vec![30, 10, 20]).with_charge(-1),
    )]
    fn test_multicenter_bond_form_permute(
        #[case] mut input: MulticenterBondForm,
        #[case] order: Vec<ParticipantPosition>,
        #[case] expected: MulticenterBondForm,
    ) {
        input.permute(&order);
        assert_eq!(input, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(MulticenterBondForm::default(), false)]
    #[case::charge_only(MulticenterBondForm::new(ElectronCountsForm::Undetermined).with_charge(0), false)]
    #[case::ground_no_atoms(MulticenterBondForm::new(ElectronCountsForm::Lit(Vec::new())).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::all_ground_three(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        true,
    )]
    #[case::ground_with_constraint(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3]))
            .with_charge(0).with_unpaired_electrons((0, 1))
            .with_constraint(MulticenterBondConstraintForm::electron_count(3)),
        true,
    )]
    fn test_multicenter_bond_form_is_ground(
        #[case] form: MulticenterBondForm,
        #[case] expected: bool,
    ) {
        assert_eq!(form.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_charge(
        MulticenterBondForm::default().with_charge(NumForm::lit_set([0])),
        Ok(MulticenterBondForm::default().with_charge(0)),
    )]
    #[case::charge_empty_litset_contradiction(
        MulticenterBondForm::default().with_charge(NumForm::lit_set(Vec::<i64>::new())),
        Err(Contradiction),
    )]
    fn test_multicenter_bond_form_normalize(
        #[case] input: MulticenterBondForm,
        #[case] expected: Result<MulticenterBondForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_default(MulticenterBondForm::default(), MulticenterBondForm::default(), true)]
    #[case::default_matches_ground(
        MulticenterBondForm::default(),
        MulticenterBondForm::new(ElectronCountsForm::Lit(Vec::new())).with_charge(0).with_unpaired_electrons((0, 1)),
        true,
    )]
    #[case::exact(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        true,
    )]
    #[case::electrons_length_mismatch(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 2])),
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        false,
    )]
    #[case::electrons_value_mismatch(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![2; 3])),
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![1; 3])).with_charge(0).with_unpaired_electrons((0, 1)),
        false,
    )]
    #[case::charge_mismatch(
        MulticenterBondForm::new(ElectronCountsForm::Undetermined).with_charge(1),
        MulticenterBondForm::new(ElectronCountsForm::Undetermined).with_charge(0),
        false,
    )]
    #[case::unpaired_electrons_mismatch(
        MulticenterBondForm::new(ElectronCountsForm::Undetermined).with_unpaired_electrons((2_u8, 3_u8)),
        MulticenterBondForm::new(ElectronCountsForm::Undetermined).with_unpaired_electrons((0_u8, 1_u8)),
        false,
    )]
    #[case::constraint_required_present(
        MulticenterBondForm::new(ElectronCountsForm::Undetermined)
            .with_constraint(MulticenterBondConstraintForm::electron_count(3)),
        MulticenterBondForm::new(ElectronCountsForm::Undetermined)
            .with_constraint(MulticenterBondConstraintForm::electron_count(3)),
        true,
    )]
    #[case::constraint_required_absent(
        MulticenterBondForm::new(ElectronCountsForm::Undetermined)
            .with_constraint(MulticenterBondConstraintForm::electron_count(3)),
        MulticenterBondForm::new(ElectronCountsForm::Undetermined),
        false,
    )]
    fn test_multicenter_bond_form_matches(
        #[case] pattern: MulticenterBondForm,
        #[case] target: MulticenterBondForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(
        MulticenterBondForm::default(),
        MulticenterBondForm::default(),
        Some(MulticenterBondForm::default())
    )]
    #[case::electrons_length_mismatch(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![2; 3])),
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![2; 4])),
        None,
    )]
    #[case::narrows_electrons(
        MulticenterBondForm::new(ElectronCountsForm::Undetermined),
        MulticenterBondForm::from_electrons(vec![1, 2]),
        Some(MulticenterBondForm::from_electrons(vec![1, 2])),
    )]
    fn test_multicenter_bond_form_meet(
        #[case] a: MulticenterBondForm,
        #[case] b: MulticenterBondForm,
        #[case] expected: Option<MulticenterBondForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::electrons_length_mismatch_widens_to_default(
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![2; 3])),
        MulticenterBondForm::new(ElectronCountsForm::Lit(vec![2; 4])),
        MulticenterBondForm::default(),
    )]
    fn test_multicenter_bond_form_join(
        #[case] a: MulticenterBondForm,
        #[case] b: MulticenterBondForm,
        #[case] expected: MulticenterBondForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }
}
