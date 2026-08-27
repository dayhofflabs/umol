//! Aromatic systems: the molecule's collection and one system's attribute form.

use std::sync::Arc;

use umol_graph_core::{
    NodeId, ParticipantPosition, RelationData, RelationId, Remapping, Unordered, VarRelationSet,
};
use umol_graph_ir_macros::{Lattice, Normalize};

use super::constraint::{AromaticSystemConstraintForm, AromaticSystemConstraintsForm};
use super::delta::EntitySpan;
use super::electrons::ElectronCountsForm;
use super::error::Contradiction;
use super::id::{AromaticSystemId, AtomId};
use super::num::NumForm;
use super::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
use super::traits::{Equiv, Lattice, Normalize, Reframe};

/// The molecule's aromatic systems.
///
/// The atoms bear the participant frame: the per-member electron counts of
/// [`AromaticSystemForm`] are read against it, position by position.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AromaticSystems(Arc<VarRelationSet<NodeId, Unordered, AromaticSystemForm>>);

impl AromaticSystems {
    pub fn new(entries: Vec<(Vec<AtomId>, AromaticSystemForm)>) -> Self {
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

    pub fn contains(&self, id: AromaticSystemId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = AromaticSystemId> {
        self.0.ids().map(AromaticSystemId::from)
    }

    /// The atoms of `id`, in their stored frame.
    pub fn atoms(&self, id: AromaticSystemId) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.0
            .participants(RelationId::from(id))
            .iter()
            .map(|&atom| AtomId::from(atom))
    }

    pub fn attributes(&self, id: AromaticSystemId) -> &AromaticSystemForm {
        self.0.data(RelationId::from(id))
    }

    pub fn attributes_mut(&mut self, id: AromaticSystemId) -> &mut AromaticSystemForm {
        Arc::make_mut(&mut self.0).data_mut(RelationId::from(id))
    }

    /// Ids of the systems `atom` belongs to. Systems are atom-disjoint, so there is at most one.
    pub fn incident_ids(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = AromaticSystemId> + '_ {
        self.0
            .incident(NodeId::from(atom))
            .iter()
            .map(|&id| AromaticSystemId::from(id))
    }

    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.0.has_incident(NodeId::from(atom))
    }

    pub fn into_entries(self) -> Vec<(Vec<AtomId>, AromaticSystemForm)> {
        Arc::try_unwrap(self.0)
            .unwrap_or_else(|shared| (*shared).clone())
            .into_entries()
            .into_iter()
            .map(|(atoms, attributes)| (atoms.into_iter().map(AtomId::from).collect(), attributes))
            .collect()
    }

    /// The atoms of `id` as graph nodes, for graph-core interop that is not yet typed in graph-IR
    /// ids. The public accessor is [`Self::atoms`].
    pub(crate) fn atom_nodes(&self, id: AromaticSystemId) -> &[NodeId] {
        self.0.participants(RelationId::from(id))
    }

    pub(crate) fn attributes_iter_mut(
        &mut self,
    ) -> impl ExactSizeIterator<Item = &mut AromaticSystemForm> {
        Arc::make_mut(&mut self.0)
            .iter_mut()
            .map(|(_, _, attributes)| attributes)
    }

    pub(crate) fn remap(&self, remapping: &Remapping) -> Self {
        Self(Arc::new(self.0.remap(remapping)))
    }

    pub(crate) fn into_arc(self) -> Arc<VarRelationSet<NodeId, Unordered, AromaticSystemForm>> {
        self.0
    }

    /// Glue `right`, relabelled into this molecule's id space, onto `self`: coinciding systems meet,
    /// non-coinciding systems are carried. `None` when a coincident meet is bottom.
    pub(crate) fn glue(&self, right: &Self, remapping: &Remapping) -> Option<Self> {
        self.0
            .pushout(
                &right.remap(remapping).0,
                // Aromatic systems anchor on their atoms: the node index.
                |set, atoms| atoms.first().and_then(|&node| set.coincident(node, atoms)),
                |(left_atoms, left), (right_atoms, right)| {
                    let left_atoms: Vec<AtomId> =
                        left_atoms.iter().map(|&atom| AtomId::from(atom)).collect();
                    let right_atoms: Vec<AtomId> =
                        right_atoms.iter().map(|&atom| AtomId::from(atom)).collect();
                    right
                        .clone()
                        .reframe_to(&right_atoms, &left_atoms)?
                        .meet(left)
                },
            )
            .map(|merged| Self(Arc::new(merged.object)))
    }

    /// Id of the system coinciding with `atoms` — the one whose atoms equal them as a multiset.
    ///
    /// The identity question, distinct from lookup: an aromatic system's uniqueness key is any
    /// member atom, which names it from a part; this names it from the whole.
    pub fn coincident_id(&self, atoms: &[AtomId]) -> Option<AromaticSystemId> {
        // Aromatic systems anchor on their atoms, so the node index is the one to scan.
        let query: Vec<NodeId> = atoms.iter().map(|&atom| NodeId::from(atom)).collect();
        let anchor = *query.first()?;
        self.0
            .coincident(anchor, &query)
            .map(AromaticSystemId::from)
    }
}

impl From<VarRelationSet<NodeId, Unordered, AromaticSystemForm>> for AromaticSystems {
    fn from(set: VarRelationSet<NodeId, Unordered, AromaticSystemForm>) -> Self {
        Self(Arc::new(set))
    }
}

impl From<Arc<VarRelationSet<NodeId, Unordered, AromaticSystemForm>>> for AromaticSystems {
    fn from(set: Arc<VarRelationSet<NodeId, Unordered, AromaticSystemForm>>) -> Self {
        Self(set)
    }
}

impl Reframe for AromaticSystems {
    type Action = (AromaticSystemId, Vec<ParticipantPosition>);

    /// Reduce every entry, then present each in its selected frame, returning the frame action selected for each system.
    ///
    /// The action is the position order taking the stored frame to the selected one, so
    /// `reframe_by(reduce(x), action)` reproduces the reframed value.
    fn reframe_with_action(&self) -> Result<(Self, Vec<Self::Action>), Contradiction> {
        let mut reframed = (*self.0).clone();
        let mut actions = Vec::with_capacity(reframed.count());
        for id in reframed.ids().collect::<Vec<_>>() {
            let stored: Vec<AtomId> = reframed
                .participants(id)
                .iter()
                .map(|&atom| AtomId::from(atom))
                .collect();
            let mut order: Vec<ParticipantPosition> = (0..stored.len())
                .map(|position| ParticipantPosition(position as u32))
                .collect();
            order.sort_by_key(|position| stored[position.index()]);
            let selected: Vec<AtomId> = order
                .iter()
                .map(|position| stored[position.index()])
                .collect();

            let attributes = reframed.data(id).clone().normalize()?;
            *reframed.data_mut(id) = attributes
                .reframe_to(&stored, &selected)
                .ok_or(Contradiction)?;
            reframed.permute_with(id, &order);
            actions.push((AromaticSystemId::from(id), order));
        }
        Ok((Self(Arc::new(reframed)), actions))
    }
}

/// The reaction span's aromatic systems, one [`EntitySpan`] per entity against a single participant frame.
///
/// The `Molecule` peer is [`AromaticSystems`]. The surface is deliberately duplicated rather than shared
/// through a payload parameter: a type parameter on the molecule-level families would complicate
/// the primary carrier to serve this one.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AromaticSystemSpans(VarRelationSet<NodeId, Unordered, EntitySpan<AromaticSystemForm>>);

impl AromaticSystemSpans {
    pub fn into_entries(self) -> Vec<(Vec<AtomId>, EntitySpan<AromaticSystemForm>)> {
        self.0
            .into_entries()
            .into_iter()
            .map(|(atoms, span)| (atoms.into_iter().map(AtomId::from).collect(), span))
            .collect()
    }

    pub fn new(entries: Vec<(Vec<AtomId>, EntitySpan<AromaticSystemForm>)>) -> Self {
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

    pub fn contains(&self, id: AromaticSystemId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = AromaticSystemId> {
        self.0.ids().map(AromaticSystemId::from)
    }

    /// The atoms of `id`, in their stored frame.
    pub fn atoms(&self, id: AromaticSystemId) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.0
            .participants(RelationId::from(id))
            .iter()
            .map(|&atom| AtomId::from(atom))
    }

    pub fn attributes(&self, id: AromaticSystemId) -> &EntitySpan<AromaticSystemForm> {
        self.0.data(RelationId::from(id))
    }

    pub(crate) fn remap(&self, remapping: &Remapping) -> Self {
        Self(self.0.remap(remapping))
    }
}

impl Reframe for AromaticSystemSpans {
    type Action = (AromaticSystemId, Vec<ParticipantPosition>);

    /// Selection never consults the payload here — the atoms sort — so a `Modified` span needs no
    /// arbitration between its sides. One action carries every side through [`EntitySpan::try_map`],
    /// and the span declines whole if any side declines.
    fn reframe_with_action(&self) -> Result<(Self, Vec<Self::Action>), Contradiction> {
        let mut reframed = self.0.clone();
        let mut actions = Vec::with_capacity(reframed.count());
        for id in reframed.ids().collect::<Vec<_>>() {
            let stored: Vec<AtomId> = reframed
                .participants(id)
                .iter()
                .map(|&atom| AtomId::from(atom))
                .collect();
            let mut order: Vec<ParticipantPosition> = (0..stored.len())
                .map(|position| ParticipantPosition(position as u32))
                .collect();
            order.sort_by_key(|position| stored[position.index()]);
            let selected: Vec<AtomId> = order
                .iter()
                .map(|position| stored[position.index()])
                .collect();

            let span = reframed.data(id).clone();
            let reduced = span
                .try_map(|form| form.normalize().ok())
                .ok_or(Contradiction)?;
            *reframed.data_mut(id) = reduced
                .try_map(|form| form.reframe_to(&stored, &selected))
                .ok_or(Contradiction)?;
            reframed.permute_with(id, &order);
            actions.push((AromaticSystemId::from(id), order));
        }
        Ok((Self(reframed), actions))
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Normalize, Lattice)]
pub struct AromaticSystemForm {
    pub electrons: ElectronCountsForm,
    pub charge: NumForm,
    pub unpaired_electrons: UnpairedElectronsForm,
    pub constraints: AromaticSystemConstraintsForm,
}

/// Attribute update for an aromatic system. Ordinary fields are optional,
/// unpaired-electron components are updated independently, and undetermined constraints remove
/// their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AromaticSystemUpdate {
    pub electrons: Option<ElectronCountsForm>,
    pub charge: Option<NumForm>,
    pub unpaired_electrons: UnpairedElectronsUpdate,
    pub constraints: AromaticSystemConstraintsForm,
}

impl From<&str> for AromaticSystemForm {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid aromatic system string")
    }
}

impl RelationData for AromaticSystemForm {
    /// The per-member electron counts are positional, so they follow a participant reorder.
    fn on_permutation(&mut self, order: &[ParticipantPosition]) {
        self.electrons.permute(order);
    }

    fn is_permutation_invariant(&self) -> bool {
        self.electrons.is_undetermined()
    }
}

impl AromaticSystemForm {
    /// Concrete: every inherent field is ground; the constraint channel does
    /// not bear on concreteness.
    pub fn is_concrete(&self) -> bool {
        let AromaticSystemForm {
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
    /// kind (last-wins per `AromaticSystemConstraintsForm::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<AromaticSystemConstraintForm>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `AromaticSystemConstraintsForm::set`).
    /// Does not clear existing constraints; use `system.constraints.clear()`
    /// or direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<AromaticSystemConstraintForm>,
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
    pub fn update(&self, update: &AromaticSystemUpdate) -> AromaticSystemForm {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        AromaticSystemForm {
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
    pub fn difference_to(&self, other: &Self) -> AromaticSystemUpdate {
        let mut constraints = AromaticSystemConstraintsForm::new();
        for new in other.constraints.iter() {
            if self
                .constraints
                .get(new.key())
                .is_none_or(|old| !old.equiv(new))
            {
                constraints.set(new.clone());
            }
        }
        for old in self.constraints.iter() {
            if other.constraints.get(old.key()).is_none() {
                constraints.set(old.as_undetermined());
            }
        }
        AromaticSystemUpdate {
            electrons: (!self.electrons.equiv(&other.electrons)).then(|| other.electrons.clone()),
            charge: (!self.charge.equiv(&other.charge)).then(|| other.charge.clone()),
            unpaired_electrons: self
                .unpaired_electrons
                .difference_to(&other.unpaired_electrons),
            constraints,
        }
    }

    /// Restate the form in the `to` participant frame, given it is stated in the `from` frame.
    /// Only `electrons` is position-indexed; charge, unpaired electrons, and constraints are
    /// positionless and carry unchanged. `None` when `electrons` declines the frame change.
    ///
    /// Destructured exhaustively on purpose: a new position-indexed field must fail to compile
    /// here rather than be silently left in the old frame.
    pub fn reframe_to(self, from: &[AtomId], to: &[AtomId]) -> Option<Self> {
        let Self {
            electrons,
            charge,
            unpaired_electrons,
            constraints,
        } = self;
        Some(Self {
            electrons: electrons.reframe_to(from, to)?,
            charge,
            unpaired_electrons,
            constraints,
        })
    }

    /// Reorder the positional `electrons` by `order`, tracking a participant
    /// reordering; charge / unpaired electrons / constraints are positionless and unchanged.
    pub fn permute(&mut self, order: &[ParticipantPosition]) {
        self.electrons.permute(order);
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ir::error::Contradiction;
    use crate::ir::traits::Normalize;

    /// A coincidence whose two sides hold the same atoms in different frames, with nonuniform
    /// electron counts stating the same per-atom fact. The glue must carry the right vector into
    /// the left frame before meeting; meeting the two vectors position-by-position without that
    /// compares different atoms' counts and yields bottom.
    ///
    /// Both sides state 1 -> 20, 4 -> 30, 7 -> 10. The left frame is permuted apart from sorted
    /// order, as the frame-preserving storage S5b introduces would leave it; the right stays sorted
    /// because `remap` reconstructs it through `new`, which still sorts an `Unordered` factor.
    #[rstest]
    fn test_aromatic_systems_glue_differing_frames() {
        let mut left = AromaticSystems::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            AromaticSystemForm::from_electrons(vec![20, 30, 10]),
        )]);
        Arc::make_mut(&mut left.0).permute_with(
            RelationId(0),
            &[
                ParticipantPosition(2),
                ParticipantPosition(0),
                ParticipantPosition(1),
            ],
        );
        // The permute moves participants and leaves the payload, so the left frame now reads
        // 7 -> 20, 1 -> 30, 4 -> 10 positionally — which is a different fact from the right's.
        Arc::make_mut(&mut left.0)
            .iter_mut()
            .for_each(|(_, _, attributes)| {
                *attributes = AromaticSystemForm::from_electrons(vec![10, 20, 30])
            });
        assert_eq!(
            left.atoms(AromaticSystemId(0)).collect::<Vec<_>>(),
            vec![AtomId(7), AtomId(1), AtomId(4)],
            "left states 7 -> 10, 1 -> 20, 4 -> 30",
        );

        let right = AromaticSystems::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            AromaticSystemForm::from_electrons(vec![20, 30, 10]),
        )]);
        assert_eq!(
            right.atoms(AromaticSystemId(0)).collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4), AtomId(7)],
            "right states the same fact in sorted order",
        );

        let glued = left
            .glue(
                &right,
                &Remapping::new((0..8).map(NodeId).collect(), vec![]),
            )
            .expect("the sides agree once the right vector is carried into the left frame");

        assert_eq!(glued.count(), 1);
        assert_eq!(
            glued.attributes(AromaticSystemId(0)),
            &AromaticSystemForm::from_electrons(vec![20, 30, 10]),
        );
    }

    /// Both sides of a `Modified` span are read against one participant list, so one action carries
    /// both. Selection never consults the payload here, so the two sides cannot disagree about it.
    #[rstest]
    fn test_aromatic_system_spans_reframe() {
        let mut spans = AromaticSystemSpans::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            EntitySpan::Modified {
                lhs: AromaticSystemForm::from_electrons(vec![10, 20, 30]),
                rhs: AromaticSystemForm::from_electrons(vec![11, 21, 31]),
            },
        )]);
        spans.0.permute_with(
            RelationId(0),
            &[
                ParticipantPosition(2),
                ParticipantPosition(0),
                ParticipantPosition(1),
            ],
        );

        let (reframed, actions) = spans
            .reframe_with_action()
            .expect("the forms are satisfiable");

        assert_eq!(
            reframed.atoms(AromaticSystemId(0)).collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4), AtomId(7)],
        );
        assert_eq!(
            reframed.attributes(AromaticSystemId(0)),
            &EntitySpan::Modified {
                lhs: AromaticSystemForm::from_electrons(vec![20, 30, 10]),
                rhs: AromaticSystemForm::from_electrons(vec![21, 31, 11]),
            },
        );
        assert_eq!(
            actions,
            vec![(
                AromaticSystemId(0),
                vec![
                    ParticipantPosition(1),
                    ParticipantPosition(2),
                    ParticipantPosition(0),
                ],
            )],
        );
    }

    #[rstest]
    fn test_aromatic_system_spans_reframe_identity() {
        let spans = AromaticSystemSpans::new(vec![(
            vec![AtomId(1), AtomId(4)],
            EntitySpan::Modified {
                lhs: AromaticSystemForm::from_electrons(vec![10, 20]),
                rhs: AromaticSystemForm::from_electrons(vec![11, 21]),
            },
        )]);
        let once = spans.reframe().expect("the forms are satisfiable");
        let twice = once.reframe().expect("the forms are satisfiable");
        assert_eq!(twice, once);
    }

    /// A side that declines the frame change takes the whole span with it: one action serves both,
    /// so there is no partial result to keep. Here the rhs electron vector disagrees in length with
    /// the participant frame.
    #[rstest]
    fn test_aromatic_system_spans_reframe_error() {
        let spans = AromaticSystemSpans::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            EntitySpan::Modified {
                lhs: AromaticSystemForm::from_electrons(vec![10, 20, 30]),
                rhs: AromaticSystemForm::from_electrons(vec![11, 21]),
            },
        )]);
        assert_eq!(spans.reframe(), Err(Contradiction));
    }

    #[rstest]
    fn test_aromatic_system_spans_framed_eq() {
        let mut unsorted = AromaticSystemSpans::new(vec![(
            vec![AtomId(1), AtomId(4)],
            EntitySpan::Modified {
                lhs: AromaticSystemForm::from_electrons(vec![10, 20]),
                rhs: AromaticSystemForm::from_electrons(vec![11, 21]),
            },
        )]);
        unsorted.0.permute_with(
            RelationId(0),
            &[ParticipantPosition(1), ParticipantPosition(0)],
        );
        // The permute leaves the payload where it was, so the stored frame now reads
        // atom 4 -> 10, atom 1 -> 20; the selected presentation states the same fact sorted.
        let selected = AromaticSystemSpans::new(vec![(
            vec![AtomId(1), AtomId(4)],
            EntitySpan::Modified {
                lhs: AromaticSystemForm::from_electrons(vec![20, 10]),
                rhs: AromaticSystemForm::from_electrons(vec![21, 11]),
            },
        )]);

        assert!(unsorted.framed_eq(&selected));
        assert!(unsorted != selected);
    }

    /// Storage sorts an `Unordered` factor on construction, so the stored frame is permuted first
    /// to model the frame-preserving storage S5 introduces.
    #[fixture]
    fn unsorted_system() -> AromaticSystems {
        let mut systems = AromaticSystems::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            AromaticSystemForm::from_electrons(vec![10, 20, 30]),
        )]);
        Arc::make_mut(&mut systems.0).permute_with(
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
    fn test_aromatic_systems_reframe(unsorted_system: AromaticSystems) {
        assert_eq!(
            unsorted_system
                .atoms(AromaticSystemId(0))
                .collect::<Vec<_>>(),
            vec![AtomId(7), AtomId(1), AtomId(4)],
        );

        let reframed = unsorted_system.reframe().expect("the form is satisfiable");

        assert_eq!(
            reframed.atoms(AromaticSystemId(0)).collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4), AtomId(7)],
        );
        assert_eq!(
            reframed.attributes(AromaticSystemId(0)),
            &AromaticSystemForm::from_electrons(vec![20, 30, 10]),
        );
    }

    #[rstest]
    fn test_aromatic_systems_reframe_identity(unsorted_system: AromaticSystems) {
        let once = unsorted_system.reframe().expect("the form is satisfiable");
        let twice = once.reframe().expect("the form is satisfiable");
        assert_eq!(twice, once);
    }

    #[rstest]
    fn test_aromatic_systems_reframe_with_action(unsorted_system: AromaticSystems) {
        let (reframed, actions) = unsorted_system
            .reframe_with_action()
            .expect("the form is satisfiable");

        assert_eq!(
            actions,
            vec![(
                AromaticSystemId(0),
                vec![
                    ParticipantPosition(1),
                    ParticipantPosition(2),
                    ParticipantPosition(0),
                ],
            )],
        );

        let (_, order) = &actions[0];
        let stored: Vec<AtomId> = unsorted_system.atoms(AromaticSystemId(0)).collect();
        let selected: Vec<AtomId> = order.iter().map(|p| stored[p.index()]).collect();
        assert_eq!(
            unsorted_system
                .attributes(AromaticSystemId(0))
                .clone()
                .reframe_to(&stored, &selected),
            Some(reframed.attributes(AromaticSystemId(0)).clone()),
        );
    }

    #[rstest]
    fn test_aromatic_systems_framed_eq(unsorted_system: AromaticSystems) {
        let selected = AromaticSystems::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            AromaticSystemForm::from_electrons(vec![20, 30, 10]),
        )]);
        assert!(unsorted_system.framed_eq(&selected));
        assert!(!unsorted_system.eq(&selected));

        let different = AromaticSystems::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            AromaticSystemForm::from_electrons(vec![10, 20, 30]),
        )]);
        assert!(!unsorted_system.framed_eq(&different));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::literal(AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 3])),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: AromaticSystemConstraintsForm::new() })]
    fn test_aromatic_system_form_new(
        #[case] actual: AromaticSystemForm,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::literal(AromaticSystemForm::from_electrons(vec![1, 1, 1]),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1; 3]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: AromaticSystemConstraintsForm::new() })]
    fn test_aromatic_system_form_from_electrons(
        #[case] actual: AromaticSystemForm,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_charge(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(-1),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Lit(-1), unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: AromaticSystemConstraintsForm::new() })]
    #[case::with_unpaired_electrons(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((0_u8, 1_u8)),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::closed_shell(),
            constraints: AromaticSystemConstraintsForm::new() })]
    #[case::with_constraint(
        AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintForm::electron_count(6)),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]),
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)) })]
    #[case::with_constraint_replaces_same_kind(
        AromaticSystemForm::default()
            .with_constraint(AromaticSystemConstraintForm::electron_count(2))
            .with_constraint(AromaticSystemConstraintForm::electron_count(6)),
        AromaticSystemForm { electrons: ElectronCountsForm::Undetermined,
            charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)) })]
    fn test_aromatic_system_form_with_methods(
        #[case] actual: AromaticSystemForm,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_ground_electrons(AromaticSystemForm::from_electrons(vec![1; 6]).into_concrete(),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1; 6]), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraintsForm::new() })]
    #[case::preserves_set_charge(AromaticSystemForm::from_electrons(vec![1; 6]).with_charge(1_i64).into_concrete(),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1; 6]), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraintsForm::new() })]
    #[case::preserves_constraints(AromaticSystemForm::from_electrons(vec![1; 6]).with_constraint(AromaticSystemConstraintForm::electron_count(6)).into_concrete(),
        AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1; 6]), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6)) })]
    fn test_aromatic_system_form_into_concrete(
        #[case] actual: AromaticSystemForm,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electrons(AromaticSystemForm::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])), ..Default::default() }, AromaticSystemForm::from_electrons(vec![2, 2, 2]))]
    #[case::electrons_undetermined(AromaticSystemForm::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { electrons: Some(ElectronCountsForm::Undetermined), ..Default::default() }, AromaticSystemForm::default())]
    #[case::charge(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(0_i64), AromaticSystemUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(-1_i64))]
    #[case::charge_undetermined(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(-1_i64), AromaticSystemUpdate { charge: Some(NumForm::Undetermined), ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]))]
    #[case::unpaired_electrons_count(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 3_u8)), AromaticSystemUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Lit(0)), multiplicity: None }, ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((0_u8, 3_u8)))]
    #[case::unpaired_electrons_multiplicity(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 3_u8)), AromaticSystemUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 1_u8)))]
    #[case::constraint_set(AromaticSystemForm::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6_i64)), ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintForm::electron_count(6_i64)))]
    #[case::constraint_replace(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintForm::electron_count(6_i64)), AromaticSystemUpdate { constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(4_i64)), ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintForm::electron_count(4_i64)))]
    #[case::constraint_remove(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintForm::electron_count(6_i64)), AromaticSystemUpdate { constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)), ..Default::default() }, AromaticSystemForm::from_electrons(vec![1, 1, 1]))]
    fn test_aromatic_system_form_update(
        #[case] system: AromaticSystemForm,
        #[case] update: AromaticSystemUpdate,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(system.update(&update), expected);
    }

    #[rstest]
    #[case::empty(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(-1_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(AromaticSystemConstraintForm::electron_count(6_i64)))]
    fn test_aromatic_system_form_update_identity(#[case] system: AromaticSystemForm) {
        assert_eq!(system.update(&AromaticSystemUpdate::default()), system);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraints(
        AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(0_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(AromaticSystemConstraintForm::electron_count(6_i64)),
        AromaticSystemForm::from_electrons(vec![2, 2, 2]).with_unpaired_electrons((2_u8, 1_u8)),
        AromaticSystemUpdate {
            electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])),
            charge: Some(NumForm::Undetermined),
            unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) },
            constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)),
        },
    )]
    fn test_aromatic_system_form_difference_to(
        #[case] system: AromaticSystemForm,
        #[case] other: AromaticSystemForm,
        #[case] expected: AromaticSystemUpdate,
    ) {
        assert_eq!(system.difference_to(&other), expected);
    }

    #[rstest]
    #[case::normalized(
        AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(1_i64),
        AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(NumForm::lit_set([1])),
    )]
    fn test_aromatic_system_form_difference_to_identity(
        #[case] system: AromaticSystemForm,
        #[case] other: AromaticSystemForm,
    ) {
        assert_eq!(
            system.difference_to(&other),
            AromaticSystemUpdate::default()
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::reorder(
        AromaticSystemForm::from_electrons(vec![10, 20, 30]).with_charge(-1),
        &[AtomId(4), AtomId(7), AtomId(9)], &[AtomId(9), AtomId(4), AtomId(7)],
        Some(AromaticSystemForm::from_electrons(vec![30, 10, 20]).with_charge(-1)),
    )]
    #[case::undetermined_electrons(
        AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_charge(-1),
        &[AtomId(4), AtomId(7)], &[AtomId(7), AtomId(4)],
        Some(AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_charge(-1)),
    )]
    #[case::constraint_carries(
        AromaticSystemForm::from_electrons(vec![10, 20]).with_constraint(AromaticSystemConstraintForm::electron_count(30)),
        &[AtomId(4), AtomId(7)], &[AtomId(7), AtomId(4)],
        Some(AromaticSystemForm::from_electrons(vec![20, 10]).with_constraint(AromaticSystemConstraintForm::electron_count(30))),
    )]
    #[case::not_a_reordering(
        AromaticSystemForm::from_electrons(vec![10, 20]),
        &[AtomId(4), AtomId(7)], &[AtomId(7), AtomId(5)],
        None,
    )]
    #[case::frames_differ_in_length(
        AromaticSystemForm::from_electrons(vec![10, 20]),
        &[AtomId(4), AtomId(7)], &[AtomId(7)],
        None,
    )]
    fn test_aromatic_system_form_reframe_to(
        #[case] input: AromaticSystemForm,
        #[case] from: &[AtomId],
        #[case] to: &[AtomId],
        #[case] expected: Option<AromaticSystemForm>,
    ) {
        assert_eq!(input.reframe_to(from, to), expected);
    }

    #[rstest]
    #[case::identity_frame(AromaticSystemForm::from_electrons(vec![10, 20, 30]).with_charge(-1))]
    #[case::undetermined_electrons(AromaticSystemForm::new(ElectronCountsForm::Undetermined))]
    fn test_aromatic_system_form_reframe_to_identity(#[case] input: AromaticSystemForm) {
        let frame = [AtomId(4), AtomId(7), AtomId(9)];
        assert_eq!(input.clone().reframe_to(&frame, &frame), Some(input));
    }

    #[rstest]
    #[case::three_members(
        AromaticSystemForm::from_electrons(vec![10, 20, 30]).with_charge(-1),
        vec![
            ParticipantPosition(2),
            ParticipantPosition(0),
            ParticipantPosition(1),
        ],
        AromaticSystemForm::from_electrons(vec![30, 10, 20]).with_charge(-1),
    )]
    fn test_aromatic_system_form_permute(
        #[case] mut input: AromaticSystemForm,
        #[case] order: Vec<ParticipantPosition>,
        #[case] expected: AromaticSystemForm,
    ) {
        input.permute(&order);
        assert_eq!(input, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::all_undetermined(AromaticSystemForm::default(), false)]
    #[case::charge_only(AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_charge(0), false)]
    #[case::ground_no_atoms(AromaticSystemForm::new(ElectronCountsForm::Lit(Vec::new())).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::all_ground_six(AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::ground_with_constraint(AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)).with_constraint(AromaticSystemConstraintForm::electron_count(6)), true)]
    fn test_aromatic_system_form_is_ground(
        #[case] form: AromaticSystemForm,
        #[case] expected: bool,
    ) {
        assert_eq!(form.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_charge(
        AromaticSystemForm::default().with_charge(NumForm::lit_set([0])),
        Ok(AromaticSystemForm::default().with_charge(0)),
    )]
    #[case::charge_empty_litset_contradiction(
        AromaticSystemForm::default().with_charge(NumForm::lit_set(Vec::<i64>::new())),
        Err(Contradiction),
    )]
    fn test_aromatic_system_form_normalize(
        #[case] input: AromaticSystemForm,
        #[case] expected: Result<AromaticSystemForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_default(AromaticSystemForm::default(), AromaticSystemForm::default(), true)]
    #[case::default_matches_ground(AromaticSystemForm::default(), AromaticSystemForm::new(ElectronCountsForm::Lit(Vec::new())).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::exact(AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)),
        AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::electrons_length_mismatch(AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 5])),
        AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)), false)]
    #[case::electrons_value_mismatch(AromaticSystemForm::new(ElectronCountsForm::Lit(vec![2; 6])),
        AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)), false)]
    #[case::pattern_undetermined_electron_matches_lit(AromaticSystemForm::new(ElectronCountsForm::Undetermined),
      AromaticSystemForm::new(ElectronCountsForm::Lit(vec![1; 6])).with_charge(0).with_unpaired_electrons((0, 1)), true)]
    #[case::charge_mismatch(AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_charge(1),
        AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_charge(0), false)]
    #[case::unpaired_electrons_mismatch(AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_unpaired_electrons((2_u8, 3_u8)),
        AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_unpaired_electrons((0_u8, 1_u8)), false)]
    #[case::constraint_required_present(
        AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_constraint(AromaticSystemConstraintForm::electron_count(6)),
        AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_constraint(AromaticSystemConstraintForm::electron_count(6)),
        true)]
    #[case::constraint_required_absent(
        AromaticSystemForm::new(ElectronCountsForm::Undetermined).with_constraint(AromaticSystemConstraintForm::electron_count(6)),
        AromaticSystemForm::new(ElectronCountsForm::Undetermined),
        false)]
    fn test_aromatic_system_form_matches(
        #[case] pattern: AromaticSystemForm,
        #[case] target: AromaticSystemForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(
        AromaticSystemForm::default(),
        AromaticSystemForm::default(),
        Some(AromaticSystemForm::default())
    )]
    #[case::electrons_length_mismatch(
        AromaticSystemForm::from_electrons(vec![1; 6]),
        AromaticSystemForm::from_electrons(vec![1; 5]),
        None,
    )]
    #[case::narrows_electrons(
        AromaticSystemForm::new(ElectronCountsForm::Undetermined),
        AromaticSystemForm::from_electrons(vec![1; 3]),
        Some(AromaticSystemForm::from_electrons(vec![1; 3])),
    )]
    fn test_aromatic_system_form_meet(
        #[case] a: AromaticSystemForm,
        #[case] b: AromaticSystemForm,
        #[case] expected: Option<AromaticSystemForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::electrons_length_mismatch_widens_to_default(
        AromaticSystemForm::from_electrons(vec![1; 6]),
        AromaticSystemForm::from_electrons(vec![1; 5]),
        AromaticSystemForm::default(),
    )]
    fn test_aromatic_system_form_join(
        #[case] a: AromaticSystemForm,
        #[case] b: AromaticSystemForm,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }
}
