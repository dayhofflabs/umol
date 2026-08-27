//! Dative bond form.

use std::sync::Arc;

use umol_graph_core::{
    BiRelationData, FixedVarBirelationSet, NodeId, Ordered, ParticipantPosition, RelationId,
    Remapping, Unordered,
};
use umol_graph_ir_macros::{Lattice, Normalize};

use super::constraint::{DativeBondConstraintForm, DativeBondConstraintsForm};
use super::delta::EntitySpan;
use super::error::Contradiction;
use super::id::{AtomId, DativeBondId};
use super::num::NumForm;
use super::traits::{Equiv, Lattice, Normalize, Reframe};

/// The site this entry is borne by.
/// The molecule's dative bonds.
///
/// The donors bear the participant frame; the acceptor is a single distinguished atom. Neither
/// factor is position-sensitive for [`DativeBondForm`], whose order is a scalar.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct DativeBonds(
    Arc<FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>>,
);

impl DativeBonds {
    pub fn new(entries: Vec<(Vec<AtomId>, AtomId, DativeBondForm)>) -> Self {
        Self(Arc::new(FixedVarBirelationSet::new(
            entries
                .into_iter()
                .map(|(donors, acceptor, attributes)| {
                    (
                        [NodeId::from(acceptor)],
                        donors.into_iter().map(NodeId::from).collect(),
                        attributes,
                    )
                })
                .collect(),
        )))
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: DativeBondId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = DativeBondId> {
        self.0.ids().map(DativeBondId::from)
    }

    /// The atom accepting the donated pairs.
    pub fn acceptor(&self, id: DativeBondId) -> AtomId {
        AtomId::from(self.0.participants_1(RelationId::from(id))[0])
    }

    /// The donating atoms of `id`, in their stored frame.
    pub fn donors(&self, id: DativeBondId) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.0
            .participants_2(RelationId::from(id))
            .iter()
            .map(|&atom| AtomId::from(atom))
    }

    pub fn attributes(&self, id: DativeBondId) -> &DativeBondForm {
        self.0.data(RelationId::from(id))
    }

    pub fn attributes_mut(&mut self, id: DativeBondId) -> &mut DativeBondForm {
        Arc::make_mut(&mut self.0).data_mut(RelationId::from(id))
    }

    /// Ids of the dative bonds `atom` takes part in, as acceptor or donor.
    pub fn incident_ids(&self, atom: AtomId) -> impl ExactSizeIterator<Item = DativeBondId> + '_ {
        self.0
            .incident(NodeId::from(atom))
            .iter()
            .map(|&id| DativeBondId::from(id))
    }

    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.0.has_incident(NodeId::from(atom))
    }

    pub fn into_entries(self) -> Vec<(Vec<AtomId>, AtomId, DativeBondForm)> {
        Arc::try_unwrap(self.0)
            .unwrap_or_else(|shared| (*shared).clone())
            .into_entries()
            .into_iter()
            .map(|(acceptor, donors, attributes)| {
                (
                    donors.into_iter().map(AtomId::from).collect(),
                    AtomId::from(acceptor[0]),
                    attributes,
                )
            })
            .collect()
    }

    /// The acceptor of `id` as a graph node, for graph-core interop that is not yet typed in
    /// graph-IR ids. The public accessor is [`Self::acceptor`].
    pub(crate) fn acceptor_node(&self, id: DativeBondId) -> NodeId {
        self.0.participants_1(RelationId::from(id))[0]
    }

    /// The donors of `id` as graph nodes, for graph-core interop that is not yet typed in graph-IR
    /// ids. The public accessor is [`Self::donors`].
    pub(crate) fn donor_nodes(&self, id: DativeBondId) -> &[NodeId] {
        self.0.participants_2(RelationId::from(id))
    }

    pub(crate) fn attributes_iter_mut(
        &mut self,
    ) -> impl ExactSizeIterator<Item = &mut DativeBondForm> {
        Arc::make_mut(&mut self.0)
            .iter_mut()
            .map(|(_, _, _, attributes)| attributes)
    }

    pub(crate) fn remap(&self, remapping: &Remapping) -> Self {
        Self(Arc::new(self.0.remap(remapping)))
    }

    pub(crate) fn into_arc(
        self,
    ) -> Arc<FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>> {
        self.0
    }

    /// Glue `right`, relabelled into this molecule's id space, onto `self`: coinciding bonds meet,
    /// non-coinciding bonds are carried. `None` when a coincident meet is bottom.
    pub(crate) fn glue(&self, right: &Self, remapping: &Remapping) -> Option<Self> {
        self.0
            .pushout(
                &right.remap(remapping).0,
                // The acceptor is one atom, the sharpest node anchor a dative bond has.
                |set, acceptor, donors| {
                    acceptor
                        .first()
                        .and_then(|&node| set.coincident(node, acceptor, donors))
                },
                |(_, _, left), (_, _, right)| {
                    // The payload is frame-invariant, so the donor presentation cannot affect it.
                    right.clone().meet(left)
                },
            )
            .map(|merged| Self(Arc::new(merged.object)))
    }

    /// Id of the entity coinciding with these participants — the one whose participants equal
    /// them as a multiset. The identity question, distinct from lookup.
    pub fn coincident_id(&self, acceptor: AtomId, donors: &[AtomId]) -> Option<DativeBondId> {
        // The acceptor is a single atom, so it is the sharpest node anchor available.
        let donors: Vec<NodeId> = donors.iter().map(|&atom| NodeId::from(atom)).collect();
        self.0
            .coincident(NodeId::from(acceptor), &[NodeId::from(acceptor)], &donors)
            .map(DativeBondId::from)
    }
}

impl From<FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>>
    for DativeBonds
{
    fn from(
        set: FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>,
    ) -> Self {
        Self(Arc::new(set))
    }
}

impl From<Arc<FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>>>
    for DativeBonds
{
    fn from(
        set: Arc<FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondForm>>,
    ) -> Self {
        Self(set)
    }
}

impl Reframe for DativeBonds {
    type Action = (DativeBondId, Vec<ParticipantPosition>);

    /// Reduce every entry, then present each in its selected frame, returning the frame action selected for each bond.
    ///
    /// The payload is frame-invariant, so the action records the reordering of the donors without
    /// changing what the bond means.
    fn reframe_with_action(&self) -> Result<(Self, Vec<Self::Action>), Contradiction> {
        let mut reframed = (*self.0).clone();
        let mut actions = Vec::with_capacity(reframed.count());
        for id in reframed.ids().collect::<Vec<_>>() {
            let stored: Vec<AtomId> = reframed
                .participants_2(id)
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
            reframed.permute_2_with(id, &order);
            actions.push((DativeBondId::from(id), order));
        }
        Ok((Self(Arc::new(reframed)), actions))
    }
}

/// The reaction span's dative bonds, one [`EntitySpan`] per entity against a single donor frame.
/// The `Molecule` peer is [`DativeBonds`].
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct DativeBondSpans(
    FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, EntitySpan<DativeBondForm>>,
);

impl DativeBondSpans {
    pub fn into_entries(self) -> Vec<(Vec<AtomId>, AtomId, EntitySpan<DativeBondForm>)> {
        self.0
            .into_entries()
            .into_iter()
            .map(|(acceptor, donors, span)| {
                (
                    donors.into_iter().map(AtomId::from).collect(),
                    AtomId::from(acceptor[0]),
                    span,
                )
            })
            .collect()
    }

    pub fn new(entries: Vec<(Vec<AtomId>, AtomId, EntitySpan<DativeBondForm>)>) -> Self {
        Self(FixedVarBirelationSet::new(
            entries
                .into_iter()
                .map(|(donors, acceptor, span)| {
                    (
                        [NodeId::from(acceptor)],
                        donors.into_iter().map(NodeId::from).collect(),
                        span,
                    )
                })
                .collect(),
        ))
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: DativeBondId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = DativeBondId> {
        self.0.ids().map(DativeBondId::from)
    }

    /// The atom accepting the donated pairs. Not frame-bearing.
    pub fn acceptor(&self, id: DativeBondId) -> AtomId {
        AtomId::from(self.0.participants_1(RelationId::from(id))[0])
    }

    /// The donors of `id`, in their stored frame.
    pub fn donors(&self, id: DativeBondId) -> impl ExactSizeIterator<Item = AtomId> + '_ {
        self.0
            .participants_2(RelationId::from(id))
            .iter()
            .map(|&atom| AtomId::from(atom))
    }

    pub fn attributes(&self, id: DativeBondId) -> &EntitySpan<DativeBondForm> {
        self.0.data(RelationId::from(id))
    }

    pub(crate) fn remap(&self, remapping: &Remapping) -> Self {
        Self(self.0.remap(remapping))
    }
}

impl Reframe for DativeBondSpans {
    type Action = (DativeBondId, Vec<ParticipantPosition>);

    /// The payload is frame-invariant and selection sorts the donors, so a `Modified` span needs no
    /// arbitration between its sides.
    fn reframe_with_action(&self) -> Result<(Self, Vec<Self::Action>), Contradiction> {
        let mut reframed = self.0.clone();
        let mut actions = Vec::with_capacity(reframed.count());
        for id in reframed.ids().collect::<Vec<_>>() {
            let stored: Vec<AtomId> = reframed
                .participants_2(id)
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
            reframed.permute_2_with(id, &order);
            actions.push((DativeBondId::from(id), order));
        }
        Ok((Self(reframed), actions))
    }
}

/// Dative bond data: bond order (number of electron pairs donated) and
/// constraints. The acceptor and donor atoms are the participants of the
/// owning molecule's dative relation; this value holds no participant refs.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Normalize, Lattice)]
pub struct DativeBondForm {
    /// Bond order — number of electron pairs donated.
    pub order: NumForm,
    pub constraints: DativeBondConstraintsForm,
}

/// Attribute update for a dative bond. The order is optional, and
/// undetermined constraints remove their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DativeBondUpdate {
    pub order: Option<NumForm>,
    pub constraints: DativeBondConstraintsForm,
}

impl DativeBondForm {
    /// Concrete: every inherent field is ground; the constraint channel does
    /// not bear on concreteness.
    pub fn is_concrete(&self) -> bool {
        let DativeBondForm {
            order,
            constraints: _,
        } = self;
        order.is_ground()
    }
    pub fn new(order: NumForm) -> Self {
        Self {
            order,
            constraints: DativeBondConstraintsForm::new(),
        }
    }

    pub fn from_order(order: u8) -> Self {
        Self::new(NumForm::Lit(order as i64))
    }

    pub fn with_order(mut self, order: impl Into<NumForm>) -> Self {
        self.order = order.into();
        self
    }

    /// Add a single constraint, replacing any existing entry of the same
    /// key (last-wins per `DativeBondConstraintsForm::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<DativeBondConstraintForm>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same key (last-wins per `DativeBondConstraintsForm::set`). Does
    /// not clear existing constraints; use `bond.constraints.clear()` or
    /// direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<DativeBondConstraintForm>,
    {
        self.constraints
            .extend(constraints.into_iter().map(Into::into));
        self
    }

    /// No-op: `DativeBondForm` has no value-bearing fields besides `order`,
    /// which is essential and never filled. Provided for API symmetry.
    pub fn into_concrete(self) -> Self {
        self
    }

    /// Apply an attribute update, leaving omitted leaves and constraint keys unchanged.
    pub fn update(&self, update: &DativeBondUpdate) -> DativeBondForm {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        DativeBondForm {
            order: update.order.clone().unwrap_or_else(|| self.order.clone()),
            constraints,
        }
    }

    /// Derive the minimal normalized attribute update carrying `self` to `other`.
    /// Frame-invariant: no field is position-indexed, so a frame change carries the form
    /// unchanged.
    ///
    /// Destructured exhaustively on purpose: a new position-indexed field must fail to compile
    /// here rather than be silently left in the old frame.
    pub fn reframe_to(self, _from: &[AtomId], _to: &[AtomId]) -> Option<Self> {
        let Self { order, constraints } = self;
        Some(Self { order, constraints })
    }

    pub fn difference_to(&self, other: &Self) -> DativeBondUpdate {
        let mut constraints = DativeBondConstraintsForm::new();
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
        DativeBondUpdate {
            order: (!self.order.equiv(&other.order)).then(|| other.order.clone()),
            constraints,
        }
    }
}

impl From<&str> for DativeBondForm {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid dative bond string")
    }
}

impl BiRelationData for DativeBondForm {
    /// `order` is a scalar; neither the acceptor nor the donor factor is position-indexed.
    fn on_permutation(
        &mut self,
        _order_1: &[ParticipantPosition],
        _order_2: &[ParticipantPosition],
    ) {
    }

    fn is_permutation_invariant(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ir::boolean::BooleanForm;
    use crate::ir::constraint::RingScope;
    use crate::ir::error::Contradiction;
    use crate::ir::traits::{Lattice, Normalize};

    #[rustfmt::skip]
    /// The payload is frame-invariant, so a `Modified` span's two sides carry unchanged through the
    /// selected action while the donors sort. The acceptor bears no frame and does not move.
    #[rstest]
    fn test_dative_bond_spans_reframe() {
        let span = EntitySpan::Modified {
            lhs: DativeBondForm::from_order(1),
            rhs: DativeBondForm::from_order(2),
        };
        let mut spans = DativeBondSpans::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            AtomId(0),
            span.clone(),
        )]);
        spans.0.permute_2_with(
            RelationId(0),
            &[
                ParticipantPosition(2),
                ParticipantPosition(0),
                ParticipantPosition(1),
            ],
        );

        let (reframed, actions) = spans.reframe_with_action().expect("the forms are satisfiable");

        assert_eq!(
            reframed.donors(DativeBondId(0)).collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4), AtomId(7)],
        );
        assert_eq!(reframed.acceptor(DativeBondId(0)), AtomId(0));
        assert_eq!(reframed.attributes(DativeBondId(0)), &span);
        assert_eq!(
            actions,
            vec![(
                DativeBondId(0),
                vec![
                    ParticipantPosition(1),
                    ParticipantPosition(2),
                    ParticipantPosition(0),
                ],
            )],
        );
    }

    #[rstest]
    fn test_dative_bond_spans_reframe_identity() {
        let spans = DativeBondSpans::new(vec![(
            vec![AtomId(1), AtomId(4)],
            AtomId(0),
            EntitySpan::Modified {
                lhs: DativeBondForm::from_order(1),
                rhs: DativeBondForm::from_order(2),
            },
        )]);
        let once = spans.reframe().expect("the forms are satisfiable");
        let twice = once.reframe().expect("the forms are satisfiable");
        assert_eq!(twice, once);
    }

    /// The donors are an `Unordered` factor, so storage sorts them on construction; the stored
    /// frame is permuted first to model the frame-preserving storage S5 introduces.
    #[fixture]
    fn unsorted_bond() -> DativeBonds {
        let mut bonds = DativeBonds::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            AtomId(0),
            DativeBondForm::from_order(1),
        )]);
        Arc::make_mut(&mut bonds.0).permute_2_with(
            RelationId(0),
            &[
                ParticipantPosition(2),
                ParticipantPosition(0),
                ParticipantPosition(1),
            ],
        );
        bonds
    }

    #[rstest]
    fn test_dative_bonds_reframe(unsorted_bond: DativeBonds) {
        assert_eq!(
            unsorted_bond.donors(DativeBondId(0)).collect::<Vec<_>>(),
            vec![AtomId(7), AtomId(1), AtomId(4)],
        );

        let reframed = unsorted_bond.reframe().expect("the form is satisfiable");

        assert_eq!(
            reframed.donors(DativeBondId(0)).collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4), AtomId(7)],
        );
        assert_eq!(reframed.acceptor(DativeBondId(0)), AtomId(0));
        assert_eq!(
            reframed.attributes(DativeBondId(0)),
            &DativeBondForm::from_order(1),
        );
    }

    #[rstest]
    fn test_dative_bonds_reframe_identity(unsorted_bond: DativeBonds) {
        let once = unsorted_bond.reframe().expect("the form is satisfiable");
        let twice = once.reframe().expect("the form is satisfiable");
        assert_eq!(twice, once);
    }

    #[rstest]
    fn test_dative_bonds_reframe_with_action(unsorted_bond: DativeBonds) {
        let (_, actions) = unsorted_bond
            .reframe_with_action()
            .expect("the form is satisfiable");
        assert_eq!(
            actions,
            vec![(
                DativeBondId(0),
                vec![
                    ParticipantPosition(1),
                    ParticipantPosition(2),
                    ParticipantPosition(0),
                ],
            )],
        );
    }

    #[rstest]
    fn test_dative_bonds_framed_eq(unsorted_bond: DativeBonds) {
        let selected = DativeBonds::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            AtomId(0),
            DativeBondForm::from_order(1),
        )]);
        assert!(unsorted_bond.framed_eq(&selected));

        let different = DativeBonds::new(vec![(
            vec![AtomId(1), AtomId(4), AtomId(7)],
            AtomId(0),
            DativeBondForm::from_order(2),
        )]);
        assert!(!unsorted_bond.framed_eq(&different));
    }

    #[rstest]
    #[case::new(DativeBondForm::new(NumForm::Lit(2)),
        DativeBondForm { order: NumForm::Lit(2), constraints: DativeBondConstraintsForm::new() })]
    #[case::from_order(DativeBondForm::from_order(3),
        DativeBondForm { order: NumForm::Lit(3), constraints: DativeBondConstraintsForm::new() })]
    fn test_dative_bond_form_new(#[case] actual: DativeBondForm, #[case] expected: DativeBondForm) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_order(DativeBondForm::default().with_order(2),
        DativeBondForm { order: NumForm::Lit(2), constraints: DativeBondConstraintsForm::new() })]
    #[case::with_constraint(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        DativeBondForm { order: NumForm::Lit(1),
            constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))) })]
    #[case::with_constraints_extends(
        DativeBondForm::from_order(1)
            .with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)))
            .with_constraints([DativeBondConstraintForm::ring_membership(RingScope::All, 1), DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)]),
        DativeBondForm { order: NumForm::Lit(1),
            constraints: DativeBondConstraintsForm::from_iter([
                DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                DativeBondConstraintForm::ring_membership(RingScope::All, 1),
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ]) })]
    #[case::with_constraint_appends_different_scopes(
        DativeBondForm::from_order(1)
            .with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1))
            .with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)),
        DativeBondForm { order: NumForm::Lit(1),
            constraints: DativeBondConstraintsForm::from_iter([
                DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1),
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ]) })]
    fn test_dative_bond_form_with_methods(#[case] actual: DativeBondForm, #[case] expected: DativeBondForm) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::from_order(DativeBondForm::from_order(1))]
    #[case::with_constraint(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    fn test_dative_bond_form_into_concrete(#[case] bond: DativeBondForm) {
        assert_eq!(bond.clone().into_concrete(), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::order(DativeBondForm::from_order(1), DativeBondUpdate { order: Some(NumForm::Lit(2)), ..Default::default() }, DativeBondForm::from_order(2))]
    #[case::order_undetermined(DativeBondForm::from_order(1), DativeBondUpdate { order: Some(NumForm::Undetermined), ..Default::default() }, DativeBondForm::default())]
    #[case::constraint_set(DativeBondForm::from_order(1), DativeBondUpdate { constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))), ..Default::default() }, DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    #[case::constraint_replace(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1_i64)), DativeBondUpdate { constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 2_i64)), ..Default::default() }, DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 2_i64)))]
    #[case::constraint_remove(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1_i64)), DativeBondUpdate { constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() }, DativeBondForm::from_order(1))]
    fn test_dative_bond_form_update(
        #[case] bond: DativeBondForm,
        #[case] update: DativeBondUpdate,
        #[case] expected: DativeBondForm,
    ) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    fn test_dative_bond_form_update_identity(#[case] bond: DativeBondForm) {
        assert_eq!(bond.update(&DativeBondUpdate::default()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ground(DativeBondForm::from_order(1))]
    #[case::with_constraint(DativeBondForm::from_order(1)
        .with_constraint(DativeBondConstraintForm::aromatic(true)))]
    #[case::undetermined(DativeBondForm::default())]
    fn test_dative_bond_form_reframe_to(#[case] input: DativeBondForm) {
        assert_eq!(
            input.clone().reframe_to(&[AtomId(4), AtomId(7)], &[AtomId(7), AtomId(4)]),
            Some(input),
        );
    }

    #[rstest]
    #[case::fields_and_constraints(
        DativeBondForm::from_order(1).with_constraints([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1_i64),
        ]),
        DativeBondForm::from_order(2).with_constraints([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false)),
            DativeBondConstraintForm::ring_membership(RingScope::All, 2_i64),
        ]),
        DativeBondUpdate {
            order: Some(NumForm::Lit(2)),
            constraints: DativeBondConstraintsForm::from_iter([
                DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false)),
                DativeBondConstraintForm::ring_membership(RingScope::All, 2_i64),
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined),
            ]),
        },
    )]
    fn test_dative_bond_form_difference_to(
        #[case] bond: DativeBondForm,
        #[case] other: DativeBondForm,
        #[case] expected: DativeBondUpdate,
    ) {
        assert_eq!(bond.difference_to(&other), expected);
    }

    #[rstest]
    #[case::normalized(
        DativeBondForm::from_order(1),
        DativeBondForm::new(NumForm::lit_set([1])),
    )]
    fn test_dative_bond_form_difference_to_identity(
        #[case] bond: DativeBondForm,
        #[case] other: DativeBondForm,
    ) {
        assert_eq!(bond.difference_to(&other), DativeBondUpdate::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(DativeBondForm::default(), false)]
    #[case::order_lit(DativeBondForm::from_order(1), true)]
    #[case::order_undetermined(DativeBondForm::new(NumForm::Undetermined), false)]
    #[case::ground_with_constraint(DativeBondForm { order: NumForm::Lit(1),
        constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)) }, true)]
    fn test_dative_bond_form_is_ground(#[case] form: DativeBondForm, #[case] expected: bool) {
        assert_eq!(form.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_order(
        DativeBondForm::new(NumForm::lit_set([1])),
        Ok(DativeBondForm::from_order(1)),
    )]
    #[case::order_empty_litset_contradiction(
        DativeBondForm::new(NumForm::lit_set(Vec::<i64>::new())),
        Err(Contradiction),
    )]
    fn test_dative_bond_form_normalize(
        #[case] input: DativeBondForm,
        #[case] expected: Result<DativeBondForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rstest]
    #[case::same_order(DativeBondForm::from_order(1), DativeBondForm::from_order(1), true)]
    #[case::order_mismatch(DativeBondForm::from_order(1), DativeBondForm::from_order(2), false)]
    #[case::pattern_wildcard_matches(
        DativeBondForm::new(NumForm::Undetermined),
        DativeBondForm::from_order(2),
        true
    )]
    #[case::pattern_more_specific(
        DativeBondForm::from_order(2),
        DativeBondForm::new(NumForm::Undetermined),
        false
    )]
    #[case::constraint_required_present(
        DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        true
    )]
    #[case::constraint_required_absent(
        DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        DativeBondForm::from_order(1),
        false
    )]
    fn test_dative_bond_form_matches(
        #[case] pattern: DativeBondForm,
        #[case] target: DativeBondForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(
        DativeBondForm::default(),
        DativeBondForm::default(),
        Some(DativeBondForm::default())
    )]
    #[case::order_mismatch(DativeBondForm::from_order(1), DativeBondForm::from_order(2), None)]
    #[case::narrows_order(
        DativeBondForm { order: NumForm::Undetermined, constraints: DativeBondConstraintsForm::new() },
        DativeBondForm { order: NumForm::Lit(1), constraints: DativeBondConstraintsForm::new() },
        Some(DativeBondForm { order: NumForm::Lit(1), constraints: DativeBondConstraintsForm::new() }),
    )]
    fn test_dative_bond_form_meet(
        #[case] a: DativeBondForm,
        #[case] b: DativeBondForm,
        #[case] expected: Option<DativeBondForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::order_widens_to_lit_set(
        DativeBondForm::from_order(1),
        DativeBondForm::from_order(2),
        DativeBondForm { order: NumForm::lit_set([1, 2]), constraints: DativeBondConstraintsForm::new() },
    )]
    fn test_dative_bond_form_join(
        #[case] a: DativeBondForm,
        #[case] b: DativeBondForm,
        #[case] expected: DativeBondForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }
}
