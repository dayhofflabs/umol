//! Noncovalent bond form.

use std::borrow::Cow;
use std::sync::Arc;

use umol_graph_core::{
    FixedRelationSet, NodeId, ParticipantPosition, RelationData, RelationId, Remapping, Unordered,
};
use umol_graph_ir_macros::{Lattice, Normalize};

use super::constraint::{NoncovalentBondConstraintForm, NoncovalentBondConstraintsForm};
use super::delta::EntitySpan;
use super::error::{Contradiction, NoJoin};
use super::id::{AtomId, NoncovalentBondId};
use super::traits::{AsLit, Equiv, Lattice, Normalize, Reframe};

/// The molecule's noncovalent bonds.
///
/// The atom pair bears the participant frame, but [`NoncovalentBondForm`] carries no
/// position-sensitive field, so a reordering of the pair leaves the attributes unchanged.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct NoncovalentBonds(Arc<FixedRelationSet<NodeId, Unordered, NoncovalentBondForm, 2>>);

impl NoncovalentBonds {
    pub fn new(entries: Vec<(AtomId, AtomId, NoncovalentBondForm)>) -> Self {
        Self(Arc::new(FixedRelationSet::new(
            entries
                .into_iter()
                .map(|(first, second, attributes)| {
                    ([NodeId::from(first), NodeId::from(second)], attributes)
                })
                .collect(),
        )))
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: NoncovalentBondId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = NoncovalentBondId> {
        self.0.ids().map(NoncovalentBondId::from)
    }

    /// The bonded pair of `id`, in its stored frame.
    pub fn atoms(&self, id: NoncovalentBondId) -> [AtomId; 2] {
        self.0.participants(RelationId::from(id)).map(AtomId::from)
    }

    pub fn attributes(&self, id: NoncovalentBondId) -> &NoncovalentBondForm {
        self.0.data(RelationId::from(id))
    }

    pub fn attributes_mut(&mut self, id: NoncovalentBondId) -> &mut NoncovalentBondForm {
        Arc::make_mut(&mut self.0).data_mut(RelationId::from(id))
    }

    /// Ids of the noncovalent bonds `atom` takes part in. Integrity rejects a parallel pair, so an
    /// atom pairs at most once with any given partner but may bond several partners.
    pub fn incident_ids(
        &self,
        atom: AtomId,
    ) -> impl ExactSizeIterator<Item = NoncovalentBondId> + '_ {
        self.0
            .incident(NodeId::from(atom))
            .iter()
            .map(|&id| NoncovalentBondId::from(id))
    }

    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.0.has_incident(NodeId::from(atom))
    }

    pub fn into_entries(self) -> Vec<(AtomId, AtomId, NoncovalentBondForm)> {
        Arc::try_unwrap(self.0)
            .unwrap_or_else(|shared| (*shared).clone())
            .into_entries()
            .into_iter()
            .map(|([first, second], attributes)| {
                (AtomId::from(first), AtomId::from(second), attributes)
            })
            .collect()
    }

    pub(crate) fn attributes_iter_mut(
        &mut self,
    ) -> impl ExactSizeIterator<Item = &mut NoncovalentBondForm> {
        Arc::make_mut(&mut self.0)
            .iter_mut()
            .map(|(_, _, attributes)| attributes)
    }

    pub(crate) fn remap(&self, remapping: &Remapping) -> Self {
        Self(Arc::new(self.0.remap(remapping)))
    }

    pub(crate) fn into_arc(
        self,
    ) -> Arc<FixedRelationSet<NodeId, Unordered, NoncovalentBondForm, 2>> {
        self.0
    }

    /// Glue `right`, relabelled into this molecule's id space, onto `self`: coinciding bonds meet,
    /// non-coinciding bonds are carried. `None` when a coincident meet is bottom.
    pub(crate) fn glue(&self, right: &Self, remapping: &Remapping) -> Option<Self> {
        self.0
            .pushout(&right.remap(remapping).0, |(_, left), (_, right)| {
                // The payload is frame-invariant, so the pair's presentation cannot affect it.
                right.clone().meet(left)
            })
            .map(|merged| Self(Arc::new(merged.object)))
    }

    /// Id of the entity coinciding with these participants — the one whose participants equal
    /// them as a multiset. The identity question, distinct from lookup.
    pub fn coincident_id(&self, first: AtomId, second: AtomId) -> Option<NoncovalentBondId> {
        self.0
            .coincident(&[NodeId::from(first), NodeId::from(second)])
            .map(NoncovalentBondId::from)
    }

    /// Superseded once frame alignment moves onto `reframe_to`.
    pub(crate) fn participant_permutation(
        &self,
        id: NoncovalentBondId,
        atoms: &[AtomId],
    ) -> Option<Vec<ParticipantPosition>> {
        let query: Vec<NodeId> = atoms.iter().map(|&atom| NodeId::from(atom)).collect();
        self.0.participant_permutation(RelationId::from(id), &query)
    }
}

impl From<FixedRelationSet<NodeId, Unordered, NoncovalentBondForm, 2>> for NoncovalentBonds {
    fn from(set: FixedRelationSet<NodeId, Unordered, NoncovalentBondForm, 2>) -> Self {
        Self(Arc::new(set))
    }
}

impl From<Arc<FixedRelationSet<NodeId, Unordered, NoncovalentBondForm, 2>>> for NoncovalentBonds {
    fn from(set: Arc<FixedRelationSet<NodeId, Unordered, NoncovalentBondForm, 2>>) -> Self {
        Self(set)
    }
}

impl Reframe for NoncovalentBonds {
    type Action = (NoncovalentBondId, Vec<ParticipantPosition>);

    /// Reduce every entry, then present each in its selected frame, returning the frame action selected for each bond.
    ///
    /// The payload is frame-invariant, so the action records the reordering without changing what
    /// the bond means.
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
            actions.push((NoncovalentBondId::from(id), order));
        }
        Ok((Self(Arc::new(reframed)), actions))
    }
}

/// The reaction span's noncovalent bonds, one [`EntitySpan`] per entity against a single
/// participant frame. The `Molecule` peer is [`NoncovalentBonds`].
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct NoncovalentBondSpans(
    FixedRelationSet<NodeId, Unordered, EntitySpan<NoncovalentBondForm>, 2>,
);

impl NoncovalentBondSpans {
    pub fn into_entries(self) -> Vec<(AtomId, AtomId, EntitySpan<NoncovalentBondForm>)> {
        self.0
            .into_entries()
            .into_iter()
            .map(|([first, second], span)| (AtomId::from(first), AtomId::from(second), span))
            .collect()
    }

    pub fn new(entries: Vec<(AtomId, AtomId, EntitySpan<NoncovalentBondForm>)>) -> Self {
        Self(FixedRelationSet::new(
            entries
                .into_iter()
                .map(|(first, second, span)| ([NodeId::from(first), NodeId::from(second)], span))
                .collect(),
        ))
    }

    pub fn count(&self) -> usize {
        self.0.count()
    }

    pub fn contains(&self, id: NoncovalentBondId) -> bool {
        self.0.contains(RelationId::from(id))
    }

    pub fn ids(&self) -> impl ExactSizeIterator<Item = NoncovalentBondId> {
        self.0.ids().map(NoncovalentBondId::from)
    }

    /// The two atoms of `id`, in their stored frame.
    pub fn atoms(&self, id: NoncovalentBondId) -> [AtomId; 2] {
        self.0.participants(RelationId::from(id)).map(AtomId::from)
    }

    pub fn attributes(&self, id: NoncovalentBondId) -> &EntitySpan<NoncovalentBondForm> {
        self.0.data(RelationId::from(id))
    }

    pub(crate) fn remap(&self, remapping: &Remapping) -> Self {
        Self(self.0.remap(remapping))
    }
}

impl Reframe for NoncovalentBondSpans {
    type Action = (NoncovalentBondId, Vec<ParticipantPosition>);

    /// The payload is frame-invariant and selection sorts the pair, so a `Modified` span needs no
    /// arbitration between its sides.
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
            actions.push((NoncovalentBondId::from(id), order));
        }
        Ok((Self(reframed), actions))
    }
}

/// Noncovalent bond: two-atom non-bonded interaction tagged by an
/// interaction kind. No bond order, no charge or spin — these do not apply
/// to noncovalent interactions.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Normalize, Lattice)]
pub struct NoncovalentBondForm {
    pub kind: NoncovalentBondKindForm,
    pub constraints: NoncovalentBondConstraintsForm,
}

/// Attribute update for a noncovalent bond. The kind is optional, and an
/// undetermined constraint removes its key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NoncovalentBondUpdate {
    pub kind: Option<NoncovalentBondKindForm>,
    pub constraints: NoncovalentBondConstraintsForm,
}

impl From<&str> for NoncovalentBondForm {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid noncovalent bond string")
    }
}

impl RelationData for NoncovalentBondForm {
    /// The kind is not position-indexed — reordering the two participants leaves it unchanged.
    fn on_permutation(&mut self, _order: &[ParticipantPosition]) {}

    fn is_permutation_invariant(&self) -> bool {
        true
    }
}

impl NoncovalentBondForm {
    /// Concrete: every inherent field is ground; the constraint channel does
    /// not bear on concreteness.
    pub fn is_concrete(&self) -> bool {
        let NoncovalentBondForm {
            kind,
            constraints: _,
        } = self;
        kind.is_ground()
    }
    pub fn new(kind: NoncovalentBondKindForm) -> Self {
        Self {
            kind,
            constraints: NoncovalentBondConstraintsForm::new(),
        }
    }

    pub fn from_kind(kind: NoncovalentBondKind) -> Self {
        Self::new(NoncovalentBondKindForm::Lit(kind))
    }

    pub fn with_kind(mut self, kind: impl Into<NoncovalentBondKindForm>) -> Self {
        self.kind = kind.into();
        self
    }

    /// Add a single constraint, replacing any existing entry of the same
    /// kind (last-wins per `NoncovalentBondConstraintsForm::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<NoncovalentBondConstraintForm>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `NoncovalentBondConstraintsForm::set`).
    /// Does not clear existing constraints; use `bond.constraints.clear()`
    /// or direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<NoncovalentBondConstraintForm>,
    {
        for c in constraints {
            self.constraints.set(c.into());
        }
        self
    }

    /// No-op on value fields: `NoncovalentBondForm` has no value-bearing field
    /// besides `kind`, which is essential and never filled. Constraints are
    /// preserved. Provided for API symmetry.
    pub fn into_concrete(self) -> Self {
        self
    }

    /// Apply an attribute update, leaving an omitted kind and constraint keys unchanged.
    pub fn update(&self, update: &NoncovalentBondUpdate) -> NoncovalentBondForm {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        NoncovalentBondForm {
            kind: update.kind.clone().unwrap_or_else(|| self.kind.clone()),
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
        let Self { kind, constraints } = self;
        Some(Self { kind, constraints })
    }

    pub fn difference_to(&self, other: &Self) -> NoncovalentBondUpdate {
        let mut constraints = NoncovalentBondConstraintsForm::new();
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
        NoncovalentBondUpdate {
            kind: (!self.kind.equiv(&other.kind)).then(|| other.kind.clone()),
            constraints,
        }
    }
}

/// Noncovalent interaction kind.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NoncovalentBondKindForm {
    #[default]
    Undetermined,
    Lit(NoncovalentBondKind),
}

impl NoncovalentBondKindForm {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn lit(kind: NoncovalentBondKind) -> Self {
        Self::Lit(kind)
    }
}

impl From<NoncovalentBondKind> for NoncovalentBondKindForm {
    fn from(kind: NoncovalentBondKind) -> Self {
        Self::Lit(kind)
    }
}

impl Normalize for NoncovalentBondKindForm {
    /// Both variants are already normalized — nothing folds.
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(self)
    }

    fn normalized(&self) -> Result<Cow<'_, Self>, Contradiction> {
        Ok(Cow::Borrowed(self))
    }
}

impl AsLit for NoncovalentBondKindForm {
    type Lit = NoncovalentBondKind;

    /// The specific interaction kind, only when it is a literal.
    #[inline]
    fn as_lit(&self) -> Option<NoncovalentBondKind> {
        match self {
            Self::Lit(k) => Some(*k),
            Self::Undetermined => None,
        }
    }
}

impl Lattice for NoncovalentBondKindForm {
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    fn is_ground(&self) -> bool {
        matches!(self, Self::Lit(_))
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::Lit(a), Self::Lit(b)) if a == b => Some(Self::Lit(*a)),
            _ => None,
        }
    }

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        Ok(match (self, other) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::Lit(a), Self::Lit(b)) if a == b => Self::Lit(*a),
            _ => Self::Undetermined,
        })
    }
}

/// Fundamental kind of a noncovalent interaction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NoncovalentBondKind {
    HydrogenBond,
    HalogenBond,
    ChalcogenBond,
    Ionic,
    VanDerWaals,
}

#[cfg(test)]
mod tests {
    use std::iter;

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ir::boolean::BooleanForm;

    #[rustfmt::skip]
    /// The payload is frame-invariant, so a `Modified` span's two sides carry unchanged through the
    /// selected action while the pair itself sorts.
    #[rstest]
    fn test_noncovalent_bond_spans_reframe() {
        let span = EntitySpan::Modified {
            lhs: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            rhs: NoncovalentBondForm::from_kind(NoncovalentBondKind::HalogenBond),
        };
        let mut spans = NoncovalentBondSpans::new(vec![(AtomId(2), AtomId(5), span.clone())]);
        spans
            .0
            .permute_with(RelationId(0), &[ParticipantPosition(1), ParticipantPosition(0)]);

        let (reframed, actions) = spans.reframe_with_action().expect("the forms are satisfiable");

        assert_eq!(reframed.atoms(NoncovalentBondId(0)), [AtomId(2), AtomId(5)]);
        assert_eq!(reframed.attributes(NoncovalentBondId(0)), &span);
        assert_eq!(
            actions,
            vec![(
                NoncovalentBondId(0),
                vec![ParticipantPosition(1), ParticipantPosition(0)],
            )],
        );
    }

    #[rstest]
    fn test_noncovalent_bond_spans_reframe_identity() {
        let spans = NoncovalentBondSpans::new(vec![(
            AtomId(2),
            AtomId(5),
            EntitySpan::Modified {
                lhs: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                rhs: NoncovalentBondForm::from_kind(NoncovalentBondKind::HalogenBond),
            },
        )]);
        let once = spans.reframe().expect("the forms are satisfiable");
        let twice = once.reframe().expect("the forms are satisfiable");
        assert_eq!(twice, once);
    }

    /// Storage sorts an `Unordered` factor on construction, so the stored pair is permuted first to
    /// model the frame-preserving storage S5 introduces.
    #[fixture]
    fn unsorted_bond() -> NoncovalentBonds {
        let mut bonds = NoncovalentBonds::new(vec![(
            AtomId(2),
            AtomId(5),
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )]);
        Arc::make_mut(&mut bonds.0).permute_with(
            RelationId(0),
            &[ParticipantPosition(1), ParticipantPosition(0)],
        );
        bonds
    }

    #[rstest]
    fn test_noncovalent_bonds_reframe(unsorted_bond: NoncovalentBonds) {
        assert_eq!(
            unsorted_bond.atoms(NoncovalentBondId(0)),
            [AtomId(5), AtomId(2)],
        );

        let reframed = unsorted_bond.reframe().expect("the form is satisfiable");

        assert_eq!(reframed.atoms(NoncovalentBondId(0)), [AtomId(2), AtomId(5)]);
        assert_eq!(
            reframed.attributes(NoncovalentBondId(0)),
            &NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
    }

    #[rstest]
    fn test_noncovalent_bonds_reframe_identity(unsorted_bond: NoncovalentBonds) {
        let once = unsorted_bond.reframe().expect("the form is satisfiable");
        let twice = once.reframe().expect("the form is satisfiable");
        assert_eq!(twice, once);
    }

    #[rstest]
    fn test_noncovalent_bonds_reframe_with_action(unsorted_bond: NoncovalentBonds) {
        let (_, actions) = unsorted_bond
            .reframe_with_action()
            .expect("the form is satisfiable");
        assert_eq!(
            actions,
            vec![(
                NoncovalentBondId(0),
                vec![ParticipantPosition(1), ParticipantPosition(0)],
            )],
        );
    }

    #[rstest]
    fn test_noncovalent_bonds_framed_eq(unsorted_bond: NoncovalentBonds) {
        let selected = NoncovalentBonds::new(vec![(
            AtomId(2),
            AtomId(5),
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )]);
        assert!(unsorted_bond.framed_eq(&selected));

        let different = NoncovalentBonds::new(vec![(
            AtomId(2),
            AtomId(6),
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )]);
        assert!(!unsorted_bond.framed_eq(&different));
    }

    #[rstest]
    #[case::new(NoncovalentBondForm::new(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond)),
        NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraintsForm::new() })]
    #[case::from_kind(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraintsForm::new() })]
    fn test_noncovalent_bond_form_new(
        #[case] actual: NoncovalentBondForm,
        #[case] expected: NoncovalentBondForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_kind_primitive(
        NoncovalentBondForm::default().with_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraintsForm::new() })]
    #[case::with_kind_form(
        NoncovalentBondForm::default().with_kind(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond)),
        NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraintsForm::new() })]
    #[case::with_constraints_empty(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond)
            .with_constraints(iter::empty::<NoncovalentBondConstraintForm>()),
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))]
    #[case::with_constraint(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond)
            .with_constraint(NoncovalentBondConstraintForm::intramolecular(true)),
        NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)) })]
    #[case::with_constraints_populated(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond)
            .with_constraints([NoncovalentBondConstraintForm::intramolecular(false)]),
        NoncovalentBondForm { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(false)) })]
    fn test_noncovalent_bond_form_with_methods(
        #[case] actual: NoncovalentBondForm,
        #[case] expected: NoncovalentBondForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::default_(NoncovalentBondForm::default())]
    #[case::ground(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_form_into_concrete(#[case] bond: NoncovalentBondForm) {
        assert_eq!(bond.clone().into_concrete(), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondUpdate { kind: Some(NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic)), ..Default::default() },
        NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic))]
    #[case::kind_undetermined(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondUpdate { kind: Some(NoncovalentBondKindForm::Undetermined), ..Default::default() },
        NoncovalentBondForm::default())]
    #[case::constraint_set(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(true)), ..Default::default() },
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintForm::intramolecular(true)))]
    #[case::constraint_replace(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintForm::intramolecular(true)),
        NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(false)), ..Default::default() },
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintForm::intramolecular(false)))]
    #[case::constraint_remove(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintForm::intramolecular(true)),
        NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Undetermined)), ..Default::default() },
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_form_update(
        #[case] bond: NoncovalentBondForm,
        #[case] update: NoncovalentBondUpdate,
        #[case] expected: NoncovalentBondForm,
    ) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintForm::intramolecular(true)))]
    fn test_noncovalent_bond_form_update_identity(#[case] bond: NoncovalentBondForm) {
        assert_eq!(bond.update(&NoncovalentBondUpdate::default()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ground(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))]
    #[case::with_constraint(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond)
        .with_constraint(NoncovalentBondConstraintForm::intramolecular(true)))]
    #[case::undetermined(NoncovalentBondForm::default())]
    fn test_noncovalent_bond_form_reframe_to(#[case] input: NoncovalentBondForm) {
        assert_eq!(
            input.clone().reframe_to(&[AtomId(4), AtomId(7)], &[AtomId(7), AtomId(4)]),
            Some(input),
        );
    }

    #[rstest]
    #[case::kind_and_constraint(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintForm::intramolecular(true)),
        NoncovalentBondForm::default(),
        NoncovalentBondUpdate {
            kind: Some(NoncovalentBondKindForm::Undetermined),
            constraints: NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::Intramolecular(BooleanForm::Undetermined)),
        },
    )]
    fn test_noncovalent_bond_form_difference_to(
        #[case] bond: NoncovalentBondForm,
        #[case] other: NoncovalentBondForm,
        #[case] expected: NoncovalentBondUpdate,
    ) {
        assert_eq!(bond.difference_to(&other), expected);
    }

    #[rstest]
    #[case::same(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_form_difference_to_identity(#[case] bond: NoncovalentBondForm) {
        assert_eq!(bond.difference_to(&bond), NoncovalentBondUpdate::default());
    }

    #[rstest]
    #[case::default_(NoncovalentBondForm::default(), false)]
    #[case::ground_lit(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        true
    )]
    fn test_noncovalent_bond_form_is_ground(
        #[case] form: NoncovalentBondForm,
        #[case] expected: bool,
    ) {
        assert_eq!(form.is_ground(), expected);
    }

    #[rstest]
    #[case::ground(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))]
    #[case::undetermined(NoncovalentBondForm::default())]
    fn test_noncovalent_bond_form_normalize_identity(#[case] input: NoncovalentBondForm) {
        assert_eq!(input.clone().normalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_ground(NoncovalentBondForm::default(), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    #[case::same(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    #[case::different(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic), false)]
    #[case::pattern_specific_target_undetermined(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondForm::default(), false)]
    fn test_noncovalent_bond_form_matches(
        #[case] pattern: NoncovalentBondForm,
        #[case] target: NoncovalentBondForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(NoncovalentBondKindForm::Undetermined, NoncovalentBondKindForm::Undetermined)]
    #[case::lit(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_kind_form_constructors(
        #[case] actual: NoncovalentBondKindForm,
        #[case] expected: NoncovalentBondKindForm,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hydrogen(NoncovalentBondKind::HydrogenBond, NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond))]
    #[case::ionic(NoncovalentBondKind::Ionic, NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic))]
    fn test_noncovalent_bond_kind_form_from(
        #[case] kind: NoncovalentBondKind,
        #[case] expected: NoncovalentBondKindForm,
    ) {
        assert_eq!(NoncovalentBondKindForm::from(kind), expected);
    }

    #[rstest]
    #[case::undetermined(NoncovalentBondKindForm::Undetermined)]
    #[case::lit(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_kind_form_normalize_identity(#[case] input: NoncovalentBondKindForm) {
        assert_eq!(input.clone().normalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HalogenBond), Some(NoncovalentBondKind::HalogenBond))]
    #[case::undetermined(NoncovalentBondKindForm::Undetermined, None)]
    fn test_noncovalent_bond_kind_form_as_lit(
        #[case] form: NoncovalentBondKindForm,
        #[case] expected: Option<NoncovalentBondKind>,
    ) {
        assert_eq!(form.as_lit(), expected);
        assert_eq!(form.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(NoncovalentBondKindForm::Undetermined, true)]
    #[case::lit(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), false)]
    fn test_noncovalent_bond_kind_form_is_undetermined(
        #[case] form: NoncovalentBondKindForm,
        #[case] expected: bool,
    ) {
        assert_eq!(form.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(NoncovalentBondKindForm::Undetermined, NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), Some(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond)))]
    #[case::lit_und(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Undetermined, Some(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond)))]
    #[case::und_und(NoncovalentBondKindForm::Undetermined, NoncovalentBondKindForm::Undetermined, Some(NoncovalentBondKindForm::Undetermined))]
    #[case::lit_lit_eq(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), Some(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond)))]
    #[case::lit_lit_neq(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic), None)]
    fn test_noncovalent_bond_kind_form_meet(
        #[case] a: NoncovalentBondKindForm,
        #[case] b: NoncovalentBondKindForm,
        #[case] expected: Option<NoncovalentBondKindForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(NoncovalentBondKindForm::Undetermined, NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Undetermined)]
    #[case::und_und(NoncovalentBondKindForm::Undetermined, NoncovalentBondKindForm::Undetermined, NoncovalentBondKindForm::Undetermined)]
    #[case::lit_lit_eq(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond))]
    #[case::lit_lit_neq(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic), NoncovalentBondKindForm::Undetermined)]
    fn test_noncovalent_bond_kind_form_join(
        #[case] a: NoncovalentBondKindForm,
        #[case] b: NoncovalentBondKindForm,
        #[case] expected: NoncovalentBondKindForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_lit(NoncovalentBondKindForm::Undetermined, NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), true)]
    #[case::undetermined_undetermined(NoncovalentBondKindForm::Undetermined, NoncovalentBondKindForm::Undetermined, true)]
    #[case::lit_undetermined(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Undetermined, false)]
    #[case::lit_lit_match(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), true)]
    #[case::lit_lit_mismatch(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic), false)]
    fn test_noncovalent_bond_kind_form_matches(
        #[case] pattern: NoncovalentBondKindForm,
        #[case] target: NoncovalentBondKindForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(
        NoncovalentBondForm::default(),
        NoncovalentBondForm::default(),
        Some(NoncovalentBondForm::default())
    )]
    #[case::kind_mismatch(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HalogenBond),
        None
    )]
    #[case::kind_narrows_from_undetermined(
        NoncovalentBondForm::default(),
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        Some(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))
    )]
    fn test_noncovalent_bond_form_meet(
        #[case] a: NoncovalentBondForm,
        #[case] b: NoncovalentBondForm,
        #[case] expected: Option<NoncovalentBondForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::kind_mismatch_widens_to_default(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HalogenBond),
        NoncovalentBondForm::default()
    )]
    fn test_noncovalent_bond_form_join(
        #[case] a: NoncovalentBondForm,
        #[case] b: NoncovalentBondForm,
        #[case] expected: NoncovalentBondForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }
}
