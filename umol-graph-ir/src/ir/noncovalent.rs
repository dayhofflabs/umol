//! Noncovalent bond form.

use std::borrow::Cow;

use umol_graph_core::{ParticipantPosition, RelationData};
use umol_graph_ir_macros::{Lattice, Normalize};

use super::constraint::{NoncovalentBondConstraintForm, NoncovalentBondConstraintsForm};
use super::error::{Contradiction, NoJoin};
use super::traits::{AsLit, Equiv, Lattice, Normalize};

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
    pub fn into_ground(self) -> Self {
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
    use std::iter::empty;

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ir::boolean::BooleanForm;

    #[rustfmt::skip]
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
            .with_constraints(empty::<NoncovalentBondConstraintForm>()),
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
    fn test_noncovalent_bond_form_into_ground(#[case] bond: NoncovalentBondForm) {
        assert_eq!(bond.clone().into_ground(), bond);
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
