//! Noncovalent bond AST.

use std::borrow::Cow;

use umol_graph_core::{ParticipantPosition, RelationData};
use umol_graph_ir_macros::{Canonicalize, Lattice};

use super::constraint::{NoncovalentBondConstraintAst, NoncovalentBondConstraintsAst};
use super::error::{Contradiction, NoJoin};
use super::traits::{AsLit, Canonicalize, Lattice};

/// Noncovalent bond: two-atom non-bonded interaction tagged by an
/// interaction kind. No bond order, no charge or spin — these do not apply
/// to noncovalent interactions.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Canonicalize, Lattice)]
pub struct NoncovalentBondAst {
    pub kind: NoncovalentBondKindForm,
    pub constraints: NoncovalentBondConstraintsAst,
}

/// Attribute update for a noncovalent bond. The kind is optional, and an
/// undetermined constraint removes its key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NoncovalentBondUpdate {
    pub kind: Option<NoncovalentBondKindForm>,
    pub constraints: NoncovalentBondConstraintsAst,
}

impl From<&str> for NoncovalentBondAst {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid noncovalent bond string")
    }
}

impl RelationData for NoncovalentBondAst {
    /// The kind is not position-indexed — reordering the two participants leaves it unchanged.
    fn on_permutation(&mut self, _order: &[ParticipantPosition]) {}

    fn is_permutation_invariant(&self) -> bool {
        true
    }
}

impl NoncovalentBondAst {
    pub fn new(kind: NoncovalentBondKindForm) -> Self {
        Self {
            kind,
            constraints: NoncovalentBondConstraintsAst::new(),
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
    /// kind (last-wins per `NoncovalentBondConstraintsAst::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<NoncovalentBondConstraintAst>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `NoncovalentBondConstraintsAst::set`).
    /// Does not clear existing constraints; use `bond.constraints.clear()`
    /// or direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<NoncovalentBondConstraintAst>,
    {
        for c in constraints {
            self.constraints.set(c.into());
        }
        self
    }

    /// No-op on value fields: `NoncovalentBondAst` has no value-bearing field
    /// besides `kind`, which is essential and never filled. Constraints are
    /// preserved. Provided for API symmetry.
    pub fn into_ground(self) -> Self {
        self
    }

    /// Apply an attribute update, leaving an omitted kind and constraint keys unchanged.
    pub fn update(&self, update: &NoncovalentBondUpdate) -> NoncovalentBondAst {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        NoncovalentBondAst {
            kind: update.kind.clone().unwrap_or_else(|| self.kind.clone()),
            constraints,
        }
    }

    /// Derive the minimal canonical attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> NoncovalentBondUpdate {
        let mut constraints = NoncovalentBondConstraintsAst::new();
        for new in other.constraints.iter() {
            if self
                .constraints
                .get(new.key())
                .is_none_or(|old| !old.canonical_eq(new))
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
            kind: (!self.kind.canonical_eq(&other.kind)).then(|| other.kind.clone()),
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

impl Canonicalize for NoncovalentBondKindForm {
    /// Both variants are already canonical — nothing folds.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(self)
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
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
    #[case::new(NoncovalentBondAst::new(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond)),
        NoncovalentBondAst { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraintsAst::new() })]
    #[case::from_kind(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondAst { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraintsAst::new() })]
    fn test_noncovalent_bond_form_new(
        #[case] actual: NoncovalentBondAst,
        #[case] expected: NoncovalentBondAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_kind_primitive(
        NoncovalentBondAst::default().with_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondAst { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraintsAst::new() })]
    #[case::with_kind_ast(
        NoncovalentBondAst::default().with_kind(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond)),
        NoncovalentBondAst { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraintsAst::new() })]
    #[case::with_constraints_empty(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)
            .with_constraints(empty::<NoncovalentBondConstraintAst>()),
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))]
    #[case::with_constraint(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)
            .with_constraint(NoncovalentBondConstraintAst::intramolecular(true)),
        NoncovalentBondAst { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)) })]
    #[case::with_constraints_populated(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)
            .with_constraints([NoncovalentBondConstraintAst::intramolecular(false)]),
        NoncovalentBondAst { kind: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(false)) })]
    fn test_noncovalent_bond_form_with_methods(
        #[case] actual: NoncovalentBondAst,
        #[case] expected: NoncovalentBondAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::default_(NoncovalentBondAst::default())]
    #[case::ground(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_form_into_ground(#[case] bond: NoncovalentBondAst) {
        assert_eq!(bond.clone().into_ground(), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondUpdate { kind: Some(NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic)), ..Default::default() },
        NoncovalentBondAst::from_kind(NoncovalentBondKind::Ionic))]
    #[case::kind_undetermined(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondUpdate { kind: Some(NoncovalentBondKindForm::Undetermined), ..Default::default() },
        NoncovalentBondAst::default())]
    #[case::constraint_set(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), ..Default::default() },
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(true)))]
    #[case::constraint_replace(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(true)),
        NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(false)), ..Default::default() },
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(false)))]
    #[case::constraint_remove(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(true)),
        NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Undetermined)), ..Default::default() },
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_form_update(
        #[case] bond: NoncovalentBondAst,
        #[case] update: NoncovalentBondUpdate,
        #[case] expected: NoncovalentBondAst,
    ) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(true)))]
    fn test_noncovalent_bond_form_update_identity(#[case] bond: NoncovalentBondAst) {
        assert_eq!(bond.update(&NoncovalentBondUpdate::default()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_and_constraint(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(true)),
        NoncovalentBondAst::default(),
        NoncovalentBondUpdate {
            kind: Some(NoncovalentBondKindForm::Undetermined),
            constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::Intramolecular(BooleanForm::Undetermined)),
        },
    )]
    fn test_noncovalent_bond_form_difference_to(
        #[case] bond: NoncovalentBondAst,
        #[case] other: NoncovalentBondAst,
        #[case] expected: NoncovalentBondUpdate,
    ) {
        assert_eq!(bond.difference_to(&other), expected);
    }

    #[rstest]
    #[case::same(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_form_difference_to_identity(#[case] bond: NoncovalentBondAst) {
        assert_eq!(bond.difference_to(&bond), NoncovalentBondUpdate::default());
    }

    #[rstest]
    #[case::default_(NoncovalentBondAst::default(), false)]
    #[case::ground_lit(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    fn test_noncovalent_bond_form_is_ground(
        #[case] ast: NoncovalentBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::ground(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))]
    #[case::undetermined(NoncovalentBondAst::default())]
    fn test_noncovalent_bond_form_canonicalize_identity(#[case] input: NoncovalentBondAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_ground(NoncovalentBondAst::default(), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    #[case::same(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    #[case::different(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondAst::from_kind(NoncovalentBondKind::Ionic), false)]
    #[case::pattern_specific_target_undetermined(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondAst::default(), false)]
    fn test_noncovalent_bond_form_matches(
        #[case] pattern: NoncovalentBondAst,
        #[case] target: NoncovalentBondAst,
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
    fn test_noncovalent_bond_kind_form_canonicalize_identity(
        #[case] input: NoncovalentBondKindForm,
    ) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HalogenBond), Some(NoncovalentBondKind::HalogenBond))]
    #[case::undetermined(NoncovalentBondKindForm::Undetermined, None)]
    fn test_noncovalent_bond_kind_form_as_lit(
        #[case] ast: NoncovalentBondKindForm,
        #[case] expected: Option<NoncovalentBondKind>,
    ) {
        assert_eq!(ast.as_lit(), expected);
        assert_eq!(ast.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(NoncovalentBondKindForm::Undetermined, true)]
    #[case::lit(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), false)]
    fn test_noncovalent_bond_kind_form_is_undetermined(
        #[case] ast: NoncovalentBondKindForm,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_undetermined(), expected);
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
        NoncovalentBondAst::default(),
        NoncovalentBondAst::default(),
        Some(NoncovalentBondAst::default())
    )]
    #[case::kind_mismatch(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HalogenBond),
        None
    )]
    #[case::kind_narrows_from_undetermined(
        NoncovalentBondAst::default(),
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        Some(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))
    )]
    fn test_noncovalent_bond_form_meet(
        #[case] a: NoncovalentBondAst,
        #[case] b: NoncovalentBondAst,
        #[case] expected: Option<NoncovalentBondAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::kind_mismatch_widens_to_default(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HalogenBond),
        NoncovalentBondAst::default()
    )]
    fn test_noncovalent_bond_form_join(
        #[case] a: NoncovalentBondAst,
        #[case] b: NoncovalentBondAst,
        #[case] expected: NoncovalentBondAst,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }
}
