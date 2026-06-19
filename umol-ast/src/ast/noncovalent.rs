//! Noncovalent bond AST.

use std::borrow::Cow;

use umol_ast_macros::{Canonicalize, Lattice};

use super::constraint::{NoncovalentBondConstraint, NoncovalentBondConstraints};
use super::error::Contradiction;
use super::traits::{AsLit, Canonicalize, Lattice};

/// Noncovalent bond: two-atom non-bonded interaction tagged by an
/// interaction kind. No bond order, no charge or spin — these do not apply
/// to noncovalent interactions.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Canonicalize, Lattice)]
pub struct NoncovalentBondAst {
    pub kind: NoncovalentBondKindAst,
    pub constraints: NoncovalentBondConstraints,
}

impl NoncovalentBondAst {
    pub fn new(kind: NoncovalentBondKindAst) -> Self {
        Self {
            kind,
            constraints: NoncovalentBondConstraints::new(),
        }
    }

    pub fn from_kind(kind: NoncovalentBondKind) -> Self {
        Self::new(NoncovalentBondKindAst::Lit(kind))
    }

    pub fn with_kind(mut self, kind: impl Into<NoncovalentBondKindAst>) -> Self {
        self.kind = kind.into();
        self
    }

    /// Add each constraint from the iterator. Vacuous today since
    /// `NoncovalentBondConstraint` is uninhabited; the signature locks in
    /// the extend semantic so that callers don't develop replace-style
    /// expectations before inhabited variants land.
    pub fn with_constraints<I>(self, _constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<NoncovalentBondConstraint>,
    {
        self
    }

    /// No-op: `NoncovalentBondAst` has no value-bearing fields besides
    /// `kind`, which is essential and never filled. Provided for API
    /// symmetry.
    pub fn into_ground(self) -> Self {
        self
    }

    /// Equivalent to `into_ground()`. `NoncovalentBondAst` has no constraint
    /// defaults.
    pub fn into_zeroed(self) -> Self {
        self.into_ground()
    }
}

/// Noncovalent interaction kind.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NoncovalentBondKindAst {
    #[default]
    Undetermined,
    Lit(NoncovalentBondKind),
}

impl NoncovalentBondKindAst {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn lit(kind: NoncovalentBondKind) -> Self {
        Self::Lit(kind)
    }
}

impl From<NoncovalentBondKind> for NoncovalentBondKindAst {
    fn from(kind: NoncovalentBondKind) -> Self {
        Self::Lit(kind)
    }
}

impl Canonicalize for NoncovalentBondKindAst {
    /// Both variants are already canonical — nothing folds.
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(self)
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        Ok(Cow::Borrowed(self))
    }
}

impl AsLit for NoncovalentBondKindAst {
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

impl Lattice for NoncovalentBondKindAst {
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

    fn join(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::Lit(a), Self::Lit(b)) if a == b => Self::Lit(*a),
            _ => Self::Undetermined,
        }
    }

    /// `target` refines `self`: `self.meet(target) == canonical(target)`.
    fn matches(&self, target: &Self) -> bool {
        match (self.meet(target), target.canonical()) {
            (Some(meet), Ok(target)) => meet == *target,
            _ => false,
        }
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

    #[rustfmt::skip]
    #[rstest]
    #[case::new(NoncovalentBondAst::new(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond)),
        NoncovalentBondAst { kind: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraints::new() })]
    #[case::from_kind(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondAst { kind: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraints::new() })]
    fn test_noncovalent_bond_ast_new(
        #[case] actual: NoncovalentBondAst,
        #[case] expected: NoncovalentBondAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_kind_primitive(
        NoncovalentBondAst::default().with_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondAst { kind: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraints::new() })]
    #[case::with_kind_ast(
        NoncovalentBondAst::default().with_kind(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond)),
        NoncovalentBondAst { kind: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
            constraints: NoncovalentBondConstraints::new() })]
    #[case::with_constraints_empty(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)
            .with_constraints(empty::<NoncovalentBondConstraint>()),
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_ast_with_methods(
        #[case] actual: NoncovalentBondAst,
        #[case] expected: NoncovalentBondAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::default_(NoncovalentBondAst::default())]
    #[case::ground(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_ast_into_ground(#[case] bond: NoncovalentBondAst) {
        assert_eq!(bond.clone().into_ground(), bond);
    }

    #[rstest]
    #[case::default_(NoncovalentBondAst::default())]
    #[case::ground(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_ast_into_zeroed(#[case] bond: NoncovalentBondAst) {
        assert_eq!(bond.clone().into_zeroed(), bond.into_ground());
    }

    #[rstest]
    #[case::default_(NoncovalentBondAst::default(), false)]
    #[case::ground_lit(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    fn test_noncovalent_bond_ast_is_ground(
        #[case] ast: NoncovalentBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::ground(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))]
    #[case::undetermined(NoncovalentBondAst::default())]
    fn test_noncovalent_bond_ast_canonicalize_identity(#[case] input: NoncovalentBondAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_ground(NoncovalentBondAst::default(), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    #[case::same(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    #[case::different(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondAst::from_kind(NoncovalentBondKind::Ionic), false)]
    #[case::pattern_specific_target_undetermined(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondAst::default(), false)]
    fn test_noncovalent_bond_ast_matches(
        #[case] pattern: NoncovalentBondAst,
        #[case] target: NoncovalentBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(NoncovalentBondKindAst::Undetermined, NoncovalentBondKindAst::Undetermined)]
    #[case::lit(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_kind_ast_constructors(
        #[case] actual: NoncovalentBondKindAst,
        #[case] expected: NoncovalentBondKindAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hydrogen(NoncovalentBondKind::HydrogenBond, NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond))]
    #[case::ionic(NoncovalentBondKind::Ionic, NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic))]
    fn test_noncovalent_bond_kind_ast_from(
        #[case] kind: NoncovalentBondKind,
        #[case] expected: NoncovalentBondKindAst,
    ) {
        assert_eq!(NoncovalentBondKindAst::from(kind), expected);
    }

    #[rstest]
    #[case::undetermined(NoncovalentBondKindAst::Undetermined)]
    #[case::lit(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond))]
    fn test_noncovalent_bond_kind_ast_canonicalize_identity(#[case] input: NoncovalentBondKindAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HalogenBond), Some(NoncovalentBondKind::HalogenBond))]
    #[case::undetermined(NoncovalentBondKindAst::Undetermined, None)]
    fn test_noncovalent_bond_kind_ast_as_lit(
        #[case] ast: NoncovalentBondKindAst,
        #[case] expected: Option<NoncovalentBondKind>,
    ) {
        assert_eq!(ast.as_lit(), expected);
        assert_eq!(ast.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_match(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKind::HydrogenBond, true)]
    #[case::lit_mismatch(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKind::Ionic, false)]
    #[case::undetermined(NoncovalentBondKindAst::Undetermined, NoncovalentBondKind::HydrogenBond, false)]
    fn test_noncovalent_bond_kind_ast_as_lit_matches(
        #[case] ast: NoncovalentBondKindAst,
        #[case] value: NoncovalentBondKind,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.as_lit_matches(value), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(NoncovalentBondKindAst::Undetermined, true)]
    #[case::lit(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), false)]
    fn test_noncovalent_bond_kind_ast_is_undetermined(
        #[case] ast: NoncovalentBondKindAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(NoncovalentBondKindAst::Undetermined, NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), Some(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond)))]
    #[case::lit_und(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Undetermined, Some(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond)))]
    #[case::und_und(NoncovalentBondKindAst::Undetermined, NoncovalentBondKindAst::Undetermined, Some(NoncovalentBondKindAst::Undetermined))]
    #[case::lit_lit_eq(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), Some(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond)))]
    #[case::lit_lit_neq(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic), None)]
    fn test_noncovalent_bond_kind_ast_meet(
        #[case] a: NoncovalentBondKindAst,
        #[case] b: NoncovalentBondKindAst,
        #[case] expected: Option<NoncovalentBondKindAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(NoncovalentBondKindAst::Undetermined, NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Undetermined)]
    #[case::und_und(NoncovalentBondKindAst::Undetermined, NoncovalentBondKindAst::Undetermined, NoncovalentBondKindAst::Undetermined)]
    #[case::lit_lit_eq(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond))]
    #[case::lit_lit_neq(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic), NoncovalentBondKindAst::Undetermined)]
    fn test_noncovalent_bond_kind_ast_join(
        #[case] a: NoncovalentBondKindAst,
        #[case] b: NoncovalentBondKindAst,
        #[case] expected: NoncovalentBondKindAst,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_lit(NoncovalentBondKindAst::Undetermined, NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), true)]
    #[case::undetermined_undetermined(NoncovalentBondKindAst::Undetermined, NoncovalentBondKindAst::Undetermined, true)]
    #[case::lit_undetermined(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Undetermined, false)]
    #[case::lit_lit_match(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), true)]
    #[case::lit_lit_mismatch(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic), false)]
    fn test_noncovalent_bond_kind_ast_matches(
        #[case] pattern: NoncovalentBondKindAst,
        #[case] target: NoncovalentBondKindAst,
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
    fn test_noncovalent_bond_ast_meet(
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
    fn test_noncovalent_bond_ast_join(
        #[case] a: NoncovalentBondAst,
        #[case] b: NoncovalentBondAst,
        #[case] expected: NoncovalentBondAst,
    ) {
        assert_eq!(a.join(&b), expected);
    }
}
