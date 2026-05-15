//! Noncovalent bond AST.

use super::constraint::{NoncovalentBondConstraint, NoncovalentBondConstraints};
use super::traits::Lattice;

/// Noncovalent bond: two-atom non-bonded interaction tagged by an
/// interaction kind. No bond order, no charge or spin — these do not apply
/// to noncovalent interactions.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
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

    pub fn matches(&self, target: &NoncovalentBondAst) -> bool {
        self.kind.matches(&target.kind)
    }

    /// Simplify every constraint's inner value in place. `kind` carries no
    /// `ValueAst`, so it is unchanged.
    pub fn simplify_values(&mut self) {
        self.constraints.simplify_each();
    }
}

impl Lattice for NoncovalentBondAst {
    fn is_undetermined(&self) -> bool {
        matches!(self.kind, NoncovalentBondKindAst::Undetermined)
            && self.constraints.is_undetermined()
    }

    fn is_ground(&self) -> bool {
        self.kind.is_ground() && self.constraints.is_ground()
    }

    /// `kind` is treated as a structural anchor (equality-only): mismatched
    /// kinds → `None`.
    fn meet(&self, other: &Self) -> Option<Self> {
        let kind = match (&self.kind, &other.kind) {
            (NoncovalentBondKindAst::Undetermined, k) | (k, NoncovalentBondKindAst::Undetermined) => {
                k.clone()
            }
            (a, b) if a == b => a.clone(),
            _ => return None,
        };
        Some(Self {
            kind,
            constraints: self.constraints.meet(&other.constraints)?,
        })
    }

    /// `kind` is equality-only: mismatched kinds widen the entity to
    /// `Default::default()` (entity-level undetermined).
    fn join(&self, other: &Self) -> Self {
        match (&self.kind, &other.kind) {
            (NoncovalentBondKindAst::Undetermined, _) | (_, NoncovalentBondKindAst::Undetermined) => {
                Self::default()
            }
            (a, b) if a == b => Self {
                kind: a.clone(),
                constraints: self.constraints.join(&other.constraints),
            },
            _ => Self::default(),
        }
    }
}

/// Noncovalent interaction kind expressions. Mirrors `ElementAst`:
/// wildcard, literal, set, bind, ref.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum NoncovalentBondKindAst {
    #[default]
    Undetermined,
    Lit(NoncovalentBondKind),
    Set(Vec<NoncovalentBondKind>),
    Bind {
        id: String,
        set: Vec<NoncovalentBondKind>,
    },
    Ref(String),
}

impl From<NoncovalentBondKind> for NoncovalentBondKindAst {
    fn from(kind: NoncovalentBondKind) -> Self {
        Self::Lit(kind)
    }
}

impl NoncovalentBondKindAst {
    pub fn is_ground(&self) -> bool {
        matches!(self, Self::Lit(_))
    }

    pub fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::Ref(_), _) | (_, Self::Ref(_)) => false,
            (Self::Lit(p), Self::Lit(t)) => p == t,
            (Self::Lit(p), Self::Set(ts) | Self::Bind { set: ts, .. }) => ts.iter().all(|t| t == p),
            (Self::Set(ps) | Self::Bind { set: ps, .. }, Self::Lit(t)) => ps.contains(t),
            (
                Self::Set(ps) | Self::Bind { set: ps, .. },
                Self::Set(ts) | Self::Bind { set: ts, .. },
            ) => ts.iter().all(|t| ps.contains(t)),
        }
    }
}

/// Fundamental kind of a noncovalent interaction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
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
        NoncovalentBondAst::default().with_kind(NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond])),
        NoncovalentBondAst { kind: NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond]),
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

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_ground(NoncovalentBondAst::default(), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    #[case::same(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), true)]
    #[case::different(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondAst::from_kind(NoncovalentBondKind::Ionic), false)]
    fn test_noncovalent_bond_ast_matches(
        #[case] pattern: NoncovalentBondAst,
        #[case] target: NoncovalentBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    fn test_noncovalent_bond_ast_simplify_values() {
        let mut bond = NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond);
        let original = bond.clone();
        bond.simplify_values();
        assert_eq!(bond, original);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), true)]
    #[case::undetermined(NoncovalentBondKindAst::Undetermined, false)]
    #[case::set(NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond, NoncovalentBondKind::Ionic]), false)]
    #[case::bind(NoncovalentBondKindAst::Bind { id: "k".into(), set: vec![NoncovalentBondKind::HydrogenBond] }, false)]
    #[case::reference(NoncovalentBondKindAst::Ref("k".into()), false)]
    fn test_noncovalent_bond_kind_ast_is_ground(
        #[case] ast: NoncovalentBondKindAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_lit(NoncovalentBondKindAst::Undetermined, NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), true)]
    #[case::undetermined_undetermined(NoncovalentBondKindAst::Undetermined, NoncovalentBondKindAst::Undetermined, true)]
    #[case::lit_undetermined(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Undetermined, false)]
    #[case::lit_lit_match(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), true)]
    #[case::lit_lit_mismatch(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic), false)]
    #[case::set_lit_in(NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond, NoncovalentBondKind::Ionic]), NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic), true)]
    #[case::set_lit_out(NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond]), NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic), false)]
    #[case::set_set_subset(NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond, NoncovalentBondKind::Ionic, NoncovalentBondKind::VanDerWaals]),
        NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond, NoncovalentBondKind::Ionic]), true)]
    #[case::set_set_superset(NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond]),
        NoncovalentBondKindAst::Set(vec![NoncovalentBondKind::HydrogenBond, NoncovalentBondKind::Ionic]), false)]
    #[case::bind_lit_match(NoncovalentBondKindAst::Bind { id: "k".into(), set: vec![NoncovalentBondKind::HydrogenBond] },
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), true)]
    #[case::ref_lit(NoncovalentBondKindAst::Ref("k".into()), NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), false)]
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
        Some(NoncovalentBondAst::default()),
    )]
    #[case::kind_mismatch(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HalogenBond),
        None,
    )]
    #[case::kind_narrows_from_undetermined(
        NoncovalentBondAst::default(),
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        Some(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)),
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
        NoncovalentBondAst::default(),
    )]
    fn test_noncovalent_bond_ast_join(
        #[case] a: NoncovalentBondAst,
        #[case] b: NoncovalentBondAst,
        #[case] expected: NoncovalentBondAst,
    ) {
        assert_eq!(a.join(&b), expected);
    }
}
