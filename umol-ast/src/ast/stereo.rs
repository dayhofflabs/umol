//! Stereochemistry AST: the configuration value and the operator-expression
//! tree over it.
//!
//! A configuration value is a dense coset index per stereo class — the
//! OpenSMILES arrangement number of that class's coset space (`umol-perm`).
//! `~` and `^` are group actions on the index; [`StereoConfigurationAst::simplify`]
//! folds closed operator-expressions against the coset algebra. The class
//! ([`StereoKind`]) is the interpretation context that the operators consume, so
//! it is passed to `simplify` rather than carried in the value.

use std::mem;

use umol_perm::{space, ClassKey, Permutation};

use super::constraint::{
    StereoAtomConstraint, StereoAtomConstraints, StereoBondConstraint, StereoBondConstraints,
};
use super::traits::{AsLit, Lattice};

/// A stereo class: the atom-centered coordination geometries and the bond
/// cis/trans class, all sharing the configuration machinery.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum StereoKind {
    Tetrahedral,
    CisTrans,
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
}

impl StereoKind {
    /// The `~` involution for this kind: swap a trans ligand pair — the axial
    /// pair for trigonal-bipyramidal/octahedral, a diagonal for square-planar.
    /// Tetrahedral has no trans pair, so it swaps the first two ligands;
    /// cis/trans swaps the two configurations.
    fn involution(self) -> Permutation {
        match self {
            StereoKind::Tetrahedral => Permutation::from_image(4, &[1, 0, 2, 3]),
            StereoKind::CisTrans => Permutation::from_image(2, &[1, 0]),
            StereoKind::SquarePlanar => Permutation::from_image(4, &[2, 1, 0, 3]),
            StereoKind::TrigonalBipyramidal => Permutation::from_image(5, &[4, 1, 2, 3, 0]),
            StereoKind::Octahedral => Permutation::from_image(6, &[5, 1, 2, 3, 4, 0]),
        }
    }

    /// Act on coset index `index` by operation `operation`, through the class's
    /// coset algebra.
    fn act(self, index: u32, operation: Permutation) -> u32 {
        let class = match self {
            StereoKind::Tetrahedral => ClassKey::Tetrahedral,
            StereoKind::CisTrans => ClassKey::CisTrans,
            StereoKind::SquarePlanar => ClassKey::SquarePlanar,
            StereoKind::TrigonalBipyramidal => ClassKey::TrigonalBipyramidal,
            StereoKind::Octahedral => ClassKey::Octahedral,
        };
        space(class).reindex(index, operation)
    }
}

/// A stereo configuration: undetermined (pattern wildcard), explicitly not a
/// stereocenter, or a stereo index.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoConfigurationAst {
    #[default]
    Undetermined,
    NotStereo,
    Stereo(StereoIndexAst),
}

impl StereoConfigurationAst {
    /// Reduce to canonical form under the class context: lift trivial `Expr`
    /// wrappers and fold closed operator-expressions via the coset algebra.
    /// Free-variable expressions are left as-is.
    pub fn simplify(self, kind: StereoKind) -> Self {
        match self {
            Self::Stereo(index) => Self::Stereo(index.simplify(kind)),
            other => other,
        }
    }
}

impl AsLit for StereoConfigurationAst {
    type Lit = u32;

    /// The coset index of a resolved stereocenter. `NotStereo` and
    /// `Undetermined` (and unfolded expressions) carry no coset index → `None`.
    fn as_lit(&self) -> Option<u32> {
        match self {
            Self::Stereo(index) => index.as_lit(),
            Self::Undetermined | Self::NotStereo => None,
        }
    }
}

impl Lattice for StereoConfigurationAst {
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    fn is_ground(&self) -> bool {
        match self {
            Self::NotStereo => true,
            Self::Stereo(index) => index.is_ground(),
            Self::Undetermined => false,
        }
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::NotStereo, Self::NotStereo) => Some(Self::NotStereo),
            (Self::Stereo(a), Self::Stereo(b)) => a.meet(b).map(Self::Stereo),
            (Self::NotStereo, Self::Stereo(_)) | (Self::Stereo(_), Self::NotStereo) => None,
        }
    }

    fn join(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::NotStereo, Self::NotStereo) => Self::NotStereo,
            (Self::Stereo(a), Self::Stereo(b)) => Self::Stereo(a.join(b)),
            (Self::NotStereo, Self::Stereo(_)) | (Self::Stereo(_), Self::NotStereo) => {
                Self::Undetermined
            }
        }
    }

    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::NotStereo, Self::NotStereo) => true,
            (Self::Stereo(p), Self::Stereo(t)) => p.matches(t),
            (Self::NotStereo, Self::Stereo(_)) | (Self::Stereo(_), Self::NotStereo) => false,
        }
    }
}

/// The stereo index: undetermined, a literal coset index, or an
/// operator-expression over indices.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoIndexAst {
    #[default]
    Undetermined,
    Lit(u32),
    Expr(Box<Expr>),
}

impl StereoIndexAst {
    /// Simplify the inner expression and lift a folded `Expr(Lit)` to `Lit`.
    pub fn simplify(self, kind: StereoKind) -> Self {
        match self {
            Self::Expr(e) => match e.simplify(kind) {
                Expr::Lit(index) => Self::Lit(index),
                other => Self::Expr(Box::new(other)),
            },
            other => other,
        }
    }
}

impl AsLit for StereoIndexAst {
    type Lit = u32;

    fn as_lit(&self) -> Option<u32> {
        match self {
            Self::Lit(index) => Some(*index),
            Self::Undetermined | Self::Expr(_) => None,
        }
    }
}

impl Lattice for StereoIndexAst {
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    fn is_ground(&self) -> bool {
        matches!(self, Self::Lit(_))
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::Lit(a), Self::Lit(b)) => (a == b).then_some(Self::Lit(*a)),
            (Self::Expr(e), Self::Expr(f)) => (e == f).then(|| Self::Expr(e.clone())),
            (Self::Expr(_), _) | (_, Self::Expr(_)) => None,
        }
    }

    fn join(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::Lit(a), Self::Lit(b)) if a == b => Self::Lit(*a),
            (Self::Expr(e), Self::Expr(f)) if e == f => Self::Expr(e.clone()),
            _ => Self::Undetermined,
        }
    }

    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::Lit(p), Self::Lit(t)) => p == t,
            (Self::Expr(_), _) | (_, Self::Expr(_)) => false,
        }
    }
}

/// An operator-expression over coset indices: a literal, a bound variable, the
/// `~` involution, the generic `^` action by a permutation, or a (deferred)
/// literal set / variable domain.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Expr {
    Lit(u32),
    Var(String),
    SwapOp(Box<Expr>),
    ApplyOp(Box<Expr>, Permutation),
    LitSet(Vec<u32>),
    VarDomain(String, Vec<u32>),
}

impl Expr {
    /// Recursively simplify, folding closed operator nodes via the coset
    /// algebra: `~~e → e` (involution); `~(Lit k) → Lit(k · involution)`;
    /// `(Lit k) ^ g → Lit(k · g)`. Free `Var`, `Set`, and `VarDomain` are inert.
    pub fn simplify(self, kind: StereoKind) -> Self {
        match self {
            Expr::Lit(_) | Expr::Var(_) | Expr::LitSet(_) | Expr::VarDomain(..) => self,
            Expr::SwapOp(inner) => match inner.simplify(kind) {
                Expr::SwapOp(inner2) => *inner2,
                Expr::Lit(index) => Expr::Lit(kind.act(index, kind.involution())),
                other => Expr::SwapOp(Box::new(other)),
            },
            Expr::ApplyOp(inner, operation) => match inner.simplify(kind) {
                Expr::Lit(index) => Expr::Lit(kind.act(index, operation)),
                other => Expr::ApplyOp(Box::new(other), operation),
            },
        }
    }
}

impl From<u32> for StereoIndexAst {
    fn from(index: u32) -> Self {
        Self::Lit(index)
    }
}

impl From<StereoIndexAst> for StereoConfigurationAst {
    fn from(index: StereoIndexAst) -> Self {
        Self::Stereo(index)
    }
}

impl From<u32> for StereoConfigurationAst {
    fn from(index: u32) -> Self {
        Self::Stereo(StereoIndexAst::Lit(index))
    }
}

/// An atom-centered stereo element: its geometry class, configuration, and
/// per-site constraints. Predicate only — the site atom and ligands live in the
/// relation overlay, not here.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct StereoAtomAst {
    pub kind: StereoKind,
    pub configuration: StereoConfigurationAst,
    pub constraints: StereoAtomConstraints,
}

impl StereoAtomAst {
    pub fn new(kind: StereoKind) -> Self {
        Self {
            kind,
            configuration: StereoConfigurationAst::Undetermined,
            constraints: StereoAtomConstraints::new(),
        }
    }

    pub fn with_configuration(mut self, configuration: impl Into<StereoConfigurationAst>) -> Self {
        self.configuration = configuration.into();
        self
    }

    /// Add each constraint from the iterator. Vacuous today since
    /// `StereoAtomConstraint` is uninhabited; the signature locks in the extend
    /// semantic for when inhabited variants land.
    pub fn with_constraints<I>(self, _constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<StereoAtomConstraint>,
    {
        self
    }

    /// No-op: an unspecified configuration has no zero default (it is either a
    /// concrete coset or `NotStereo`). Provided for API symmetry.
    pub fn into_ground(self) -> Self {
        self
    }

    /// Equivalent to `into_ground()`. `StereoAtomAst` has no constraint defaults.
    pub fn into_zeroed(self) -> Self {
        self.into_ground()
    }

    /// Fold the configuration's closed operator-expressions (under `kind`) and
    /// simplify each constraint's value in place.
    pub fn simplify_values(&mut self) {
        self.configuration = mem::take(&mut self.configuration).simplify(self.kind);
        self.constraints.simplify_each();
    }
}

impl Lattice for StereoAtomAst {
    fn is_undetermined(&self) -> bool {
        self.configuration.is_undetermined() && self.constraints.is_undetermined()
    }

    fn is_ground(&self) -> bool {
        self.configuration.is_ground() && self.constraints.is_ground()
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        if self.kind != other.kind {
            return None;
        }
        Some(Self {
            kind: self.kind,
            configuration: self.configuration.meet(&other.configuration)?,
            constraints: self.constraints.meet(&other.constraints)?,
        })
    }

    /// Same-site payloads always share a kind, so the kind-mismatch arm is a
    /// defensive identity rather than a meaningful join.
    fn join(&self, other: &Self) -> Self {
        if self.kind != other.kind {
            return self.clone();
        }
        Self {
            kind: self.kind,
            configuration: self.configuration.join(&other.configuration),
            constraints: self.constraints.join(&other.constraints),
        }
    }

    fn matches(&self, target: &Self) -> bool {
        self.kind == target.kind
            && self.configuration.matches(&target.configuration)
            && self.constraints.matches(&target.constraints)
    }
}

/// A bond-centered stereo element (cis/trans). Predicate only — the site bond
/// and ligands live in the relation overlay.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct StereoBondAst {
    pub kind: StereoKind,
    pub configuration: StereoConfigurationAst,
    pub constraints: StereoBondConstraints,
}

impl StereoBondAst {
    pub fn new(kind: StereoKind) -> Self {
        Self {
            kind,
            configuration: StereoConfigurationAst::Undetermined,
            constraints: StereoBondConstraints::new(),
        }
    }

    pub fn with_configuration(mut self, configuration: impl Into<StereoConfigurationAst>) -> Self {
        self.configuration = configuration.into();
        self
    }

    /// Add each constraint from the iterator. Vacuous today since
    /// `StereoBondConstraint` is uninhabited; the signature locks in the extend
    /// semantic for when inhabited variants land.
    pub fn with_constraints<I>(self, _constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<StereoBondConstraint>,
    {
        self
    }

    /// No-op: an unspecified configuration has no zero default. Provided for API
    /// symmetry.
    pub fn into_ground(self) -> Self {
        self
    }

    /// Equivalent to `into_ground()`. `StereoBondAst` has no constraint defaults.
    pub fn into_zeroed(self) -> Self {
        self.into_ground()
    }

    /// Fold the configuration's closed operator-expressions (under `kind`) and
    /// simplify each constraint's value in place.
    pub fn simplify_values(&mut self) {
        self.configuration = mem::take(&mut self.configuration).simplify(self.kind);
        self.constraints.simplify_each();
    }
}

impl Lattice for StereoBondAst {
    fn is_undetermined(&self) -> bool {
        self.configuration.is_undetermined() && self.constraints.is_undetermined()
    }

    fn is_ground(&self) -> bool {
        self.configuration.is_ground() && self.constraints.is_ground()
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        if self.kind != other.kind {
            return None;
        }
        Some(Self {
            kind: self.kind,
            configuration: self.configuration.meet(&other.configuration)?,
            constraints: self.constraints.meet(&other.constraints)?,
        })
    }

    /// Same-site payloads always share a kind, so the kind-mismatch arm is a
    /// defensive identity rather than a meaningful join.
    fn join(&self, other: &Self) -> Self {
        if self.kind != other.kind {
            return self.clone();
        }
        Self {
            kind: self.kind,
            configuration: self.configuration.join(&other.configuration),
            constraints: self.constraints.join(&other.constraints),
        }
    }

    fn matches(&self, target: &Self) -> bool {
        self.kind == target.kind
            && self.configuration.matches(&target.configuration)
            && self.constraints.matches(&target.constraints)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoKind::Tetrahedral, Expr::Lit(1), Expr::Lit(1))]
    #[case::var(StereoKind::Tetrahedral, Expr::Var("o".into()), Expr::Var("o".into()))]
    #[case::swap_lit_even(StereoKind::Tetrahedral, Expr::SwapOp(Box::new(Expr::Lit(0))), Expr::Lit(1))]
    #[case::swap_lit_odd(StereoKind::Tetrahedral, Expr::SwapOp(Box::new(Expr::Lit(1))), Expr::Lit(0))]
    #[case::double_swap_lit(StereoKind::Tetrahedral, Expr::SwapOp(Box::new(Expr::SwapOp(Box::new(Expr::Lit(1))))), Expr::Lit(1))]
    #[case::double_swap_var(StereoKind::Tetrahedral, Expr::SwapOp(Box::new(Expr::SwapOp(Box::new(Expr::Var("o".into()))))), Expr::Var("o".into()))]
    #[case::swap_var_stays(StereoKind::Tetrahedral, Expr::SwapOp(Box::new(Expr::Var("o".into()))), Expr::SwapOp(Box::new(Expr::Var("o".into()))))]
    #[case::apply_lit(StereoKind::Tetrahedral, Expr::ApplyOp(Box::new(Expr::Lit(0)), Permutation::from_image(4, &[1, 0, 2, 3])), Expr::Lit(1))]
    #[case::apply_identity(StereoKind::Tetrahedral, Expr::ApplyOp(Box::new(Expr::Lit(1)), Permutation::from_image(4, &[0, 1, 2, 3])), Expr::Lit(1))]
    #[case::apply_var_stays(StereoKind::Tetrahedral, Expr::ApplyOp(Box::new(Expr::Var("o".into())), Permutation::from_image(4, &[1, 0, 2, 3])), Expr::ApplyOp(Box::new(Expr::Var("o".into())), Permutation::from_image(4, &[1, 0, 2, 3])))]
    #[case::cistrans_swap(StereoKind::CisTrans, Expr::SwapOp(Box::new(Expr::Lit(0))), Expr::Lit(1))]
    #[case::sp_swap_u_fixed(StereoKind::SquarePlanar, Expr::SwapOp(Box::new(Expr::Lit(0))), Expr::Lit(0))]
    #[case::sp_swap_four_z(StereoKind::SquarePlanar, Expr::SwapOp(Box::new(Expr::Lit(1))), Expr::Lit(2))]
    #[case::tb_swap_axial(StereoKind::TrigonalBipyramidal, Expr::SwapOp(Box::new(Expr::Lit(0))), Expr::Lit(1))]
    #[case::tb_swap_other(StereoKind::TrigonalBipyramidal, Expr::SwapOp(Box::new(Expr::Lit(2))), Expr::Lit(17))]
    #[case::oh_swap_axial(StereoKind::Octahedral, Expr::SwapOp(Box::new(Expr::Lit(0))), Expr::Lit(1))]
    #[case::oh_swap_other(StereoKind::Octahedral, Expr::SwapOp(Box::new(Expr::Lit(2))), Expr::Lit(21))]
    fn test_expr_simplify(#[case] kind: StereoKind, #[case] input: Expr, #[case] expected: Expr) {
        assert_eq!(input.simplify(kind), expected);
    }

    #[rstest]
    #[case::swap_var(Expr::SwapOp(Box::new(Expr::Var("o".into()))))]
    #[case::apply_var(Expr::ApplyOp(Box::new(Expr::Var("o".into())), Permutation::from_image(4, &[1, 0, 2, 3])))]
    #[case::double_swap_lit(Expr::SwapOp(Box::new(Expr::SwapOp(Box::new(Expr::Lit(0))))))]
    fn test_expr_simplify_idempotent(#[case] input: Expr) {
        let once = input.simplify(StereoKind::Tetrahedral);
        let twice = once.clone().simplify(StereoKind::Tetrahedral);
        assert_eq!(once, twice);
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, 0)]
    #[case::cis_trans(StereoKind::CisTrans, 0)]
    #[case::square_planar(StereoKind::SquarePlanar, 1)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, 2)]
    #[case::octahedral(StereoKind::Octahedral, 2)]
    fn test_expr_simplify_involution(#[case] kind: StereoKind, #[case] index: u32) {
        let double = Expr::SwapOp(Box::new(Expr::SwapOp(Box::new(Expr::Lit(index)))));
        assert_eq!(double.simplify(kind), Expr::Lit(index));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoIndexAst::Lit(1), StereoIndexAst::Lit(1))]
    #[case::undetermined(StereoIndexAst::Undetermined, StereoIndexAst::Undetermined)]
    #[case::expr_lit_lifts(StereoIndexAst::Expr(Box::new(Expr::Lit(2))), StereoIndexAst::Lit(2))]
    #[case::expr_swap_lifts(StereoIndexAst::Expr(Box::new(Expr::SwapOp(Box::new(Expr::Lit(0))))), StereoIndexAst::Lit(1))]
    #[case::expr_var_stays(StereoIndexAst::Expr(Box::new(Expr::Var("o".into()))), StereoIndexAst::Expr(Box::new(Expr::Var("o".into()))))]
    fn test_stereo_index_ast_simplify(#[case] input: StereoIndexAst, #[case] expected: StereoIndexAst) {
        assert_eq!(input.simplify(StereoKind::Tetrahedral), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoConfigurationAst::Undetermined, StereoConfigurationAst::Undetermined)]
    #[case::not_stereo(StereoConfigurationAst::NotStereo, StereoConfigurationAst::NotStereo)]
    #[case::stereo_lit(StereoConfigurationAst::Stereo(StereoIndexAst::Lit(1)), StereoConfigurationAst::Stereo(StereoIndexAst::Lit(1)))]
    #[case::stereo_expr_lifts(
        StereoConfigurationAst::Stereo(StereoIndexAst::Expr(Box::new(Expr::SwapOp(Box::new(Expr::Lit(0)))))),
        StereoConfigurationAst::Stereo(StereoIndexAst::Lit(1)),
    )]
    fn test_stereo_configuration_ast_simplify(
        #[case] input: StereoConfigurationAst,
        #[case] expected: StereoConfigurationAst,
    ) {
        assert_eq!(input.simplify(StereoKind::Tetrahedral), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::u32(StereoConfigurationAst::from(2u32), StereoConfigurationAst::Stereo(StereoIndexAst::Lit(2)))]
    #[case::index(StereoConfigurationAst::from(StereoIndexAst::Lit(3)), StereoConfigurationAst::Stereo(StereoIndexAst::Lit(3)))]
    fn test_stereo_configuration_ast_from(
        #[case] actual: StereoConfigurationAst,
        #[case] expected: StereoConfigurationAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stereo_lit(StereoConfigurationAst::Stereo(StereoIndexAst::Lit(2)), Some(2))]
    #[case::not_stereo(StereoConfigurationAst::NotStereo, None)]
    #[case::undetermined(StereoConfigurationAst::Undetermined, None)]
    #[case::stereo_undetermined(StereoConfigurationAst::Stereo(StereoIndexAst::Undetermined), None)]
    fn test_stereo_configuration_ast_as_lit(
        #[case] config: StereoConfigurationAst,
        #[case] expected: Option<u32>,
    ) {
        assert_eq!(config.as_lit(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoConfigurationAst::Undetermined, false)]
    #[case::not_stereo(StereoConfigurationAst::NotStereo, true)]
    #[case::stereo_lit(StereoConfigurationAst::from(1u32), true)]
    #[case::stereo_undetermined(StereoConfigurationAst::Stereo(StereoIndexAst::Undetermined), false)]
    fn test_stereo_configuration_ast_is_ground(
        #[case] config: StereoConfigurationAst,
        #[case] expected: bool,
    ) {
        assert_eq!(config.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_narrows(StereoConfigurationAst::Undetermined, StereoConfigurationAst::from(1u32), Some(StereoConfigurationAst::from(1u32)))]
    #[case::lit_same(StereoConfigurationAst::from(1u32), StereoConfigurationAst::from(1u32), Some(StereoConfigurationAst::from(1u32)))]
    #[case::lit_conflict(StereoConfigurationAst::from(1u32), StereoConfigurationAst::from(2u32), None)]
    #[case::not_stereo_same(StereoConfigurationAst::NotStereo, StereoConfigurationAst::NotStereo, Some(StereoConfigurationAst::NotStereo))]
    #[case::not_stereo_vs_stereo(StereoConfigurationAst::NotStereo, StereoConfigurationAst::from(0u32), None)]
    fn test_stereo_configuration_ast_meet(
        #[case] a: StereoConfigurationAst,
        #[case] b: StereoConfigurationAst,
        #[case] expected: Option<StereoConfigurationAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_matches_any(StereoConfigurationAst::Undetermined, StereoConfigurationAst::from(1u32), true)]
    #[case::specific_vs_undetermined(StereoConfigurationAst::from(1u32), StereoConfigurationAst::Undetermined, false)]
    #[case::lit_match(StereoConfigurationAst::from(1u32), StereoConfigurationAst::from(1u32), true)]
    #[case::lit_mismatch(StereoConfigurationAst::from(1u32), StereoConfigurationAst::from(2u32), false)]
    #[case::not_stereo_match(StereoConfigurationAst::NotStereo, StereoConfigurationAst::NotStereo, true)]
    #[case::not_stereo_vs_stereo(StereoConfigurationAst::NotStereo, StereoConfigurationAst::from(0u32), false)]
    fn test_stereo_configuration_ast_matches(
        #[case] pattern: StereoConfigurationAst,
        #[case] target: StereoConfigurationAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoIndexAst::Lit(2), Some(2))]
    #[case::undetermined(StereoIndexAst::Undetermined, None)]
    #[case::expr(StereoIndexAst::Expr(Box::new(Expr::Var("o".into()))), None)]
    fn test_stereo_index_ast_as_lit(#[case] index: StereoIndexAst, #[case] expected: Option<u32>) {
        assert_eq!(index.as_lit(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_narrows(StereoIndexAst::Undetermined, StereoIndexAst::Lit(1), Some(StereoIndexAst::Lit(1)))]
    #[case::lit_same(StereoIndexAst::Lit(1), StereoIndexAst::Lit(1), Some(StereoIndexAst::Lit(1)))]
    #[case::lit_conflict(StereoIndexAst::Lit(1), StereoIndexAst::Lit(2), None)]
    #[case::expr_vs_lit(StereoIndexAst::Expr(Box::new(Expr::Var("o".into()))), StereoIndexAst::Lit(1), None)]
    fn test_stereo_index_ast_meet(
        #[case] a: StereoIndexAst,
        #[case] b: StereoIndexAst,
        #[case] expected: Option<StereoIndexAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    fn test_stereo_atom_ast_new() {
        assert_eq!(
            StereoAtomAst::new(StereoKind::Tetrahedral),
            StereoAtomAst {
                kind: StereoKind::Tetrahedral,
                configuration: StereoConfigurationAst::Undetermined,
                constraints: StereoAtomConstraints::new(),
            }
        );
    }

    #[rstest]
    fn test_stereo_atom_ast_with_configuration() {
        assert_eq!(
            StereoAtomAst::new(StereoKind::Tetrahedral).with_configuration(1u32),
            StereoAtomAst {
                kind: StereoKind::Tetrahedral,
                configuration: StereoConfigurationAst::Stereo(StereoIndexAst::Lit(1)),
                constraints: StereoAtomConstraints::new(),
            }
        );
    }

    #[rstest]
    fn test_stereo_atom_ast_simplify_values() {
        let mut atom = StereoAtomAst::new(StereoKind::Tetrahedral).with_configuration(
            StereoConfigurationAst::Stereo(StereoIndexAst::Expr(Box::new(Expr::SwapOp(Box::new(
                Expr::Lit(0),
            ))))),
        );
        atom.simplify_values();
        assert_eq!(
            atom.configuration,
            StereoConfigurationAst::Stereo(StereoIndexAst::Lit(1))
        );
    }

    #[rstest]
    #[case::undetermined(StereoAtomAst::new(StereoKind::Tetrahedral), false)]
    #[case::ground(StereoAtomAst::new(StereoKind::Tetrahedral).with_configuration(1u32), true)]
    fn test_stereo_atom_ast_is_ground(#[case] atom: StereoAtomAst, #[case] expected: bool) {
        assert_eq!(atom.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_narrows(StereoAtomAst::new(StereoKind::Tetrahedral), StereoAtomAst::new(StereoKind::Tetrahedral).with_configuration(1u32), Some(StereoAtomAst::new(StereoKind::Tetrahedral).with_configuration(1u32)))]
    #[case::different_kind(StereoAtomAst::new(StereoKind::Tetrahedral), StereoAtomAst::new(StereoKind::SquarePlanar), None)]
    #[case::config_conflict(StereoAtomAst::new(StereoKind::Tetrahedral).with_configuration(0u32), StereoAtomAst::new(StereoKind::Tetrahedral).with_configuration(1u32), None)]
    fn test_stereo_atom_ast_meet(
        #[case] a: StereoAtomAst,
        #[case] b: StereoAtomAst,
        #[case] expected: Option<StereoAtomAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_match(StereoAtomAst::new(StereoKind::Tetrahedral), StereoAtomAst::new(StereoKind::Tetrahedral).with_configuration(1u32), true)]
    #[case::different_kind(StereoAtomAst::new(StereoKind::Tetrahedral).with_configuration(1u32), StereoAtomAst::new(StereoKind::SquarePlanar).with_configuration(1u32), false)]
    fn test_stereo_atom_ast_matches(
        #[case] pattern: StereoAtomAst,
        #[case] target: StereoAtomAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    fn test_stereo_bond_ast_new() {
        assert_eq!(
            StereoBondAst::new(StereoKind::CisTrans),
            StereoBondAst {
                kind: StereoKind::CisTrans,
                configuration: StereoConfigurationAst::Undetermined,
                constraints: StereoBondConstraints::new(),
            }
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_narrows(StereoBondAst::new(StereoKind::CisTrans), StereoBondAst::new(StereoKind::CisTrans).with_configuration(1u32), Some(StereoBondAst::new(StereoKind::CisTrans).with_configuration(1u32)))]
    #[case::config_conflict(StereoBondAst::new(StereoKind::CisTrans).with_configuration(0u32), StereoBondAst::new(StereoKind::CisTrans).with_configuration(1u32), None)]
    fn test_stereo_bond_ast_meet(
        #[case] a: StereoBondAst,
        #[case] b: StereoBondAst,
        #[case] expected: Option<StereoBondAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }
}
