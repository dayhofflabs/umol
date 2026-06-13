//! Stereochemistry AST: the configuration value and the operator-expression
//! tree over it.
//!
//! A configuration value is a dense coset index per stereo kind, corresponds to OpenSMILES
//! numbering for SP, TB, and OH.
//! `~` and `^` are group actions on the index; [`StereoConfigurationAst::simplify`]
//! folds closed operator-expressions against the coset algebra.

use std::mem;

use umol_perm::{space, ClassKey, Permutation};

use super::constraint::{
    StereoAtomConstraint, StereoAtomConstraints, StereoBondConstraint, StereoBondConstraints,
};
use super::traits::{AsLit, Lattice};

/// Stereo kind: the atom-centered coordination geometries and the bond-centered cis/trans kind.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum StereoKind {
    Tetrahedral,
    CisTrans,
    Axial,
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
}

impl StereoKind {
    /// The `umol-perm` class key for this stereo kind.
    pub fn class_key(self) -> ClassKey {
        match self {
            StereoKind::Tetrahedral => ClassKey::Tetrahedral,
            StereoKind::CisTrans => ClassKey::CisTrans,
            StereoKind::Axial => ClassKey::Axial,
            StereoKind::SquarePlanar => ClassKey::SquarePlanar,
            StereoKind::TrigonalBipyramidal => ClassKey::TrigonalBipyramidal,
            StereoKind::Octahedral => ClassKey::Octahedral,
        }
    }

    /// Number of ligand positions in this stereo kind.
    pub fn degree(self) -> usize {
        space(self.class_key()).degree()
    }

    /// Number of cosets/configurations in this stereo kind.
    pub fn count(self) -> usize {
        space(self.class_key()).count()
    }

    /// Whether this stereo kind can encode local handedness.
    pub fn is_chiral_class(self) -> bool {
        space(self.class_key()).is_chiral()
    }

    /// Kind-specific `~` involution. Chiral kinds borrow the orientation-reversing
    /// generator from umol-perm; achiral kinds use a chosen ligand swap (no improper
    /// generator to borrow — theirs is the identity):
    /// - cis/trans: swap the two configurations
    /// - square-planar: swap the diagonal ligand pair
    fn involution(self) -> Permutation {
        let coset_space = space(self.class_key());
        if coset_space.is_chiral() {
            coset_space.improper()
        } else {
            match self {
                StereoKind::CisTrans => Permutation::from_image(4, &[1, 0, 2, 3]),
                StereoKind::SquarePlanar => Permutation::from_image(4, &[2, 1, 0, 3]),
                _ => unreachable!("only achiral kinds reach the chosen-swap branch"),
            }
        }
    }

    /// Act on coset index `index` by `perm`, through the class's coset algebra.
    fn act(self, index: u32, perm: Permutation) -> u32 {
        space(self.class_key()).reindex(index, perm)
    }
}

/// Topicity of two ligand positions of a stereo carrier (derived ground value).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, strum::VariantArray)]
pub enum Topicity {
    Homotopic,
    Enantiotopic,
    Diastereotopic,
}

/// Stereogenicity classification of a stereo carrier (derived ground value).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, strum::VariantArray)]
pub enum Stereogenicity {
    Symmetric,
    Prochiral,
    Stereogenic,
}

/// Stereo configuration AST.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoConfigurationAst {
    #[default]
    Undetermined,
    NotStereo,
    Stereo(StereoCosetAst),
}

impl StereoConfigurationAst {
    pub fn stereo(v: impl Into<StereoCosetAst>) -> Self {
        Self::Stereo(v.into())
    }

    pub fn is_stereo(&self) -> bool {
        matches!(self, Self::Stereo(_))
    }

    /// Reduce to canonical form under the kind context: lift trivial `StereoExpr`
    /// wrappers and fold closed operator-expressions via the coset algebra.
    pub fn simplify(self, kind: StereoKind) -> Self {
        match self {
            Self::Stereo(index) => Self::Stereo(index.simplify(kind)),
            other => other,
        }
    }

    /// Matches literal coset index `value` under `kind`.
    pub fn matches_value(&self, value: u32, kind: StereoKind) -> bool {
        match self {
            Self::Stereo(v) => v.matches_value(value, kind),
            Self::NotStereo => false,
            Self::Undetermined => true,
        }
    }
}

impl AsLit for StereoConfigurationAst {
    type Lit = u32;

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

/// Dense coset index AST. 1/2 for TH, CT, follows OpenSMILES numbering for SP, TB, OH.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoCosetAst {
    #[default]
    Undetermined,
    Lit(u32),
    Expr(Box<StereoExpr>),
}

impl StereoCosetAst {
    /// Wrap operator-expression as coset index.
    pub fn expr(expr: StereoExpr) -> Self {
        Self::Expr(Box::new(expr))
    }

    /// Simplify the inner expression under `kind`'s coset algebra.
    pub fn simplify(self, kind: StereoKind) -> Self {
        match self {
            Self::Expr(e) => match e.simplify(kind) {
                StereoExpr::Lit(index) => Self::Lit(index),
                other => Self::Expr(Box::new(other)),
            },
            other => other,
        }
    }

    /// Apply a ligand-order permutation to this coset under `kind`.
    pub fn apply_permutation(&self, kind: StereoKind, perm: Permutation) -> Self {
        match self {
            Self::Undetermined => Self::Undetermined,
            Self::Lit(index) => Self::Lit(kind.act(*index, perm)),
            Self::Expr(expr) => {
                Self::Expr(Box::new(StereoExpr::ApplyOp(expr.clone(), perm))).simplify(kind)
            }
        }
    }

    /// Matches literal coset index `value` under `kind`.
    pub fn matches_value(&self, value: u32, kind: StereoKind) -> bool {
        match self {
            Self::Undetermined => true,
            Self::Lit(v) => *v == value,
            Self::Expr(e) => e.matches_value(value, kind),
        }
    }
}

impl AsLit for StereoCosetAst {
    type Lit = u32;

    fn as_lit(&self) -> Option<u32> {
        match self {
            Self::Lit(index) => Some(*index),
            Self::Undetermined | Self::Expr(_) => None,
        }
    }
}

impl Lattice for StereoCosetAst {
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

/// Operator expression over coset algebra for `kind`. Includes ~ and ^ operators, and literals.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoExpr {
    Lit(u32),
    Var(String),
    SwapOp(Box<StereoExpr>),
    MirrorOp(Box<StereoExpr>),
    ApplyOp(Box<StereoExpr>, Permutation),
    LitSet(Vec<u32>),
    VarDomain(String, Vec<u32>),
}

impl StereoExpr {
    /// `~inner` — the involution operator applied to `inner`.
    pub fn swap(inner: Self) -> Self {
        Self::SwapOp(Box::new(inner))
    }

    /// `'inner` — the improper (mirror) operator: the enantiomer of `inner`.
    pub fn mirror(inner: Self) -> Self {
        Self::MirrorOp(Box::new(inner))
    }

    /// `inner ^ perm` — the generic group action of `perm` on `inner`.
    pub fn apply(inner: Self, perm: Permutation) -> Self {
        Self::ApplyOp(Box::new(inner), perm)
    }

    /// Recursively simplify, folding closed operator nodes via the coset
    /// algebra: `~~e → e` (involution); `~(Lit k) → Lit(k · involution)`;
    /// `(Lit k) ^ g → Lit(k · g)`.
    pub fn simplify(self, kind: StereoKind) -> Self {
        match self {
            StereoExpr::Lit(_)
            | StereoExpr::Var(_)
            | StereoExpr::LitSet(_)
            | StereoExpr::VarDomain(..) => self,
            StereoExpr::SwapOp(inner) => match inner.simplify(kind) {
                StereoExpr::SwapOp(inner2) => *inner2,
                StereoExpr::Lit(index) => StereoExpr::Lit(kind.act(index, kind.involution())),
                other => StereoExpr::SwapOp(Box::new(other)),
            },
            StereoExpr::MirrorOp(inner) => match inner.simplify(kind) {
                StereoExpr::MirrorOp(inner2) => *inner2,
                StereoExpr::Lit(index) => {
                    StereoExpr::Lit(space(kind.class_key()).enantiomer(index))
                }
                other => StereoExpr::MirrorOp(Box::new(other)),
            },
            StereoExpr::ApplyOp(inner, perm) => match inner.simplify(kind) {
                StereoExpr::Lit(index) => StereoExpr::Lit(kind.act(index, perm)),
                other => StereoExpr::ApplyOp(Box::new(other), perm),
            },
        }
    }

    /// Matches literal coset index `value` under `kind`'s coset algebra.
    pub fn matches_value(&self, value: u32, kind: StereoKind) -> bool {
        match self {
            StereoExpr::Lit(n) => *n == value,
            StereoExpr::Var(_) => true,
            StereoExpr::LitSet(set) | StereoExpr::VarDomain(_, set) => set.contains(&value),
            StereoExpr::SwapOp(inner) => {
                inner.matches_value(kind.act(value, kind.involution()), kind)
            }
            StereoExpr::MirrorOp(inner) => {
                inner.matches_value(space(kind.class_key()).enantiomer(value), kind)
            }
            StereoExpr::ApplyOp(inner, perm) => {
                inner.matches_value(kind.act(value, perm.inverse()), kind)
            }
        }
    }
}

impl From<StereoCosetAst> for StereoConfigurationAst {
    fn from(index: StereoCosetAst) -> Self {
        Self::Stereo(index)
    }
}

impl From<u32> for StereoConfigurationAst {
    fn from(index: u32) -> Self {
        Self::Stereo(StereoCosetAst::Lit(index))
    }
}

impl From<Vec<u32>> for StereoConfigurationAst {
    fn from(values: Vec<u32>) -> Self {
        Self::Stereo(StereoCosetAst::Expr(Box::new(StereoExpr::LitSet(values))))
    }
}

impl From<u32> for StereoCosetAst {
    fn from(index: u32) -> Self {
        Self::Lit(index)
    }
}

impl From<Vec<u32>> for StereoCosetAst {
    fn from(values: Vec<u32>) -> Self {
        Self::Expr(Box::new(StereoExpr::LitSet(values)))
    }
}

/// StereoAtomAst and StereoBondAst generator
macro_rules! stereo_element {
    (
        $(#[doc = $doc:literal])+
        $name:ident, $constraints:ident, $constraint:ident
    ) => {
        $(#[doc = $doc])+
        #[derive(Clone, Debug, PartialEq, Eq, Hash)]
        pub struct $name {
            pub kind: StereoKind,
            pub coset: StereoCosetAst,
            pub constraints: $constraints,
        }

        impl $name {
            pub fn new(kind: StereoKind, coset: impl Into<StereoCosetAst>) -> Self {
                Self {
                    kind,
                    coset: coset.into(),
                    constraints: $constraints::new(),
                }
            }

            /// Add a single constraint.
            pub fn with_constraint(self, _constraint: impl Into<$constraint>) -> Self {
                self
            }

            /// Add each constraint from iterator.
            pub fn with_constraints<I>(mut self, constraints: I) -> Self
            where
                I: IntoIterator,
                I::Item: Into<$constraint>,
            {
                self.constraints.extend(constraints.into_iter().map(Into::into));
                self
            }

            /// No-op. A stereo element is always stereogenic, so its coset has no
            /// zero default — an unspecified coset has no canonical ground term.
            /// The element is ground iff its coset is ground.
            pub fn into_ground(self) -> Self {
                self
            }

            /// Equivalent to `into_ground()`; there are no constraint defaults.
            pub fn into_zeroed(self) -> Self {
                self.into_ground()
            }

            /// Fold the coset's closed operator-expressions (under `kind`) and
            /// simplify each constraint's value in place.
            pub fn simplify_values(&mut self) {
                self.coset = mem::take(&mut self.coset).simplify(self.kind);
                self.constraints.simplify_each();
            }
        }

        impl Lattice for $name {
            fn is_undetermined(&self) -> bool {
                self.coset.is_undetermined() && self.constraints.is_undetermined()
            }

            fn is_ground(&self) -> bool {
                self.coset.is_ground() && self.constraints.is_ground()
            }

            fn meet(&self, other: &Self) -> Option<Self> {
                if self.kind != other.kind {
                    return None;
                }
                Some(Self {
                    kind: self.kind,
                    coset: self.coset.meet(&other.coset)?,
                    constraints: self.constraints.meet(&other.constraints)?,
                })
            }

            /// Per-kind lattices. Cross-kind join is a precondition violation.
            fn join(&self, other: &Self) -> Self {
                debug_assert_eq!(
                    self.kind, other.kind,
                    "stereo elements join only within a kind"
                );
                Self {
                    kind: self.kind,
                    coset: self.coset.join(&other.coset),
                    constraints: self.constraints.join(&other.constraints),
                }
            }

            fn matches(&self, target: &Self) -> bool {
                self.kind == target.kind
                    && self.coset.matches(&target.coset)
                    && self.constraints.matches(&target.constraints)
            }
        }
    };
}

stereo_element! {
    /// StereoAtomAst with geometry class, configuration, and per-site constraints.
    StereoAtomAst, StereoAtomConstraints, StereoAtomConstraint
}

stereo_element! {
    /// StereoBondAst with cis/trans configuration and per-site constraints.
    StereoBondAst, StereoBondConstraints, StereoBondConstraint
}

/// Default implementation for StereoAtomAst.
impl Default for StereoAtomAst {
    fn default() -> Self {
        Self {
            kind: StereoKind::Tetrahedral,
            coset: StereoCosetAst::Undetermined,
            constraints: StereoAtomConstraints::new(),
        }
    }
}

/// Default implementation for StereoBondAst.
impl Default for StereoBondAst {
    fn default() -> Self {
        Self {
            kind: StereoKind::CisTrans,
            coset: StereoCosetAst::Undetermined,
            constraints: StereoBondConstraints::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, ClassKey::Tetrahedral)]
    #[case::cis_trans(StereoKind::CisTrans, ClassKey::CisTrans)]
    #[case::axial(StereoKind::Axial, ClassKey::Axial)]
    #[case::square_planar(StereoKind::SquarePlanar, ClassKey::SquarePlanar)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, ClassKey::TrigonalBipyramidal)]
    #[case::octahedral(StereoKind::Octahedral, ClassKey::Octahedral)]
    fn test_stereo_kind_class_key(#[case] kind: StereoKind, #[case] expected: ClassKey) {
        assert_eq!(kind.class_key(), expected);
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, 4)]
    #[case::cis_trans(StereoKind::CisTrans, 4)]
    #[case::axial(StereoKind::Axial, 4)]
    #[case::square_planar(StereoKind::SquarePlanar, 4)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, 5)]
    #[case::octahedral(StereoKind::Octahedral, 6)]
    fn test_stereo_kind_degree(#[case] kind: StereoKind, #[case] expected: usize) {
        assert_eq!(kind.degree(), expected);
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, 2)]
    #[case::cis_trans(StereoKind::CisTrans, 2)]
    #[case::axial(StereoKind::Axial, 2)]
    #[case::square_planar(StereoKind::SquarePlanar, 3)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, 20)]
    #[case::octahedral(StereoKind::Octahedral, 30)]
    fn test_stereo_kind_count(#[case] kind: StereoKind, #[case] expected: usize) {
        assert_eq!(kind.count(), expected);
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, true)]
    #[case::cis_trans(StereoKind::CisTrans, false)]
    #[case::axial(StereoKind::Axial, true)]
    #[case::square_planar(StereoKind::SquarePlanar, false)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, true)]
    #[case::octahedral(StereoKind::Octahedral, true)]
    fn test_stereo_kind_is_chiral_class(#[case] kind: StereoKind, #[case] expected: bool) {
        assert_eq!(kind.is_chiral_class(), expected);
    }

    #[rstest]
    fn test_stereo_expr_swap() {
        assert_eq!(
            StereoExpr::swap(StereoExpr::Lit(0)).simplify(StereoKind::Tetrahedral),
            StereoExpr::Lit(1),
        );
    }

    #[rstest]
    fn test_stereo_expr_apply() {
        assert_eq!(
            StereoExpr::apply(
                StereoExpr::Lit(0),
                Permutation::from_image(4, &[1, 0, 2, 3])
            )
            .simplify(StereoKind::Tetrahedral),
            StereoExpr::Lit(1),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoKind::Tetrahedral, StereoExpr::Lit(1), StereoExpr::Lit(1))]
    #[case::var(StereoKind::Tetrahedral, StereoExpr::Var("o".into()), StereoExpr::Var("o".into()))]
    #[case::swap_lit_even(StereoKind::Tetrahedral, StereoExpr::SwapOp(Box::new(StereoExpr::Lit(0))), StereoExpr::Lit(1))]
    #[case::swap_lit_odd(StereoKind::Tetrahedral, StereoExpr::SwapOp(Box::new(StereoExpr::Lit(1))), StereoExpr::Lit(0))]
    #[case::double_swap_lit(StereoKind::Tetrahedral, StereoExpr::SwapOp(Box::new(StereoExpr::SwapOp(Box::new(StereoExpr::Lit(1))))), StereoExpr::Lit(1))]
    #[case::double_swap_var(StereoKind::Tetrahedral, StereoExpr::SwapOp(Box::new(StereoExpr::SwapOp(Box::new(StereoExpr::Var("o".into()))))), StereoExpr::Var("o".into()))]
    #[case::swap_var_stays(StereoKind::Tetrahedral, StereoExpr::SwapOp(Box::new(StereoExpr::Var("o".into()))), StereoExpr::SwapOp(Box::new(StereoExpr::Var("o".into()))))]
    #[case::apply_lit(StereoKind::Tetrahedral, StereoExpr::ApplyOp(Box::new(StereoExpr::Lit(0)), Permutation::from_image(4, &[1, 0, 2, 3])), StereoExpr::Lit(1))]
    #[case::apply_identity(StereoKind::Tetrahedral, StereoExpr::ApplyOp(Box::new(StereoExpr::Lit(1)), Permutation::from_image(4, &[0, 1, 2, 3])), StereoExpr::Lit(1))]
    #[case::apply_var_stays(StereoKind::Tetrahedral, StereoExpr::ApplyOp(Box::new(StereoExpr::Var("o".into())), Permutation::from_image(4, &[1, 0, 2, 3])), StereoExpr::ApplyOp(Box::new(StereoExpr::Var("o".into())), Permutation::from_image(4, &[1, 0, 2, 3])))]
    #[case::cistrans_swap(StereoKind::CisTrans, StereoExpr::SwapOp(Box::new(StereoExpr::Lit(0))), StereoExpr::Lit(1))]
    #[case::sp_swap_u_fixed(StereoKind::SquarePlanar, StereoExpr::SwapOp(Box::new(StereoExpr::Lit(0))), StereoExpr::Lit(0))]
    #[case::sp_swap_four_z(StereoKind::SquarePlanar, StereoExpr::SwapOp(Box::new(StereoExpr::Lit(1))), StereoExpr::Lit(2))]
    #[case::tb_swap_axial(StereoKind::TrigonalBipyramidal, StereoExpr::SwapOp(Box::new(StereoExpr::Lit(0))), StereoExpr::Lit(1))]
    #[case::tb_swap_other(StereoKind::TrigonalBipyramidal, StereoExpr::SwapOp(Box::new(StereoExpr::Lit(2))), StereoExpr::Lit(17))]
    #[case::oh_swap_axial(StereoKind::Octahedral, StereoExpr::SwapOp(Box::new(StereoExpr::Lit(0))), StereoExpr::Lit(1))]
    #[case::oh_swap_other(StereoKind::Octahedral, StereoExpr::SwapOp(Box::new(StereoExpr::Lit(2))), StereoExpr::Lit(21))]
    #[case::mirror_chiral(StereoKind::Tetrahedral, StereoExpr::MirrorOp(Box::new(StereoExpr::Lit(0))), StereoExpr::Lit(1))]
    #[case::mirror_achiral_noop(StereoKind::CisTrans, StereoExpr::MirrorOp(Box::new(StereoExpr::Lit(0))), StereoExpr::Lit(0))]
    #[case::double_mirror_lit(StereoKind::Tetrahedral, StereoExpr::MirrorOp(Box::new(StereoExpr::MirrorOp(Box::new(StereoExpr::Lit(1))))), StereoExpr::Lit(1))]
    #[case::mirror_var_stays(StereoKind::Tetrahedral, StereoExpr::MirrorOp(Box::new(StereoExpr::Var("o".into()))), StereoExpr::MirrorOp(Box::new(StereoExpr::Var("o".into()))))]
    fn test_stereo_expr_simplify(#[case] kind: StereoKind, #[case] input: StereoExpr, #[case] expected: StereoExpr) {
        assert_eq!(input.simplify(kind), expected);
    }

    #[rstest]
    #[case::swap_var(StereoExpr::SwapOp(Box::new(StereoExpr::Var("o".into()))))]
    #[case::apply_var(StereoExpr::ApplyOp(Box::new(StereoExpr::Var("o".into())), Permutation::from_image(4, &[1, 0, 2, 3])))]
    #[case::double_swap_lit(StereoExpr::SwapOp(Box::new(StereoExpr::SwapOp(Box::new(
        StereoExpr::Lit(0)
    )))))]
    fn test_stereo_expr_simplify_idempotent(#[case] input: StereoExpr) {
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
    fn test_stereo_expr_simplify_involution(#[case] kind: StereoKind, #[case] index: u32) {
        let double = StereoExpr::SwapOp(Box::new(StereoExpr::SwapOp(Box::new(StereoExpr::Lit(
            index,
        )))));
        assert_eq!(double.simplify(kind), StereoExpr::Lit(index));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoCosetAst::Lit(1), StereoCosetAst::Lit(1))]
    #[case::undetermined(StereoCosetAst::Undetermined, StereoCosetAst::Undetermined)]
    #[case::expr_lit_lifts(StereoCosetAst::Expr(Box::new(StereoExpr::Lit(2))), StereoCosetAst::Lit(2))]
    #[case::expr_swap_lifts(StereoCosetAst::Expr(Box::new(StereoExpr::SwapOp(Box::new(StereoExpr::Lit(0))))), StereoCosetAst::Lit(1))]
    #[case::expr_var_stays(StereoCosetAst::Expr(Box::new(StereoExpr::Var("o".into()))), StereoCosetAst::Expr(Box::new(StereoExpr::Var("o".into()))))]
    fn test_stereo_coset_ast_simplify(#[case] input: StereoCosetAst, #[case] expected: StereoCosetAst) {
        assert_eq!(input.simplify(StereoKind::Tetrahedral), expected);
    }

    #[rstest]
    fn test_stereo_coset_ast_expr() {
        assert_eq!(
            StereoCosetAst::expr(StereoExpr::swap(StereoExpr::Lit(0)))
                .simplify(StereoKind::Tetrahedral),
            StereoCosetAst::Lit(1),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoConfigurationAst::Undetermined, StereoConfigurationAst::Undetermined)]
    #[case::not_stereo(StereoConfigurationAst::NotStereo, StereoConfigurationAst::NotStereo)]
    #[case::stereo_lit(StereoConfigurationAst::Stereo(StereoCosetAst::Lit(1)), StereoConfigurationAst::Stereo(StereoCosetAst::Lit(1)))]
    #[case::stereo_expr_lifts(
        StereoConfigurationAst::Stereo(StereoCosetAst::Expr(Box::new(StereoExpr::SwapOp(Box::new(StereoExpr::Lit(0)))))),
        StereoConfigurationAst::Stereo(StereoCosetAst::Lit(1)),
    )]
    fn test_stereo_configuration_ast_simplify(
        #[case] input: StereoConfigurationAst,
        #[case] expected: StereoConfigurationAst,
    ) {
        assert_eq!(input.simplify(StereoKind::Tetrahedral), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::u32(StereoConfigurationAst::from(2u32), StereoConfigurationAst::Stereo(StereoCosetAst::Lit(2)))]
    #[case::index(StereoConfigurationAst::from(StereoCosetAst::Lit(3)), StereoConfigurationAst::Stereo(StereoCosetAst::Lit(3)))]
    fn test_stereo_configuration_ast_from(
        #[case] actual: StereoConfigurationAst,
        #[case] expected: StereoConfigurationAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stereo_lit(StereoConfigurationAst::Stereo(StereoCosetAst::Lit(2)), Some(2))]
    #[case::not_stereo(StereoConfigurationAst::NotStereo, None)]
    #[case::undetermined(StereoConfigurationAst::Undetermined, None)]
    #[case::stereo_undetermined(StereoConfigurationAst::Stereo(StereoCosetAst::Undetermined), None)]
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
    #[case::stereo_undetermined(StereoConfigurationAst::Stereo(StereoCosetAst::Undetermined), false)]
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
    #[case::lit(StereoCosetAst::Lit(2), Some(2))]
    #[case::undetermined(StereoCosetAst::Undetermined, None)]
    #[case::expr(StereoCosetAst::Expr(Box::new(StereoExpr::Var("o".into()))), None)]
    fn test_stereo_coset_ast_as_lit(#[case] index: StereoCosetAst, #[case] expected: Option<u32>) {
        assert_eq!(index.as_lit(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_narrows(StereoCosetAst::Undetermined, StereoCosetAst::Lit(1), Some(StereoCosetAst::Lit(1)))]
    #[case::lit_same(StereoCosetAst::Lit(1), StereoCosetAst::Lit(1), Some(StereoCosetAst::Lit(1)))]
    #[case::lit_conflict(StereoCosetAst::Lit(1), StereoCosetAst::Lit(2), None)]
    #[case::expr_vs_lit(StereoCosetAst::Expr(Box::new(StereoExpr::Var("o".into()))), StereoCosetAst::Lit(1), None)]
    fn test_stereo_coset_ast_meet(
        #[case] a: StereoCosetAst,
        #[case] b: StereoCosetAst,
        #[case] expected: Option<StereoCosetAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(StereoCosetAst::Undetermined, 1, StereoKind::Tetrahedral, true)]
    #[case::lit_match(StereoCosetAst::Lit(2), 2, StereoKind::Octahedral, true)]
    #[case::lit_miss(StereoCosetAst::Lit(2), 3, StereoKind::Octahedral, false)]
    #[case::var_wildcard(StereoCosetAst::expr(StereoExpr::Var("o".into())), 4, StereoKind::Octahedral, true)]
    #[case::lit_set_member(StereoCosetAst::expr(StereoExpr::LitSet(vec![1, 3])), 3, StereoKind::Octahedral, true)]
    #[case::lit_set_nonmember(StereoCosetAst::expr(StereoExpr::LitSet(vec![1, 3])), 2, StereoKind::Octahedral, false)]
    #[case::var_domain_member(StereoCosetAst::expr(StereoExpr::VarDomain("o".into(), vec![1, 3])), 1, StereoKind::Octahedral, true)]
    #[case::var_domain_nonmember(StereoCosetAst::expr(StereoExpr::VarDomain("o".into(), vec![1, 3])), 2, StereoKind::Octahedral, false)]
    #[case::swap_pulls_back(StereoCosetAst::expr(StereoExpr::swap(StereoExpr::Lit(0))), 1, StereoKind::Tetrahedral, true)]
    #[case::swap_pulls_back_miss(StereoCosetAst::expr(StereoExpr::swap(StereoExpr::Lit(0))), 0, StereoKind::Tetrahedral, false)]
    fn test_stereo_coset_ast_matches_value(
        #[case] coset: StereoCosetAst,
        #[case] value: u32,
        #[case] kind: StereoKind,
        #[case] expected: bool,
    ) {
        assert_eq!(coset.matches_value(value, kind), expected);
    }

    #[rstest]
    fn test_stereo_atom_ast_new() {
        let stereo_atom = StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined);
        assert_eq!(stereo_atom.kind, StereoKind::Tetrahedral);
        assert_eq!(stereo_atom.coset, StereoCosetAst::Undetermined);
        assert_eq!(stereo_atom.constraints, StereoAtomConstraints::new());
    }

    #[rstest]
    fn test_stereo_atom_ast_simplify_values() {
        let mut atom = StereoAtomAst::new(
            StereoKind::Tetrahedral,
            StereoCosetAst::Expr(Box::new(StereoExpr::SwapOp(Box::new(StereoExpr::Lit(0))))),
        );
        atom.simplify_values();
        assert_eq!(atom.coset, StereoCosetAst::Lit(1));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined), false)]
    #[case::ground(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), true)]
    fn test_stereo_atom_ast_is_ground(#[case] atom: StereoAtomAst, #[case] expected: bool) {
        assert_eq!(atom.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::open_coset(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined))]
    #[case::ground(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32))]
    fn test_stereo_atom_ast_into_ground(#[case] atom: StereoAtomAst) {
        assert_eq!(atom.clone().into_ground(), atom);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_narrows(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined), StereoAtomAst::new(StereoKind::Tetrahedral, 1u32),
        Some(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32)))]
    #[case::different_kind(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined),
        StereoAtomAst::new(StereoKind::SquarePlanar, StereoCosetAst::Undetermined), None)]
    #[case::config_conflict(StereoAtomAst::new(StereoKind::Tetrahedral, 0u32), StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), None)]
    fn test_stereo_atom_ast_meet(
        #[case] a: StereoAtomAst,
        #[case] b: StereoAtomAst,
        #[case] expected: Option<StereoAtomAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_coset(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), StereoAtomAst::new(StereoKind::Tetrahedral, 1u32))]
    #[case::distinct_cosets_widen(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), StereoAtomAst::new(StereoKind::Tetrahedral, 2u32), StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined))]
    fn test_stereo_atom_ast_join(#[case] a: StereoAtomAst, #[case] b: StereoAtomAst, #[case] expected: StereoAtomAst) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_match(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined), StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), true)]
    #[case::different_kind(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), StereoAtomAst::new(StereoKind::SquarePlanar, 1u32), false)]
    fn test_stereo_atom_ast_matches(
        #[case] pattern: StereoAtomAst,
        #[case] target: StereoAtomAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    fn test_stereo_bond_ast_new() {
        let stereo_bond = StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Undetermined);
        assert_eq!(stereo_bond.kind, StereoKind::CisTrans);
        assert_eq!(stereo_bond.coset, StereoCosetAst::Undetermined);
        assert_eq!(stereo_bond.constraints, StereoBondConstraints::new())
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_narrows(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Undetermined), StereoBondAst::new(StereoKind::CisTrans, 1u32),
        Some(StereoBondAst::new(StereoKind::CisTrans, 1u32)))]
    #[case::config_conflict(StereoBondAst::new(StereoKind::CisTrans, 0u32), StereoBondAst::new(StereoKind::CisTrans, 1u32), None)]
    fn test_stereo_bond_ast_meet(
        #[case] a: StereoBondAst,
        #[case] b: StereoBondAst,
        #[case] expected: Option<StereoBondAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }
}
