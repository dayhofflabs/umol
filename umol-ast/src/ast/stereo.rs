//! Stereochemistry AST: the configuration value and the operator-expression
//! tree over it.
//!
//! A configuration value is a dense coset index per stereo class — the
//! OpenSMILES arrangement number of that class's coset space (`umol-perm`).
//! `~` and `^` are group actions on the index; [`StereoConfigurationAst::simplify`]
//! folds closed operator-expressions against the coset algebra. The class
//! ([`StereoKind`]) is the interpretation context that the operators consume, so
//! it is passed to `simplify` rather than carried in the value.

use umol_perm::{space, ClassKey, Permutation};

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
                Expr::SwapOp(grandchild) => *grandchild,
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
}
