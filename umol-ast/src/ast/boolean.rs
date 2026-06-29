//! Boolean AST value.

use std::borrow::Cow;

use super::error::Contradiction;
use super::traits::{AsLit, Canonicalize, Lattice};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BooleanAst {
    #[default]
    Undetermined,
    Lit(bool),
}

impl Canonicalize for BooleanAst {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(self)
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        Ok(Cow::Borrowed(self))
    }
}

impl AsLit for BooleanAst {
    type Lit = bool;

    fn as_lit(&self) -> Option<bool> {
        match self {
            Self::Lit(b) => Some(*b),
            Self::Undetermined => None,
        }
    }
}

impl Lattice for BooleanAst {
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    fn is_ground(&self) -> bool {
        matches!(self, Self::Lit(_))
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Undetermined, _) => Some(*other),
            (_, Self::Undetermined) => Some(*self),
            (Self::Lit(a), Self::Lit(b)) => (a == b).then_some(Self::Lit(*a)),
        }
    }

    fn join(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Lit(a), Self::Lit(b)) if a == b => Self::Lit(*a),
            _ => Self::Undetermined,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::lit_true(BooleanAst::Lit(true), Some(true))]
    #[case::lit_false(BooleanAst::Lit(false), Some(false))]
    #[case::undetermined(BooleanAst::Undetermined, None)]
    fn test_boolean_ast_as_lit(#[case] b: BooleanAst, #[case] expected: Option<bool>) {
        assert_eq!(b.as_lit(), expected);
    }

    #[rstest]
    #[case::undetermined(BooleanAst::Undetermined, true)]
    #[case::lit(BooleanAst::Lit(true), false)]
    fn test_boolean_ast_is_undetermined(#[case] b: BooleanAst, #[case] expected: bool) {
        assert_eq!(b.is_undetermined(), expected);
    }

    #[rstest]
    #[case::lit(BooleanAst::Lit(false), true)]
    #[case::undetermined(BooleanAst::Undetermined, false)]
    fn test_boolean_ast_is_ground(#[case] b: BooleanAst, #[case] expected: bool) {
        assert_eq!(b.is_ground(), expected);
    }

    #[rstest]
    #[case::top_left(
        BooleanAst::Undetermined,
        BooleanAst::Lit(true),
        Some(BooleanAst::Lit(true))
    )]
    #[case::top_right(
        BooleanAst::Lit(false),
        BooleanAst::Undetermined,
        Some(BooleanAst::Lit(false))
    )]
    #[case::same(
        BooleanAst::Lit(true),
        BooleanAst::Lit(true),
        Some(BooleanAst::Lit(true))
    )]
    #[case::incompatible(BooleanAst::Lit(true), BooleanAst::Lit(false), None)]
    fn test_boolean_ast_meet(
        #[case] a: BooleanAst,
        #[case] b: BooleanAst,
        #[case] expected: Option<BooleanAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::same(BooleanAst::Lit(true), BooleanAst::Lit(true), BooleanAst::Lit(true))]
    #[case::differ(
        BooleanAst::Lit(true),
        BooleanAst::Lit(false),
        BooleanAst::Undetermined
    )]
    #[case::top(
        BooleanAst::Undetermined,
        BooleanAst::Lit(true),
        BooleanAst::Undetermined
    )]
    fn test_boolean_ast_join(
        #[case] a: BooleanAst,
        #[case] b: BooleanAst,
        #[case] expected: BooleanAst,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rstest]
    #[case::lit(BooleanAst::Lit(true))]
    #[case::undetermined(BooleanAst::Undetermined)]
    fn test_boolean_ast_canonicalize(#[case] b: BooleanAst) {
        assert_eq!(b.canonicalize(), Ok(b));
    }
}
