//! Boolean form.

use std::borrow::Cow;

use super::error::{Contradiction, NoJoin};
use super::traits::{AsLit, Lattice, Normalize};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BooleanForm {
    #[default]
    Undetermined,
    Lit(bool),
}

impl BooleanForm {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn lit(b: bool) -> Self {
        Self::Lit(b)
    }
}

impl From<bool> for BooleanForm {
    fn from(b: bool) -> Self {
        Self::lit(b)
    }
}

impl Normalize for BooleanForm {
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(self)
    }

    fn normalized(&self) -> Result<Cow<'_, Self>, Contradiction> {
        Ok(Cow::Borrowed(self))
    }
}

impl AsLit for BooleanForm {
    type Lit = bool;

    fn as_lit(&self) -> Option<bool> {
        match self {
            Self::Lit(b) => Some(*b),
            Self::Undetermined => None,
        }
    }
}

impl Lattice for BooleanForm {
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

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        Ok(match (self, other) {
            (Self::Lit(a), Self::Lit(b)) if a == b => Self::Lit(*a),
            _ => Self::Undetermined,
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::lit_true(BooleanForm::Lit(true), Some(true))]
    #[case::lit_false(BooleanForm::Lit(false), Some(false))]
    #[case::undetermined(BooleanForm::Undetermined, None)]
    fn test_boolean_form_as_lit(#[case] b: BooleanForm, #[case] expected: Option<bool>) {
        assert_eq!(b.as_lit(), expected);
    }

    #[rstest]
    #[case::undetermined(BooleanForm::Undetermined, true)]
    #[case::lit(BooleanForm::Lit(true), false)]
    fn test_boolean_form_is_undetermined(#[case] b: BooleanForm, #[case] expected: bool) {
        assert_eq!(b.is_undetermined(), expected);
    }

    #[rstest]
    #[case::lit(BooleanForm::Lit(false), true)]
    #[case::undetermined(BooleanForm::Undetermined, false)]
    fn test_boolean_form_is_ground(#[case] b: BooleanForm, #[case] expected: bool) {
        assert_eq!(b.is_ground(), expected);
    }

    #[rstest]
    #[case::top_left(
        BooleanForm::Undetermined,
        BooleanForm::Lit(true),
        Some(BooleanForm::Lit(true))
    )]
    #[case::top_right(
        BooleanForm::Lit(false),
        BooleanForm::Undetermined,
        Some(BooleanForm::Lit(false))
    )]
    #[case::same(
        BooleanForm::Lit(true),
        BooleanForm::Lit(true),
        Some(BooleanForm::Lit(true))
    )]
    #[case::incompatible(BooleanForm::Lit(true), BooleanForm::Lit(false), None)]
    fn test_boolean_form_meet(
        #[case] a: BooleanForm,
        #[case] b: BooleanForm,
        #[case] expected: Option<BooleanForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::same(BooleanForm::Lit(true), BooleanForm::Lit(true), BooleanForm::Lit(true))]
    #[case::differ(
        BooleanForm::Lit(true),
        BooleanForm::Lit(false),
        BooleanForm::Undetermined
    )]
    #[case::top(
        BooleanForm::Undetermined,
        BooleanForm::Lit(true),
        BooleanForm::Undetermined
    )]
    fn test_boolean_form_join(
        #[case] a: BooleanForm,
        #[case] b: BooleanForm,
        #[case] expected: BooleanForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rstest]
    #[case::lit(BooleanForm::Lit(true))]
    #[case::undetermined(BooleanForm::Undetermined)]
    fn test_boolean_form_normalize(#[case] b: BooleanForm) {
        assert_eq!(b.normalize(), Ok(b));
    }
}
