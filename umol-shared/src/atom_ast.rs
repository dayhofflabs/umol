//! Atom-level AST fragments shared across crates.

use crate::element::Element;
use crate::value_ast::ValueAst;

/// Element expressions
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ElementAst {
    Lit(Element),
    Wildcard,
    Set(Vec<Element>),
    Bind { id: String, set: Vec<Element> },
    Ref(String),
}

impl ElementAst {
    pub fn new(element: Element) -> Self {
        Self::Lit(element)
    }

    pub fn is_ground(&self) -> bool {
        matches!(self, Self::Lit(_))
    }
}

/// Isotope-mass expressions (Natural = #i=)
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum IsotopeAst {
    Natural,
    Lit(u32),
    Wildcard,
    Set(Vec<u32>),
    Bind { id: String, set: Vec<u32> },
    Ref(String),
}

impl IsotopeAst {
    pub fn is_ground(&self) -> bool {
        matches!(self, Self::Natural | Self::Lit(_))
    }
}

/// Implicit hydrogen expressions (Normal = #h=)
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum HydrogenAst {
    Normal,
    Value(ValueAst),
}

impl HydrogenAst {
    pub fn from_value(value: ValueAst) -> Self {
        Self::Value(value)
    }

    pub fn is_ground(&self) -> bool {
        match self {
            Self::Normal => true,
            Self::Value(v) => v.is_ground(),
        }
    }
}

/// Aromatic valence expressions
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticValenceAst {
    Unspecified,
    NotAromatic,
    Value(ValueAst),
}

impl AromaticValenceAst {
    pub fn is_ground(&self) -> bool {
        match self {
            Self::Unspecified | Self::NotAromatic => true,
            Self::Value(v) => v.is_ground(),
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::lit(ElementAst::Lit(Element::C), true)]
    #[case::wildcard(ElementAst::Wildcard, false)]
    #[case::set(ElementAst::Set(vec![Element::C, Element::N]), false)]
    #[case::bind(ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, false)]
    #[case::reference(ElementAst::Ref("e".into()), false)]
    fn test_element_ast_is_ground(#[case] ast: ElementAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::natural(IsotopeAst::Natural, true)]
    #[case::lit(IsotopeAst::Lit(12), true)]
    #[case::wildcard(IsotopeAst::Wildcard, false)]
    #[case::set(IsotopeAst::Set(vec![12, 13]), false)]
    #[case::bind(IsotopeAst::Bind { id: "i".into(), set: vec![12] }, false)]
    #[case::reference(IsotopeAst::Ref("i".into()), false)]
    fn test_isotope_ast_is_ground(#[case] ast: IsotopeAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::normal(HydrogenAst::Normal, true)]
    #[case::lit(HydrogenAst::Value(ValueAst::Lit(2)), true)]
    #[case::wildcard(HydrogenAst::Value(ValueAst::Wildcard), false)]
    fn test_hydrogen_ast_is_ground(#[case] ast: HydrogenAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::unspecified(AromaticValenceAst::Unspecified, true)]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, true)]
    #[case::lit(AromaticValenceAst::Value(ValueAst::Lit(2)), true)]
    #[case::wildcard(AromaticValenceAst::Value(ValueAst::Wildcard), false)]
    fn test_aromatic_ast_is_ground(#[case] ast: AromaticValenceAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }
}
