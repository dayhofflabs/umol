//! Atom-level AST fragments shared across crates.

use crate::element::Element;
use crate::value_ast::ValueAst;

/// Element expressions
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ElementAst {
    Lit(Element),
    Undetermined,
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

    pub fn matches(&self, target: &Element) -> bool {
        match self {
            Self::Lit(e) => e == target,
            Self::Undetermined => true,
            Self::Set(s) => s.contains(target),
            Self::Bind { set, .. } => set.contains(target),
            Self::Ref(_) => false,
        }
    }
}

/// Isotope-mass expressions (Natural = #i=)
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum IsotopeAst {
    Natural,
    Lit(u32),
    Undetermined,
    Set(Vec<u32>),
    Bind { id: String, set: Vec<u32> },
    Ref(String),
}

impl IsotopeAst {
    pub fn is_ground(&self) -> bool {
        matches!(self, Self::Natural | Self::Lit(_))
    }

    pub fn matches(&self, target: &IsotopeAst) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (Self::Natural, Self::Natural) => true,
            (Self::Lit(a), Self::Lit(b)) => a == b,
            (Self::Set(s), Self::Lit(b)) => s.contains(b),
            (Self::Bind { set, .. }, Self::Lit(b)) => set.contains(b),
            _ => false,
        }
    }
}

/// Implicit hydrogen expressions (Normal = #h=)
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum HydrogenAst {
    Undetermined,
    Normal,
    Value(ValueAst),
}

impl HydrogenAst {
    pub fn from_value(value: ValueAst) -> Self {
        Self::Value(value)
    }

    pub fn is_ground(&self) -> bool {
        match self {
            Self::Undetermined => false,
            Self::Normal => true,
            Self::Value(v) => v.is_ground(),
        }
    }

    pub fn matches(&self, target: &HydrogenAst) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (Self::Normal, Self::Normal) => true,
            (Self::Value(pattern), Self::Value(ValueAst::Lit(n))) => pattern.matches(*n),
            _ => false,
        }
    }
}

/// Aromatic valence expressions
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticValenceAst {
    Undetermined,
    NotAromatic,
    Value(ValueAst),
}

impl AromaticValenceAst {
    pub fn is_ground(&self) -> bool {
        match self {
            Self::Undetermined => false,
            Self::NotAromatic => true,
            Self::Value(v) => v.is_ground(),
        }
    }

    pub fn matches(&self, target: &AromaticValenceAst) -> bool {
        match (self, target) {
            (Self::Undetermined, Self::Undetermined) => true,
            (Self::NotAromatic, Self::NotAromatic) => true,
            (Self::Value(pattern), Self::Value(ValueAst::Lit(n))) => pattern.matches(*n),
            _ => false,
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
    #[case::wildcard(ElementAst::Undetermined, false)]
    #[case::set(ElementAst::Set(vec![Element::C, Element::N]), false)]
    #[case::bind(ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, false)]
    #[case::reference(ElementAst::Ref("e".into()), false)]
    fn test_element_ast_is_ground(#[case] ast: ElementAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::natural(IsotopeAst::Natural, true)]
    #[case::lit(IsotopeAst::Lit(12), true)]
    #[case::wildcard(IsotopeAst::Undetermined, false)]
    #[case::set(IsotopeAst::Set(vec![12, 13]), false)]
    #[case::bind(IsotopeAst::Bind { id: "i".into(), set: vec![12] }, false)]
    #[case::reference(IsotopeAst::Ref("i".into()), false)]
    fn test_isotope_ast_is_ground(#[case] ast: IsotopeAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::normal(HydrogenAst::Normal, true)]
    #[case::lit(HydrogenAst::Value(ValueAst::Lit(2)), true)]
    #[case::wildcard(HydrogenAst::Value(ValueAst::Undetermined), false)]
    fn test_hydrogen_ast_is_ground(#[case] ast: HydrogenAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceAst::Undetermined, false)]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, true)]
    #[case::lit(AromaticValenceAst::Value(ValueAst::Lit(2)), true)]
    #[case::wildcard(AromaticValenceAst::Value(ValueAst::Undetermined), false)]
    fn test_aromatic_ast_is_ground(#[case] ast: AromaticValenceAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::lit_match(ElementAst::Lit(Element::C), Element::C, true)]
    #[case::lit_mismatch(ElementAst::Lit(Element::C), Element::N, false)]
    #[case::wildcard(ElementAst::Undetermined, Element::C, true)]
    #[case::set_match(ElementAst::Set(vec![Element::C, Element::N]), Element::N, true)]
    #[case::set_mismatch(ElementAst::Set(vec![Element::C, Element::N]), Element::O, false)]
    #[case::bind_match(ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, Element::C, true)]
    #[case::bind_mismatch(ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, Element::N, false)]
    #[case::ref_no_match(ElementAst::Ref("e".into()), Element::C, false)]
    fn test_element_ast_matches_element(#[case] ast: ElementAst, #[case] target: Element, #[case] expected: bool) {
        assert_eq!(ast.matches(&target), expected);
    }

    #[rstest]
    #[case::natural_match(IsotopeAst::Natural, IsotopeAst::Natural, true)]
    #[case::lit_match(IsotopeAst::Lit(12), IsotopeAst::Lit(12), true)]
    #[case::lit_mismatch(IsotopeAst::Lit(12), IsotopeAst::Lit(13), false)]
    #[case::wildcard(IsotopeAst::Undetermined, IsotopeAst::Lit(12), true)]
    #[case::wildcard_natural(IsotopeAst::Undetermined, IsotopeAst::Natural, true)]
    #[case::set_match(IsotopeAst::Set(vec![12, 13]), IsotopeAst::Lit(13), true)]
    #[case::set_mismatch(IsotopeAst::Set(vec![12, 13]), IsotopeAst::Lit(14), false)]
    #[case::set_vs_natural(IsotopeAst::Set(vec![12]), IsotopeAst::Natural, false)]
    #[case::natural_vs_lit(IsotopeAst::Natural, IsotopeAst::Lit(12), false)]
    fn test_isotope_ast_matches_isotope(#[case] pattern: IsotopeAst, #[case] target: IsotopeAst, #[case] expected: bool) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::normal_match(HydrogenAst::Normal, HydrogenAst::Normal, true)]
    #[case::lit_match(HydrogenAst::Value(ValueAst::Lit(2)), HydrogenAst::Value(ValueAst::Lit(2)), true)]
    #[case::lit_mismatch(HydrogenAst::Value(ValueAst::Lit(2)), HydrogenAst::Value(ValueAst::Lit(3)), false)]
    #[case::wildcard_match(HydrogenAst::Value(ValueAst::Undetermined), HydrogenAst::Value(ValueAst::Lit(2)), true)]
    #[case::normal_vs_value(HydrogenAst::Normal, HydrogenAst::Value(ValueAst::Lit(0)), false)]
    #[case::value_vs_normal(HydrogenAst::Value(ValueAst::Lit(0)), HydrogenAst::Normal, false)]
    fn test_hydrogen_ast_matches_hydrogen(#[case] pattern: HydrogenAst, #[case] target: HydrogenAst, #[case] expected: bool) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::unspecified_match(AromaticValenceAst::Undetermined, AromaticValenceAst::Undetermined, true)]
    #[case::not_aromatic_match(AromaticValenceAst::NotAromatic, AromaticValenceAst::NotAromatic, true)]
    #[case::lit_match(AromaticValenceAst::Value(ValueAst::Lit(2)), AromaticValenceAst::Value(ValueAst::Lit(2)), true)]
    #[case::lit_mismatch(AromaticValenceAst::Value(ValueAst::Lit(2)), AromaticValenceAst::Value(ValueAst::Lit(3)), false)]
    #[case::wildcard_match(AromaticValenceAst::Value(ValueAst::Undetermined), AromaticValenceAst::Value(ValueAst::Lit(2)), true)]
    #[case::unspecified_vs_not_aromatic(AromaticValenceAst::Undetermined, AromaticValenceAst::NotAromatic, false)]
    #[case::value_vs_unspecified(AromaticValenceAst::Value(ValueAst::Lit(2)), AromaticValenceAst::Undetermined, false)]
    fn test_aromatic_valence_ast_matches(#[case] pattern: AromaticValenceAst, #[case] target: AromaticValenceAst, #[case] expected: bool) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
