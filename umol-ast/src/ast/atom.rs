//! Atom-level AST fragments shared across crates.

use umol_shared::element::Element;

use super::constraint::AtomConstraint;
use super::spin::SpinStateAst;
use super::value::ValueAst;

/// Atom AST: structural representation of an atom plus the atom-level
/// constraints (valence, degree, ring membership, etc.) that pattern
/// against the surrounding topology.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AtomAst {
    pub element: ElementAst,
    pub isotope_mass: IsotopeAst,
    pub charge: ValueAst,
    pub implicit_hydrogens: ImplicitHydrogensAst,
    pub lone_pairs: ValueAst,
    pub spin: SpinStateAst,
    pub constraints: Vec<AtomConstraint>,
}

impl AtomAst {
    pub fn new(element: ElementAst) -> Self {
        Self {
            element,
            ..Default::default()
        }
    }

    pub fn from_element(element: Element) -> Self {
        Self::new(ElementAst::Lit(element))
    }

    pub fn is_ground(&self) -> bool {
        self.element.is_ground()
            && self.isotope_mass.is_ground()
            && self.charge.is_ground()
            && self.implicit_hydrogens.is_ground()
            && self.lone_pairs.is_ground()
            && self.spin.is_ground()
    }

    /// `self` (pattern) matches `target` iff every admissible assignment
    /// of `target` is also admissible by `self`, checked field-wise.
    /// See per-field `matches` for the scalar rules.
    pub fn matches(&self, target: &AtomAst) -> bool {
        self.element.matches(&target.element)
            && self.isotope_mass.matches(&target.isotope_mass)
            && self.charge.matches(&target.charge)
            && self.implicit_hydrogens.matches(&target.implicit_hydrogens)
            && self.lone_pairs.matches(&target.lone_pairs)
            && self.spin.matches(&target.spin)
    }
}

/// Element expressions
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum ElementAst {
    Lit(Element),
    #[default]
    Undetermined,
    Set(Vec<Element>),
    Bind {
        id: String,
        set: Vec<Element>,
    },
    Ref(String),
}

impl ElementAst {
    pub fn new(element: Element) -> Self {
        Self::Lit(element)
    }

    pub fn is_ground(&self) -> bool {
        matches!(self, Self::Lit(_))
    }

    /// Pattern matches target iff every element the target admits is also
    /// admitted by the pattern (superset semantics).
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

/// Isotope-mass expressions (Natural = #i=)
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum IsotopeAst {
    #[default]
    Undetermined,
    Natural,
    Value(ValueAst),
}

impl IsotopeAst {
    pub fn new(mass: u32) -> Self {
        Self::Value(ValueAst::Lit(mass as i64))
    }

    pub fn from_value(value: ValueAst) -> Self {
        Self::Value(value)
    }

    pub fn is_ground(&self) -> bool {
        match self {
            Self::Undetermined => false,
            Self::Natural => true,
            Self::Value(v) => v.is_ground(),
        }
    }

    pub fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::Natural, Self::Natural) => true,
            (Self::Value(p), Self::Value(t)) => p.matches(t),
            _ => false,
        }
    }
}

/// Implicit hydrogen expressions (Normal = #h=)
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum ImplicitHydrogensAst {
    #[default]
    Undetermined,
    Normal,
    Value(ValueAst),
}

impl ImplicitHydrogensAst {
    pub fn new(count: u8) -> Self {
        Self::Value(ValueAst::Lit(count as i64))
    }

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

    pub fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::Normal, Self::Normal) => true,
            (Self::Value(p), Self::Value(t)) => p.matches(t),
            _ => false,
        }
    }
}

/// Aromatic valence expressions
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum AromaticValenceAst {
    #[default]
    Undetermined,
    NotAromatic,
    Value(ValueAst),
}

impl AromaticValenceAst {
    pub fn new(valence: u8) -> Self {
        Self::Value(ValueAst::Lit(valence as i64))
    }

    pub fn from_value(value: ValueAst) -> Self {
        Self::Value(value)
    }

    pub fn is_ground(&self) -> bool {
        match self {
            Self::Undetermined => false,
            Self::NotAromatic => true,
            Self::Value(v) => v.is_ground(),
        }
    }

    pub fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::NotAromatic, Self::NotAromatic) => true,
            (Self::Value(p), Self::Value(t)) => p.matches(t),
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
    #[case::new_element(ElementAst::new(Element::C), ElementAst::Lit(Element::C))]
    fn test_element_ast_new(#[case] actual: ElementAst, #[case] expected: ElementAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::new(IsotopeAst::new(12), IsotopeAst::Value(ValueAst::Lit(12)))]
    #[case::from_value_lit(IsotopeAst::from_value(ValueAst::Lit(13)), IsotopeAst::Value(ValueAst::Lit(13)))]
    #[case::from_value_undetermined(IsotopeAst::from_value(ValueAst::Undetermined), IsotopeAst::Value(ValueAst::Undetermined))]
    fn test_isotope_ast_new(#[case] actual: IsotopeAst, #[case] expected: IsotopeAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::new(ImplicitHydrogensAst::new(3), ImplicitHydrogensAst::Value(ValueAst::Lit(3)))]
    #[case::from_value(ImplicitHydrogensAst::from_value(ValueAst::LitSet(vec![0, 1])), ImplicitHydrogensAst::Value(ValueAst::LitSet(vec![0, 1])))]
    fn test_implicit_hydrogens_ast_new(
        #[case] actual: ImplicitHydrogensAst,
        #[case] expected: ImplicitHydrogensAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::new(AromaticValenceAst::new(2), AromaticValenceAst::Value(ValueAst::Lit(2)))]
    #[case::from_value(AromaticValenceAst::from_value(ValueAst::Undetermined), AromaticValenceAst::Value(ValueAst::Undetermined))]
    fn test_aromatic_valence_ast_new(
        #[case] actual: AromaticValenceAst,
        #[case] expected: AromaticValenceAst,
    ) {
        assert_eq!(actual, expected);
    }

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
    #[case::lit(IsotopeAst::Value(ValueAst::Lit(12)), true)]
    #[case::wildcard(IsotopeAst::Undetermined, false)]
    #[case::value_wildcard(IsotopeAst::Value(ValueAst::Undetermined), false)]
    #[case::set(IsotopeAst::Value(ValueAst::LitSet(vec![12, 13])), false)]
    fn test_isotope_ast_is_ground(#[case] ast: IsotopeAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::normal(ImplicitHydrogensAst::Normal, true)]
    #[case::lit(ImplicitHydrogensAst::Value(ValueAst::Lit(2)), true)]
    #[case::wildcard(ImplicitHydrogensAst::Value(ValueAst::Undetermined), false)]
    fn test_implicit_hydrogens_ast_is_ground(
        #[case] ast: ImplicitHydrogensAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::undetermined(AromaticValenceAst::Undetermined, false)]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, true)]
    #[case::lit(AromaticValenceAst::Value(ValueAst::Lit(2)), true)]
    #[case::wildcard(AromaticValenceAst::Value(ValueAst::Undetermined), false)]
    fn test_aromatic_valence_ast_is_ground(
        #[case] ast: AromaticValenceAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }


    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_lit(ElementAst::Undetermined, ElementAst::Lit(Element::C), true)]
    #[case::undetermined_undetermined(ElementAst::Undetermined, ElementAst::Undetermined, true)]
    #[case::undetermined_set(ElementAst::Undetermined, ElementAst::Set(vec![Element::C, Element::N]), true)]
    #[case::lit_undetermined(ElementAst::Lit(Element::C), ElementAst::Undetermined, false)]
    #[case::set_undetermined(ElementAst::Set(vec![Element::C]), ElementAst::Undetermined, false)]
    #[case::lit_lit_match(ElementAst::Lit(Element::C), ElementAst::Lit(Element::C), true)]
    #[case::lit_lit_mismatch(ElementAst::Lit(Element::C), ElementAst::Lit(Element::N), false)]
    #[case::lit_singleton_set(ElementAst::Lit(Element::C), ElementAst::Set(vec![Element::C]), true)]
    #[case::lit_multi_set(ElementAst::Lit(Element::C), ElementAst::Set(vec![Element::C, Element::N]), false)]
    #[case::set_lit_in(ElementAst::Set(vec![Element::C, Element::N]), ElementAst::Lit(Element::N), true)]
    #[case::set_lit_out(ElementAst::Set(vec![Element::C, Element::N]), ElementAst::Lit(Element::O), false)]
    #[case::set_set_subset(ElementAst::Set(vec![Element::C, Element::N, Element::O]), ElementAst::Set(vec![Element::C, Element::N]), true)]
    #[case::set_set_equal(ElementAst::Set(vec![Element::C, Element::N]), ElementAst::Set(vec![Element::C, Element::N]), true)]
    #[case::set_set_superset(ElementAst::Set(vec![Element::C]), ElementAst::Set(vec![Element::C, Element::N]), false)]
    #[case::bind_lit_match(ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, ElementAst::Lit(Element::C), true)]
    #[case::bind_lit_mismatch(ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, ElementAst::Lit(Element::N), false)]
    #[case::ref_lit(ElementAst::Ref("e".into()), ElementAst::Lit(Element::C), false)]
    #[case::lit_ref(ElementAst::Lit(Element::C), ElementAst::Ref("e".into()), false)]
    fn test_element_ast_matches(
        #[case] pattern: ElementAst,
        #[case] target: ElementAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_natural(IsotopeAst::Undetermined, IsotopeAst::Natural, true)]
    #[case::undetermined_value(IsotopeAst::Undetermined, IsotopeAst::Value(ValueAst::Lit(12)), true)]
    #[case::undetermined_undetermined(IsotopeAst::Undetermined, IsotopeAst::Undetermined, true)]
    #[case::natural_undetermined(IsotopeAst::Natural, IsotopeAst::Undetermined, false)]
    #[case::value_undetermined(IsotopeAst::Value(ValueAst::Lit(12)), IsotopeAst::Undetermined, false)]
    #[case::natural_natural(IsotopeAst::Natural, IsotopeAst::Natural, true)]
    #[case::natural_value(IsotopeAst::Natural, IsotopeAst::Value(ValueAst::Lit(12)), false)]
    #[case::value_natural(IsotopeAst::Value(ValueAst::Lit(12)), IsotopeAst::Natural, false)]
    #[case::value_lit_match(IsotopeAst::Value(ValueAst::Lit(12)), IsotopeAst::Value(ValueAst::Lit(12)), true)]
    #[case::value_lit_mismatch(IsotopeAst::Value(ValueAst::Lit(12)), IsotopeAst::Value(ValueAst::Lit(13)), false)]
    #[case::value_wildcard_lit(IsotopeAst::Value(ValueAst::Undetermined), IsotopeAst::Value(ValueAst::Lit(12)), true)]
    #[case::value_set_lit_in(IsotopeAst::Value(ValueAst::LitSet(vec![12, 13])), IsotopeAst::Value(ValueAst::Lit(13)), true)]
    #[case::value_set_lit_out(IsotopeAst::Value(ValueAst::LitSet(vec![12, 13])), IsotopeAst::Value(ValueAst::Lit(14)), false)]
    #[case::value_set_set_subset(IsotopeAst::Value(ValueAst::LitSet(vec![12, 13, 14])), IsotopeAst::Value(ValueAst::LitSet(vec![12, 13])), true)]
    #[case::value_set_set_superset(IsotopeAst::Value(ValueAst::LitSet(vec![12])), IsotopeAst::Value(ValueAst::LitSet(vec![12, 13])), false)]
    #[case::value_lit_wildcard(IsotopeAst::Value(ValueAst::Lit(12)), IsotopeAst::Value(ValueAst::Undetermined), false)]
    fn test_isotope_ast_matches(
        #[case] pattern: IsotopeAst,
        #[case] target: IsotopeAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_normal(ImplicitHydrogensAst::Undetermined, ImplicitHydrogensAst::Normal, true)]
    #[case::undetermined_value(ImplicitHydrogensAst::Undetermined, ImplicitHydrogensAst::Value(ValueAst::Lit(3)), true)]
    #[case::normal_undetermined(ImplicitHydrogensAst::Normal, ImplicitHydrogensAst::Undetermined, false)]
    #[case::normal_normal(ImplicitHydrogensAst::Normal, ImplicitHydrogensAst::Normal, true)]
    #[case::normal_value(ImplicitHydrogensAst::Normal, ImplicitHydrogensAst::Value(ValueAst::Lit(0)), false)]
    #[case::value_normal(ImplicitHydrogensAst::Value(ValueAst::Lit(0)), ImplicitHydrogensAst::Normal, false)]
    #[case::value_lit_match(ImplicitHydrogensAst::Value(ValueAst::Lit(2)), ImplicitHydrogensAst::Value(ValueAst::Lit(2)), true)]
    #[case::value_lit_mismatch(ImplicitHydrogensAst::Value(ValueAst::Lit(2)), ImplicitHydrogensAst::Value(ValueAst::Lit(3)), false)]
    #[case::value_wildcard(ImplicitHydrogensAst::Value(ValueAst::Undetermined), ImplicitHydrogensAst::Value(ValueAst::Lit(2)), true)]
    #[case::value_set_subset(ImplicitHydrogensAst::Value(ValueAst::LitSet(vec![1, 2])), ImplicitHydrogensAst::Value(ValueAst::LitSet(vec![1])), true)]
    fn test_implicit_hydrogens_ast_matches(
        #[case] pattern: ImplicitHydrogensAst,
        #[case] target: ImplicitHydrogensAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_undetermined(AromaticValenceAst::Undetermined, AromaticValenceAst::Undetermined, true)]
    #[case::undetermined_not_aromatic(AromaticValenceAst::Undetermined, AromaticValenceAst::NotAromatic, true)]
    #[case::undetermined_value(AromaticValenceAst::Undetermined, AromaticValenceAst::Value(ValueAst::Lit(2)), true)]
    #[case::not_aromatic_undetermined(AromaticValenceAst::NotAromatic, AromaticValenceAst::Undetermined, false)]
    #[case::not_aromatic_not_aromatic(AromaticValenceAst::NotAromatic, AromaticValenceAst::NotAromatic, true)]
    #[case::not_aromatic_value(AromaticValenceAst::NotAromatic, AromaticValenceAst::Value(ValueAst::Lit(2)), false)]
    #[case::value_not_aromatic(AromaticValenceAst::Value(ValueAst::Lit(2)), AromaticValenceAst::NotAromatic, false)]
    #[case::value_lit_match(AromaticValenceAst::Value(ValueAst::Lit(2)), AromaticValenceAst::Value(ValueAst::Lit(2)), true)]
    #[case::value_lit_mismatch(AromaticValenceAst::Value(ValueAst::Lit(2)), AromaticValenceAst::Value(ValueAst::Lit(3)), false)]
    #[case::value_wildcard(AromaticValenceAst::Value(ValueAst::Undetermined), AromaticValenceAst::Value(ValueAst::Lit(2)), true)]
    fn test_aromatic_valence_ast_matches(
        #[case] pattern: AromaticValenceAst,
        #[case] target: AromaticValenceAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard_vs_ground(AtomAst::default(), AtomAst::from_element(Element::C), true)]
    #[case::same_element(AtomAst::from_element(Element::C), AtomAst::from_element(Element::C), true)]
    #[case::element_mismatch(AtomAst::from_element(Element::C), AtomAst::from_element(Element::N), false)]
    #[case::pattern_more_specific_than_target(AtomAst::from_element(Element::C), AtomAst::default(), false)]
    #[case::charge_mismatch(AtomAst { element: ElementAst::Lit(Element::C), charge: ValueAst::Lit(1),..Default::default() },
                            AtomAst { element: ElementAst::Lit(Element::C), charge: ValueAst::Lit(0),..Default::default() }, false)]
    #[case::charge_wildcard_pattern(AtomAst::from_element(Element::C), AtomAst { element: ElementAst::Lit(Element::C), charge: ValueAst::Lit(1),..Default::default() }, true)]
    fn test_atom_ast_matches(
        #[case] pattern: AtomAst,
        #[case] target: AtomAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
