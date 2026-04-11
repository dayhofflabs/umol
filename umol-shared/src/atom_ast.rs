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
}

/// Aromatic valence expressions
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticAst {
    Unspecified,
    NotAromatic,
    Value(ValueAst),
}
