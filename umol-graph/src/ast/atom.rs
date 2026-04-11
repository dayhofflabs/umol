//! Atom structural AST.

use umol_data::Element;

use crate::ast::config::AtomAstConfig;
use crate::ast::value::ValueAst;
use crate::ast::Ast;

/// Atom AST: structural representation of an atom (ground or pattern).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct AtomAst {
    pub element: ElementExpr,
    pub isotope_mass: Option<IsotopeExpr>,
    pub charge: Option<ValueAst>,
    pub implicit_hydrogens: Option<HydrogenExpr>,
    pub lone_pairs: Option<ValueAst>,
    pub unpaired_electrons: Option<ValueAst>,
    pub multiplicity: Option<ValueAst>,
    pub valence: Option<ValueAst>,
    pub donated_pairs: Option<ValueAst>,
    pub accepted_pairs: Option<ValueAst>,
    pub aromatic_valence: Option<AromaticExpr>,
    pub multicenter_valence: Option<ValueAst>,
}

impl AtomAst {
    pub fn new(element: ElementExpr) -> Self {
        Self {
            element,
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            valence: None,
            donated_pairs: None,
            accepted_pairs: None,
            aromatic_valence: None,
            multicenter_valence: None,
        }
    }

    pub fn from_element(element: Element) -> Self {
        Self::new(ElementExpr::Lit(element))
    }
}

impl Ast for AtomAst {
    type Config = AtomAstConfig;
}

/// Element expressions
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ElementExpr {
    Lit(Element),
    Wildcard,
    Set(Vec<Element>),
    Bind { id: String, set: Vec<Element> },
    Ref(String),
}

impl ElementExpr {
    pub fn new(element: Element) -> Self {
        Self::Lit(element)
    }
}

/// Isotope-mass expressions (Natural = #i=)
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum IsotopeExpr {
    Natural,
    Lit(u32),
    Wildcard,
    Set(Vec<u32>),
    Bind { id: String, set: Vec<u32> },
    Ref(String),
}

/// Implicit hydrogen expressions (Normal = #h=)
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum HydrogenExpr {
    Normal,
    Value(ValueAst),
}

impl HydrogenExpr {
    pub fn from_value(value: ValueAst) -> Self {
        Self::Value(value)
    }
}

/// Aromatic valence expressions
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticExpr {
    Unspecified,
    NotAromatic,
    Value(ValueAst),
}
