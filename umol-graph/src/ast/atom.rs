//! Atom structural AST.

use umol_shared::atom_ast::{AromaticAst, ElementAst, HydrogenAst, IsotopeAst};
use umol_shared::element::Element;
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::config::AtomAstConfig;
use crate::ast::Ast;

/// Atom AST: structural representation of an atom (ground or pattern).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct AtomAst {
    pub element: ElementAst,
    pub isotope_mass: Option<IsotopeAst>,
    pub charge: Option<ValueAst>,
    pub implicit_hydrogens: Option<HydrogenAst>,
    pub lone_pairs: Option<ValueAst>,
    pub spin: Option<SpinStateAst>,
    pub valence: Option<ValueAst>,
    pub donated_pairs: Option<ValueAst>,
    pub accepted_pairs: Option<ValueAst>,
    pub aromatic_valence: Option<AromaticAst>,
    pub multicenter_valence: Option<ValueAst>,
}

impl AtomAst {
    pub fn new(element: ElementAst) -> Self {
        Self {
            element,
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            lone_pairs: None,
            spin: None,
            valence: None,
            donated_pairs: None,
            accepted_pairs: None,
            aromatic_valence: None,
            multicenter_valence: None,
        }
    }

    pub fn from_element(element: Element) -> Self {
        Self::new(ElementAst::Lit(element))
    }
}

impl Ast for AtomAst {
    type Config = AtomAstConfig;
}
