//! Atom structural AST.

use umol_shared::atom_ast::{AromaticValenceAst, ElementAst, HydrogenAst, IsotopeAst};
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
    pub aromatic_valence: Option<AromaticValenceAst>,
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

    pub fn is_ground(&self) -> bool {
        self.element.is_ground()
            && self.isotope_mass.as_ref().map_or(true, |v| v.is_ground())
            && self.charge.as_ref().map_or(true, |v| v.is_ground())
            && self.implicit_hydrogens.as_ref().map_or(true, |v| v.is_ground())
            && self.lone_pairs.as_ref().map_or(true, |v| v.is_ground())
            && self.spin.as_ref().map_or(true, |v| v.is_ground())
            && self.valence.as_ref().map_or(true, |v| v.is_ground())
            && self.donated_pairs.as_ref().map_or(true, |v| v.is_ground())
            && self.accepted_pairs.as_ref().map_or(true, |v| v.is_ground())
            && self.aromatic_valence.as_ref().map_or(true, |v| v.is_ground())
            && self.multicenter_valence.as_ref().map_or(true, |v| v.is_ground())
    }
}

impl Ast for AtomAst {
    type Config = AtomAstConfig;
}
