//! Atom structural AST.

use umol_shared::atom_ast::{AromaticValenceAst, ElementAst, HydrogenAst, IsotopeAst};
use umol_shared::element::Element;
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::config::AtomAstConfig;
use crate::ast::Ast;

fn match_element(query: &ElementAst, target: &ElementAst) -> bool {
    match target {
        ElementAst::Lit(e) => query.matches(e),
        _ => false,
    }
}

fn match_option<Q, T>(query: &Option<Q>, target: &Option<T>, f: impl FnOnce(&Q, &T) -> bool) -> bool {
    match (query, target) {
        (None, _) => true,
        (Some(_), None) => false,
        (Some(q), Some(t)) => f(q, t),
    }
}

fn match_option_value(query: &Option<ValueAst>, target: &Option<ValueAst>) -> bool {
    match (query, target) {
        (None, _) => true,
        (Some(pattern), Some(ValueAst::Lit(n))) => pattern.matches(*n),
        (Some(_), _) => false,
    }
}

fn match_option_spin(query: &Option<SpinStateAst>, target: &Option<SpinStateAst>) -> bool {
    match (query, target) {
        (None, _) => true,
        (Some(pattern), Some(SpinStateAst::Lit(s))) => pattern.matches(*s),
        (Some(_), _) => false,
    }
}

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

    pub fn matches_ground(&self, target: &AtomAst) -> bool {
        match_element(&self.element, &target.element)
            && match_option(&self.isotope_mass, &target.isotope_mass, |q, t| q.matches(t))
            && match_option_value(&self.charge, &target.charge)
            && match_option(&self.implicit_hydrogens, &target.implicit_hydrogens, |q, t| q.matches(t))
            && match_option_value(&self.lone_pairs, &target.lone_pairs)
            && match_option_spin(&self.spin, &target.spin)
            && match_option_value(&self.valence, &target.valence)
            && match_option_value(&self.donated_pairs, &target.donated_pairs)
            && match_option_value(&self.accepted_pairs, &target.accepted_pairs)
            && match_option(&self.aromatic_valence, &target.aromatic_valence, |q, t| q.matches(t))
            && match_option_value(&self.multicenter_valence, &target.multicenter_valence)
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
