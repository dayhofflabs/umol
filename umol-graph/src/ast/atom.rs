//! Atom structural AST.

use umol_shared::atom_ast::{AromaticValenceAst, ElementAst, HydrogenAst, IsotopeAst};
use umol_shared::element::Element;
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::config::AtomAstConfig;
use crate::ast::Ast;

/// Atom AST: structural representation of an atom (ground or pattern).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AtomAst {
    pub element: ElementAst,
    pub isotope_mass: IsotopeAst,
    pub charge: ValueAst,
    pub implicit_hydrogens: HydrogenAst,
    pub lone_pairs: ValueAst,
    pub spin: SpinStateAst,
    pub valence: ValueAst,
    pub donated_pairs: ValueAst,
    pub accepted_pairs: ValueAst,
    pub aromatic_valence: AromaticValenceAst,
    pub multicenter_valence: ValueAst,
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

    pub fn matches_ground(&self, target: &AtomAst) -> bool {
        (match &target.element {
            ElementAst::Lit(e) => self.element.matches(e),
            _ => false,
        }) && self.isotope_mass.matches(&target.isotope_mass)
            && (match &target.charge {
                ValueAst::Lit(n) => self.charge.matches(*n),
                _ => matches!(self.charge, ValueAst::Undetermined),
            })
            && self.implicit_hydrogens.matches(&target.implicit_hydrogens)
            && (match &target.lone_pairs {
                ValueAst::Lit(n) => self.lone_pairs.matches(*n),
                _ => matches!(self.lone_pairs, ValueAst::Undetermined),
            })
            && (match &target.spin {
                SpinStateAst::Lit(s) => self.spin.matches(*s),
                _ => matches!(self.spin, SpinStateAst::Pair { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Undetermined }),
            })
            && (match &target.valence {
                ValueAst::Lit(n) => self.valence.matches(*n),
                _ => matches!(self.valence, ValueAst::Undetermined),
            })
            && (match &target.donated_pairs {
                ValueAst::Lit(n) => self.donated_pairs.matches(*n),
                _ => matches!(self.donated_pairs, ValueAst::Undetermined),
            })
            && (match &target.accepted_pairs {
                ValueAst::Lit(n) => self.accepted_pairs.matches(*n),
                _ => matches!(self.accepted_pairs, ValueAst::Undetermined),
            })
            && self.aromatic_valence.matches(&target.aromatic_valence)
            && (match &target.multicenter_valence {
                ValueAst::Lit(n) => self.multicenter_valence.matches(*n),
                _ => matches!(self.multicenter_valence, ValueAst::Undetermined),
            })
    }

    pub fn is_ground(&self) -> bool {
        self.element.is_ground()
            && self.isotope_mass.is_ground()
            && self.charge.is_ground()
            && self.implicit_hydrogens.is_ground()
            && self.lone_pairs.is_ground()
            && self.spin.is_ground()
            && self.valence.is_ground()
            && self.donated_pairs.is_ground()
            && self.accepted_pairs.is_ground()
            && self.aromatic_valence.is_ground()
            && self.multicenter_valence.is_ground()
    }
}

impl Ast for AtomAst {
    type Config = AtomAstConfig;
}
