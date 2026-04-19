//! Substructure-query pattern backed by [`MoleculeAst`].
//!
//! [`MoleculeAst`]: crate::ast::molecule::MoleculeAst

use std::sync::Arc;

use umol_shared::value_ast::ValueAst;

use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::config::{AromaticValenceMode, AtomAstConfig, NumericMode};
use crate::ast::constraint::{AromaticValenceConstraint, AtomConstraint, BondConstraint};
use crate::ast::molecule::MoleculeAst;

#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AtomPattern {
    pub ast: AtomAst,
    pub constraints: Vec<AtomConstraint>,
}

impl AtomPattern {
    pub fn new(ast: AtomAst) -> Self {
        Self { ast, constraints: Vec::new() }
    }

    pub fn with_constraints(ast: AtomAst, constraints: Vec<AtomConstraint>) -> Self {
        Self { ast, constraints }
    }

    pub fn coerce(&mut self, cfg: &AtomAstConfig) {
        self.ast.coerce(cfg);
        coerce_atom_constraints(&mut self.constraints, cfg);
    }

    pub fn release(&mut self, cfg: &AtomAstConfig) {
        self.ast.release(cfg);
        release_atom_constraints(&mut self.constraints, cfg);
    }
}

type NumericDefault<'a> = (fn(&AtomConstraint) -> bool, AtomConstraint, &'a NumericMode);

pub(crate) fn coerce_atom_constraints(constraints: &mut Vec<AtomConstraint>, cfg: &AtomAstConfig) {
    let numeric_defaults: &[NumericDefault] = &[
        (
            |c| matches!(c, AtomConstraint::Valence(_)),
            AtomConstraint::Valence(ValueAst::Lit(0)),
            &cfg.valence_mode,
        ),
        (
            |c| matches!(c, AtomConstraint::DonatedPairs(_)),
            AtomConstraint::DonatedPairs(ValueAst::Lit(0)),
            &cfg.donated_pairs_mode,
        ),
        (
            |c| matches!(c, AtomConstraint::AcceptedPairs(_)),
            AtomConstraint::AcceptedPairs(ValueAst::Lit(0)),
            &cfg.accepted_pairs_mode,
        ),
        (
            |c| matches!(c, AtomConstraint::MulticenterValence(_)),
            AtomConstraint::MulticenterValence(ValueAst::Lit(0)),
            &cfg.multicenter_valence_mode,
        ),
    ];
    for (pred, default, mode) in numeric_defaults {
        if !matches!(mode, NumericMode::Zero) {
            continue;
        }
        if constraints.iter().any(pred) {
            continue;
        }
        constraints.push(default.clone());
    }
    let aromatic_default = match cfg.aromatic_valence_mode {
        AromaticValenceMode::NotAromatic => Some(AromaticValenceConstraint::NotAromatic),
        AromaticValenceMode::Aromatic => {
            Some(AromaticValenceConstraint::Value(ValueAst::Undetermined))
        }
        AromaticValenceMode::Required => None,
    };
    if let Some(default) = aromatic_default {
        if !constraints.iter().any(|c| matches!(c, AtomConstraint::AromaticValence(_))) {
            constraints.push(AtomConstraint::AromaticValence(default));
        }
    }
}

pub(crate) fn release_atom_constraints(constraints: &mut Vec<AtomConstraint>, cfg: &AtomAstConfig) {
    let strip_aromatic_undetermined =
        matches!(cfg.aromatic_valence_mode, AromaticValenceMode::Aromatic);
    let strip_zero_valence = matches!(cfg.valence_mode, NumericMode::Zero);
    let strip_zero_donated = matches!(cfg.donated_pairs_mode, NumericMode::Zero);
    let strip_zero_accepted = matches!(cfg.accepted_pairs_mode, NumericMode::Zero);
    let strip_zero_multicenter = matches!(cfg.multicenter_valence_mode, NumericMode::Zero);
    constraints.retain(|c| {
        if is_undetermined_numeric(c) {
            return false;
        }
        match c {
            AtomConstraint::Valence(ValueAst::Lit(0)) if strip_zero_valence => false,
            AtomConstraint::DonatedPairs(ValueAst::Lit(0)) if strip_zero_donated => false,
            AtomConstraint::AcceptedPairs(ValueAst::Lit(0)) if strip_zero_accepted => false,
            AtomConstraint::MulticenterValence(ValueAst::Lit(0)) if strip_zero_multicenter => {
                false
            }
            AtomConstraint::AromaticValence(AromaticValenceConstraint::NotAromatic) => false,
            AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(
                ValueAst::Undetermined,
            )) if strip_aromatic_undetermined => false,
            _ => true,
        }
    });
}

fn is_undetermined_numeric(c: &AtomConstraint) -> bool {
    matches!(
        c,
        AtomConstraint::Valence(ValueAst::Undetermined)
            | AtomConstraint::DonatedPairs(ValueAst::Undetermined)
            | AtomConstraint::AcceptedPairs(ValueAst::Undetermined)
            | AtomConstraint::MulticenterValence(ValueAst::Undetermined)
            | AtomConstraint::Degree(ValueAst::Undetermined)
            | AtomConstraint::Connectivity(ValueAst::Undetermined)
            | AtomConstraint::TotalHCount(ValueAst::Undetermined)
            | AtomConstraint::RingCount(ValueAst::Undetermined)
            | AtomConstraint::RingSize(ValueAst::Undetermined)
    )
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BondPattern {
    pub ast: BondAst,
    pub constraints: Vec<BondConstraint>,
}

impl BondPattern {
    pub fn new(ast: BondAst) -> Self {
        Self { ast, constraints: Vec::new() }
    }

    pub fn with_constraints(ast: BondAst, constraints: Vec<BondConstraint>) -> Self {
        Self { ast, constraints }
    }
}

#[derive(Debug)]
struct MoleculeMoleculePatternInner {
    ast: MoleculeAst,
}

#[derive(Clone, Debug)]
pub struct MoleculePattern(Arc<MoleculeMoleculePatternInner>);

impl MoleculePattern {
    pub fn new(ast: MoleculeAst) -> Self {
        Self(Arc::new(MoleculeMoleculePatternInner { ast }))
    }

    pub fn ast(&self) -> &MoleculeAst {
        &self.0.ast
    }
}

impl PartialEq for MoleculePattern {
    fn eq(&self, other: &Self) -> bool {
        self.0.ast == other.0.ast
    }
}

impl Eq for MoleculePattern {}

#[cfg(test)]
mod tests {
    use umol_shared::atom_ast::ElementAst;
    use umol_shared::element::Element;
    use umol_shared::value_ast::ValueAst;

    use super::*;

    #[test]
    fn test_molecule_pattern_new() {
        let ast = MoleculeAst::new(
            vec![AtomAst::new(ElementAst::Undetermined)],
            vec![], vec![], vec![], vec![], vec![], vec![],
        );
        let pattern = MoleculePattern::new(ast);
        assert_eq!(pattern.ast().atoms().count(), 1);
    }

    #[test]
    fn test_atom_pattern_new() {
        let pattern = AtomPattern::new(AtomAst::from_element(Element::C));
        assert_eq!(pattern.ast.element, ElementAst::Lit(Element::C));
        assert!(pattern.constraints.is_empty());
    }

    #[test]
    fn test_atom_pattern_with_constraints() {
        let pattern = AtomPattern::with_constraints(
            AtomAst::from_element(Element::C),
            vec![AtomConstraint::Valence(ValueAst::Lit(4))],
        );
        assert_eq!(pattern.constraints.len(), 1);
        assert_eq!(
            pattern.constraints[0],
            AtomConstraint::Valence(ValueAst::Lit(4)),
        );
    }

    #[test]
    fn test_bond_pattern_new() {
        let pattern = BondPattern::new(BondAst::from_order(1));
        assert_eq!(pattern.ast.order, ValueAst::Lit(1));
        assert!(pattern.constraints.is_empty());
    }

    #[test]
    fn test_bond_pattern_with_constraints() {
        let pattern = BondPattern::with_constraints(
            BondAst::from_order(2),
            vec![BondConstraint::RingBond],
        );
        assert_eq!(pattern.constraints.len(), 1);
        assert_eq!(pattern.constraints[0], BondConstraint::RingBond);
    }
}
