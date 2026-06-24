//! Multicenter-bond resolver: fills `MulticenterBondAst` charge / spin
//! defaults parallel to `BondsResolver`.

use thiserror::Error;
use umol_ast::ast::{Lattice, MoleculeAst, SpinStateAst, ValueAst};

use umol_utils::solution::Solution;

#[derive(Clone, Debug, Default)]
pub struct MulticenterBondsResolver;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum MulticenterBondsContradiction {}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum MulticenterBondsError {}

impl MulticenterBondsResolver {
    pub fn new() -> Self {
        Self
    }

    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), MulticenterBondsContradiction>, MulticenterBondsError> {
        for bond in ast.multicenter_bonds_mut() {
            if matches!(bond.charge, ValueAst::Undetermined) {
                bond.charge = ValueAst::Lit(0);
            }
            if bond.spin.is_undetermined() {
                bond.spin = SpinStateAst::closed_shell();
            }
        }
        Ok(Solution::Determined(()))
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{
        AtomAst, AtomId, Constraints, MoleculeAst, MulticenterBondAst, MulticenterBondId,
        SpinStateAst, ValueAst,
    };
    use umol_chem::element::Element;

    use super::*;

    fn one_multicenter(charge: ValueAst, spin: SpinStateAst) -> MoleculeAst {
        let bond =
            MulticenterBondAst::from_counts(vec![1, 0, 1])
                .with_charge(charge)
                .with_spin(spin);
        MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::B),
                AtomAst::from_element(Element::H),
                AtomAst::from_element(Element::B),
            ],
            vec![],
            vec![],
            vec![],
            vec![(vec![AtomId(0), AtomId(1), AtomId(2)], bond)],
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_multicenter_bonds_resolver_fills_undetermined() {
        let mut ast = one_multicenter(ValueAst::Undetermined, SpinStateAst::default());
        MulticenterBondsResolver::new().resolve(&mut ast).unwrap();
        let bond = &ast[MulticenterBondId(0)];
        assert_eq!(bond.charge, ValueAst::Lit(0));
        assert_eq!(bond.spin, SpinStateAst::closed_shell());
    }

    #[rstest]
    fn test_multicenter_bonds_resolver_preserves_existing_charge() {
        let mut ast = one_multicenter(ValueAst::Lit(-1), SpinStateAst::default());
        MulticenterBondsResolver::new().resolve(&mut ast).unwrap();
        let bond = &ast[MulticenterBondId(0)];
        assert_eq!(bond.charge, ValueAst::Lit(-1));
    }
}
