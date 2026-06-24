//! Covalent-bond resolver: fills `BondAst` charge / spin defaults on bonds
//! whose corresponding inherent slots are still `Undetermined` after the
//! valence and aromaticity passes.

use thiserror::Error;
use umol_ast::ast::{Lattice, MoleculeAst, SpinStateAst, ValueAst};
use umol_utils::solution::Solution;

#[derive(Clone, Debug, Default)]
pub struct BondsResolver;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum BondsContradiction {}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum BondsError {}

impl BondsResolver {
    pub fn new() -> Self {
        Self
    }

    /// Fills `charge → Lit(0)` for any covalent bond whose `charge` is
    /// `Undetermined`, and `spin → closed_shell` for any bond whose spin is
    /// fully `Undetermined`. Existing literals / partially-constrained spin
    /// states pass through unchanged.
    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), BondsContradiction>, BondsError> {
        for bond in ast.bonds_mut() {
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
        AtomAst, AtomId, BondAst, BondId, Constraints, MoleculeAst, SpinStateAst, ValueAst,
    };
    use umol_chem::element::Element;

    use super::*;

    fn one_bond(charge: ValueAst, spin: SpinStateAst) -> MoleculeAst {
        let mut bond = BondAst::from_order(1);
        bond.charge = charge;
        bond.spin = spin;
        MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomId(0), AtomId(1), bond)],
            vec![],
            vec![],
            vec![],
            vec![],
            Vec::new(),
            Vec::new(),
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_bonds_resolver_fills_undetermined_charge_and_spin() {
        let mut ast = one_bond(ValueAst::Undetermined, SpinStateAst::default());
        BondsResolver::new().resolve(&mut ast).unwrap();
        let bond = ast.bond(BondId(0)).ast;
        assert_eq!(bond.charge, ValueAst::Lit(0));
        assert_eq!(bond.spin, SpinStateAst::closed_shell());
    }

    #[rstest]
    fn test_bonds_resolver_preserves_existing_charge() {
        let mut ast = one_bond(ValueAst::Lit(1), SpinStateAst::default());
        BondsResolver::new().resolve(&mut ast).unwrap();
        let bond = ast.bond(BondId(0)).ast;
        assert_eq!(bond.charge, ValueAst::Lit(1));
    }

    #[rstest]
    fn test_bonds_resolver_preserves_partial_spin() {
        let partial = SpinStateAst {
            unpaired: ValueAst::Lit(2),
            multiplicity: ValueAst::Undetermined,
        };
        let mut ast = one_bond(ValueAst::Undetermined, partial.clone());
        BondsResolver::new().resolve(&mut ast).unwrap();
        let bond = ast.bond(BondId(0)).ast;
        assert_eq!(bond.spin, partial);
    }
}
