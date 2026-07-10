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
            } else {
                bond.spin.high_spin_complete();
            }
        }
        Ok(Solution::Determined(()))
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{MulticenterBondId, SpinStateAst, ValueAst};
    use umol_ast::mol_dsl;

    use super::*;

    #[rstest]
    #[case::undetermined(mol_dsl!(r#"{:atoms ["B" "H" "B"] :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]"}]}"#), ValueAst::Lit(0), SpinStateAst::closed_shell())]
    #[case::spin_undetermined(mol_dsl!(r#"{:atoms ["B" "H" "B"] :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]#c-"}]}"#), ValueAst::Lit(-1), SpinStateAst::closed_shell())]
    #[case::charge_undetermined(mol_dsl!(r#"{:atoms ["B" "H" "B"] :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]#u2#s1"}]}"#), ValueAst::Lit(0), SpinStateAst::from((2_u8, 1_u8)))]
    #[case::unpaired_undetermined(mol_dsl!(r#"{:atoms ["B" "H" "B"] :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]#s3"}]}"#), ValueAst::Lit(0), SpinStateAst::from((2_u8, 3_u8)))]
    #[case::multiplicity_undetermined(mol_dsl!(r#"{:atoms ["B" "H" "B"] :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]#u2"}]}"#), ValueAst::Lit(0), SpinStateAst::from((2_u8, 3_u8)))]
    fn test_multicenter_bonds_resolver_fills_undetermined_charge_and_spin(
        #[case] mut mol: MoleculeAst,
        #[case] charge: ValueAst,
        #[case] spin: SpinStateAst,
    ) {
        MulticenterBondsResolver::new().resolve(&mut mol).unwrap();
        assert_eq!(mol[MulticenterBondId(0)].charge, charge);
        assert_eq!(mol[MulticenterBondId(0)].spin, spin);
    }
}
