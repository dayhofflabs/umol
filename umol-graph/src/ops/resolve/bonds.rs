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

    /// Fills charge and spin defaults. Uses maximum multiplicity for partially defined
    /// spin states, otherwise uses closed-shell singlet.
    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), BondsContradiction>, BondsError> {
        ast.map_bonds(|mut bond| {
            if matches!(bond.charge, ValueAst::Undetermined) {
                bond.charge = ValueAst::Lit(0);
            }
            if bond.spin.is_undetermined() {
                bond.spin = SpinStateAst::closed_shell();
            } else {
                bond.spin.high_spin_complete();
            }
            bond
        });
        Ok(Solution::Determined(()))
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::BondId;
    use umol_ast::mol_dsl;

    use super::*;

    #[rstest]
    #[case::undetermined(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [{:atoms [0 1] :type "1"}]}"#), ValueAst::Lit(0), SpinStateAst::closed_shell())]
    #[case::spin_undetermined(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [{:atoms [0 1] :type "1#c+"}]}"#), ValueAst::Lit(1), SpinStateAst::closed_shell())]
    #[case::charge_undetermined(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [{:atoms [0 1] :type "1#u2#s1"}]}"#), ValueAst::Lit(0), SpinStateAst::from((2_u8, 1_u8)))]
    #[case::unpaired_undetermined(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [{:atoms [0 1] :type "1#s3"}]}"#), ValueAst::Lit(0), SpinStateAst::from((2_u8, 3_u8)))]
    #[case::multiplicity_undetermined(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [{:atoms [0 1] :type "1#u2"}]}"#), ValueAst::Lit(0), SpinStateAst::from((2_u8, 3_u8)))]
    fn test_bonds_resolver_resolve(
        #[case] mut mol: MoleculeAst,
        #[case] charge: ValueAst,
        #[case] spin: SpinStateAst,
    ) {
        BondsResolver::new().resolve(&mut mol).unwrap();
        assert_eq!(mol[BondId(0)].charge, charge);
        assert_eq!(mol[BondId(0)].spin, spin);
    }
}
