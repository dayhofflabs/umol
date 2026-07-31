//! Multicenter-bond resolver: fills `MulticenterBondAst` charge and unpaired-electron
//! defaults parallel to `BondsResolver`.

use thiserror::Error;
use umol_ast::ast::{
    Edit, Lattice, MoleculeAst, MulticenterBondHandle, MulticenterBondUpdate, TransactionError,
    UnpairedElectronsAst, ValueAst,
};
use umol_utils::solution::Solution;

#[derive(Clone, Debug, Default)]
pub struct MulticenterBondsResolver;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum MulticenterBondsContradiction {}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum MulticenterBondsError {
    #[error(transparent)]
    Transaction(#[from] TransactionError),
}

impl MulticenterBondsResolver {
    pub fn new() -> Self {
        Self
    }

    /// Construct charge and unpaired-electron default edits without mutating `ast`.
    pub fn plan(&self, ast: &MoleculeAst) -> Vec<Edit> {
        let mut edits = Vec::new();
        for bond_id in ast.multicenter_bonds().ids() {
            let bond = ast.multicenter_bond(bond_id).ast;
            let mut selected_unpaired_electrons = bond.unpaired_electrons.clone();
            let mut update = MulticenterBondUpdate::default();
            if matches!(bond.charge, ValueAst::Undetermined) {
                update.charge = Some(ValueAst::Lit(0));
            }
            if selected_unpaired_electrons.is_undetermined() {
                selected_unpaired_electrons = UnpairedElectronsAst::closed_shell();
            } else {
                selected_unpaired_electrons.high_spin_complete();
            }
            update.unpaired_electrons = bond
                .unpaired_electrons
                .difference_to(&selected_unpaired_electrons);
            edits.extend(Edit::for_multicenter_bond_update(
                MulticenterBondHandle::Id(bond_id),
                bond,
                &update,
            ));
        }
        edits
    }

    /// Plan and atomically apply multicenter-bond defaults.
    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), MulticenterBondsContradiction>, MulticenterBondsError> {
        let edits = self.plan(ast);
        let mut editor = ast.edit();
        editor.transact(edits)?;
        *ast = editor.build();
        Ok(Solution::Determined(()))
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{MulticenterBondFieldChange, MulticenterBondId};
    use umol_ast::mol_dsl;

    use super::*;

    #[rstest]
    #[case::undetermined(
        mol_dsl!(r#"{:atoms ["B" "H" "B"]
                       :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]"}]}"#),
        vec![
            Edit::ModifyMulticenterBondField {
                id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                change: MulticenterBondFieldChange::Charge {
                    old: ValueAst::Undetermined,
                    new: ValueAst::Lit(0),
                },
            },
            Edit::ModifyMulticenterBondField {
                id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                change: MulticenterBondFieldChange::UnpairedElectrons {
                    old: UnpairedElectronsAst::default(),
                    new: UnpairedElectronsAst::closed_shell(),
                },
            },
        ]
    )]
    #[case::partial_unpaired_electrons(
        mol_dsl!(r#"{:atoms ["B" "H" "B"]
                       :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]#c-#s3"}]}"#),
        vec![Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            change: MulticenterBondFieldChange::UnpairedElectrons {
                old: UnpairedElectronsAst {
                    count: ValueAst::Undetermined,
                    multiplicity: ValueAst::Lit(3),
                },
                new: UnpairedElectronsAst::from((2_u8, 3_u8)),
            },
        }]
    )]
    fn test_multicenter_bonds_resolver_plan(
        #[case] molecule: MoleculeAst,
        #[case] expected: Vec<Edit>,
    ) {
        assert_eq!(MulticenterBondsResolver::new().plan(&molecule), expected);
    }

    #[rstest]
    #[case::determined(mol_dsl!(r#"{:atoms ["B" "H" "B"]
        :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]#c-#u2#s1"}]}"#))]
    fn test_multicenter_bonds_resolver_plan_identity(#[case] molecule: MoleculeAst) {
        assert_eq!(MulticenterBondsResolver::new().plan(&molecule), Vec::new());
    }

    #[rstest]
    #[case::partial_unpaired_electrons(
        mol_dsl!(r#"{:atoms ["B" "H" "B"]
                       :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]#s3"}]}"#),
        mol_dsl!(r#"{:atoms ["B" "H" "B"]
                       :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]#c0#u2#s3"}]}"#)
    )]
    fn test_multicenter_bonds_resolver_resolve(
        #[case] mut molecule: MoleculeAst,
        #[case] expected: MoleculeAst,
    ) {
        assert_eq!(
            MulticenterBondsResolver::new().resolve(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    fn test_multicenter_bonds_resolver_plan_stale() {
        let mut molecule = mol_dsl!(
            r#"{
            :atoms ["B" "H" "B" "H" "B"]
            :multicenter-bonds [
                {:atoms [0 1 2] :type "[1, 0, 1]"}
                {:atoms [2 3 4] :type "[1, 0, 1]"}
            ]
        }"#
        );
        let edits = MulticenterBondsResolver::new().plan(&molecule);
        molecule
            .multicenter_bond_mut(MulticenterBondId(1))
            .ast
            .charge = ValueAst::Lit(9);
        let expected = molecule.clone();
        let mut editor = molecule.edit();
        assert_eq!(
            editor.transact(edits),
            Err(TransactionError::OldStateMismatch)
        );
        assert_eq!(editor.build(), expected);
    }
}
