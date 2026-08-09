//! Covalent-bond resolver: fills `BondAst` charge and unpaired-electron defaults on bonds
//! whose corresponding inherent slots are still `Undetermined` after the
//! valence and aromaticity passes.

use thiserror::Error;
use umol_graph_ir::ir::{
    BondHandle, BondUpdate, Edits, Lattice, MoleculeAst, NumForm, TransactionError,
    UnpairedElectronsAst,
};
use umol_utils::solution::Solution;

#[derive(Clone, Debug, Default)]
pub struct BondsResolver;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum BondsContradiction {}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum BondsError {
    #[error(transparent)]
    Transaction(#[from] TransactionError),
}

impl BondsResolver {
    pub fn new() -> Self {
        Self
    }

    /// Construct charge and unpaired-electron default edits without mutating `ast`.
    pub fn plan(&self, ast: &MoleculeAst) -> Edits {
        let mut edits = Edits::new();
        for bond_id in ast.bonds().ids() {
            let bond = ast.bond(bond_id).ast;
            let mut selected_unpaired_electrons = bond.unpaired_electrons.clone();
            let mut update = BondUpdate::default();
            if matches!(bond.charge, NumForm::Undetermined) {
                update.charge = Some(NumForm::Lit(0));
            }
            if selected_unpaired_electrons.is_undetermined() {
                selected_unpaired_electrons = UnpairedElectronsAst::closed_shell();
            } else {
                selected_unpaired_electrons.high_spin_complete();
            }
            update.unpaired_electrons = bond
                .unpaired_electrons
                .difference_to(&selected_unpaired_electrons);
            edits.update_bond(BondHandle::Id(bond_id), bond, &update);
        }
        edits
    }

    /// Plan and atomically apply localized-bond defaults.
    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), BondsContradiction>, BondsError> {
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
    use umol_graph_ir::ir::{BondFieldChange, BondId, Edit, Edits};
    use umol_graph_ir::mol_dsl;

    use super::*;

    #[rstest]
    #[case::undetermined(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        Edits::from_iter([
            Edit::ModifyBondField {
                id: BondHandle::Id(BondId(0)),
                change: BondFieldChange::Charge {
                    old: NumForm::Undetermined,
                    new: NumForm::Lit(0),
                },
            },
            Edit::ModifyBondField {
                id: BondHandle::Id(BondId(0)),
                change: BondFieldChange::UnpairedElectrons {
                    old: UnpairedElectronsAst::default(),
                    new: UnpairedElectronsAst::closed_shell(),
                },
            },
        ])
    )]
    #[case::partial_unpaired_electrons(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c+#s3"]]}"#),
        Edits::from_iter([Edit::ModifyBondField {
            id: BondHandle::Id(BondId(0)),
            change: BondFieldChange::UnpairedElectrons {
                old: UnpairedElectronsAst {
                    count: NumForm::Undetermined,
                    multiplicity: NumForm::Lit(3),
                },
                new: UnpairedElectronsAst::from((2_u8, 3_u8)),
            },
        }])
    )]
    fn test_bonds_resolver_plan(#[case] molecule: MoleculeAst, #[case] expected: Edits) {
        assert_eq!(BondsResolver::new().plan(&molecule), expected);
    }

    #[rstest]
    #[case::determined(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c+#u2#s1"]]}"#))]
    fn test_bonds_resolver_plan_identity(#[case] molecule: MoleculeAst) {
        assert_eq!(BondsResolver::new().plan(&molecule), Edits::new());
    }

    #[rstest]
    #[case::partial_unpaired_electrons(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#s3"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c0#u2#s3"]]}"#)
    )]
    fn test_bonds_resolver_resolve(
        #[case] mut molecule: MoleculeAst,
        #[case] expected: MoleculeAst,
    ) {
        assert_eq!(
            BondsResolver::new().resolve(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    fn test_bonds_resolver_plan_stale() {
        let mut molecule = mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#);
        let edits = BondsResolver::new().plan(&molecule);
        molecule.bond_mut(BondId(1)).ast.charge = NumForm::Lit(9);
        let expected = molecule.clone();
        let mut editor = molecule.edit();
        assert_eq!(
            editor.transact(edits),
            Err(TransactionError::OldStateMismatch)
        );
        assert_eq!(editor.build(), expected);
    }
}
