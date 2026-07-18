//! Covalent-bond resolver: fills `BondAst` charge / spin defaults on bonds
//! whose corresponding inherent slots are still `Undetermined` after the
//! valence and aromaticity passes.

use thiserror::Error;
use umol_ast::ast::{
    BondHandle, BondUpdate, Edit, Lattice, MoleculeAst, SpinStateAst, TransactionError, ValueAst,
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

    /// Construct charge and spin default edits without mutating `ast`.
    pub fn plan(&self, ast: &MoleculeAst) -> Vec<Edit> {
        let mut edits = Vec::new();
        for bond_id in ast.bonds().ids() {
            let bond = ast.bond(bond_id).ast;
            let mut selected_spin = bond.spin.clone();
            let mut update = BondUpdate::default();
            if matches!(bond.charge, ValueAst::Undetermined) {
                update.charge = Some(ValueAst::Lit(0));
            }
            if selected_spin.is_undetermined() {
                selected_spin = SpinStateAst::closed_shell();
            } else {
                selected_spin.high_spin_complete();
            }
            update.spin = bond.spin.difference_to(&selected_spin);
            edits.extend(Edit::for_bond_update(
                BondHandle::Id(bond_id),
                bond,
                &update,
            ));
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
    use umol_ast::ast::{BondFieldChange, BondId};
    use umol_ast::mol_dsl;

    use super::*;

    #[rstest]
    #[case::undetermined(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        vec![
            Edit::ModifyBondField {
                id: BondHandle::Id(BondId(0)),
                change: BondFieldChange::Charge {
                    old: ValueAst::Undetermined,
                    new: ValueAst::Lit(0),
                },
            },
            Edit::ModifyBondField {
                id: BondHandle::Id(BondId(0)),
                change: BondFieldChange::Spin {
                    old: SpinStateAst::default(),
                    new: SpinStateAst::closed_shell(),
                },
            },
        ]
    )]
    #[case::partial_spin(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c+#s3"]]}"#),
        vec![Edit::ModifyBondField {
            id: BondHandle::Id(BondId(0)),
            change: BondFieldChange::Spin {
                old: SpinStateAst {
                    unpaired: ValueAst::Undetermined,
                    multiplicity: ValueAst::Lit(3),
                },
                new: SpinStateAst::from((2_u8, 3_u8)),
            },
        }]
    )]
    fn test_bonds_resolver_plan(#[case] molecule: MoleculeAst, #[case] expected: Vec<Edit>) {
        assert_eq!(BondsResolver::new().plan(&molecule), expected);
    }

    #[rstest]
    #[case::determined(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c+#u2#s1"]]}"#))]
    fn test_bonds_resolver_plan_identity(#[case] molecule: MoleculeAst) {
        assert_eq!(BondsResolver::new().plan(&molecule), Vec::new());
    }

    #[rstest]
    #[case::partial_spin(
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
        molecule.bond_mut(BondId(1)).ast.charge = ValueAst::Lit(9);
        let expected = molecule.clone();
        let mut editor = molecule.edit();
        assert_eq!(
            editor.transact(edits),
            Err(TransactionError::OldStateMismatch)
        );
        assert_eq!(editor.build(), expected);
    }
}
