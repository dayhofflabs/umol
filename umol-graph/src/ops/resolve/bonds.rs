//! Covalent-bond resolver: fills `BondForm` charge and unpaired-electron defaults on bonds
//! whose corresponding inherent fields are still `Undetermined` after the
//! valence and aromaticity passes.

use std::mem;

use thiserror::Error;
use umol_graph_ir::ir::{
    BondHandle, Edits, Lattice, Molecule, Normalize, NumForm, UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

/// Resolves localized-bond defaults from the molecule's evidence.
///
/// # Semantic properties
///
/// Plans apply to their input and preserve molecule integrity. Consuming and in-place
/// operations agree on outcomes and determined molecules for the same input.
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

    /// Construct charge and unpaired-electron default edits without mutating `molecule`.
    pub fn plan_resolve(&self, molecule: &Molecule) -> Edits {
        let mut charges = Vec::new();
        let mut unpaired_electrons = Vec::new();
        for bond_id in molecule.bonds().ids() {
            let bond = molecule.bond(bond_id).attributes();
            let mut selected_unpaired_electrons = bond.unpaired_electrons.clone();
            if matches!(bond.charge, NumForm::Undetermined) {
                charges.push((
                    BondHandle::Id(bond_id),
                    NumForm::Undetermined,
                    NumForm::Lit(0),
                ));
            }
            if selected_unpaired_electrons.is_undetermined() {
                selected_unpaired_electrons = UnpairedElectronsForm::closed_shell();
            } else {
                selected_unpaired_electrons.high_spin_complete();
            }
            if !bond
                .unpaired_electrons
                .normalized_eq(&selected_unpaired_electrons)
            {
                unpaired_electrons.push((
                    BondHandle::Id(bond_id),
                    bond.unpaired_electrons.clone(),
                    selected_unpaired_electrons,
                ));
            }
        }
        let mut edits = Edits::new();
        if !charges.is_empty() {
            edits.modify_bond_charges(charges);
        }
        if !unpaired_electrons.is_empty() {
            edits.modify_bond_unpaired_electrons(unpaired_electrons);
        }
        edits
    }

    /// Resolves the molecule, consuming its storage without a recovery journal.
    ///
    /// Uses the same plan as [BondsResolver::resolve_into]. A determined result has checked
    /// molecule integrity. Rejection drops the input.
    pub fn resolve(
        &self,
        molecule: Molecule,
    ) -> Result<Solution<Molecule, BondsContradiction, ()>, BondsError> {
        let edits = self.plan_resolve(&molecule);
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("bond resolve plan applies to its input");
        let molecule = editor
            .finish()
            .expect("bond resolve plan preserves molecule integrity");
        Ok(Solution::Determined(molecule))
    }

    /// Plans and applies localized-bond defaults without a recovery journal.
    ///
    /// # Semantic properties
    ///
    /// Planning finishes before mutation. Only Determined changes the molecule;
    /// every other outcome leaves it unchanged.
    pub fn resolve_into(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), BondsContradiction>, BondsError> {
        let edits = self.plan_resolve(molecule);
        if !edits.is_empty() {
            *molecule = mem::take(molecule)
                .apply(edits)
                .expect("bond resolve plan applies and preserves molecule integrity");
        }
        Ok(Solution::Determined(()))
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{
        BondId, Edit, Edits, MoleculeApplyError, RollbackMethod, TransactionError,
    };
    use umol_graph_ir::mol_dsl;

    use super::*;

    #[rstest]
    #[case::undetermined(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        Edits::from_iter([
            Edit::ModifyBondCharges {
                changes: vec![(
                    BondHandle::Id(BondId(0)),
                    NumForm::Undetermined,
                    NumForm::Lit(0),
                )],
            },
            Edit::ModifyBondUnpairedElectrons {
                changes: vec![(
                    BondHandle::Id(BondId(0)),
                    UnpairedElectronsForm::default(),
                    UnpairedElectronsForm::closed_shell(),
                )],
            },
        ]),
    )]
    #[case::partial_unpaired_electrons(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c+#s3"]]}"#),
        Edits::from_iter([Edit::ModifyBondUnpairedElectrons {
            changes: vec![(
                BondHandle::Id(BondId(0)),
                UnpairedElectronsForm {
                    count: NumForm::Undetermined,
                    multiplicity: NumForm::Lit(3),
                },
                UnpairedElectronsForm::from((2_u8, 3_u8)),
            )],
        }]),
    )]
    #[case::independent_fields(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1#c+#s3"] [2 3 "1#u0#s"]]}"#
        ),
        Edits::from_iter([
            Edit::ModifyBondCharges {
                changes: vec![
                    (
                        BondHandle::Id(BondId(0)),
                        NumForm::Undetermined,
                        NumForm::Lit(0),
                    ),
                    (
                        BondHandle::Id(BondId(2)),
                        NumForm::Undetermined,
                        NumForm::Lit(0),
                    ),
                ],
            },
            Edit::ModifyBondUnpairedElectrons {
                changes: vec![
                    (
                        BondHandle::Id(BondId(0)),
                        UnpairedElectronsForm::default(),
                        UnpairedElectronsForm::closed_shell(),
                    ),
                    (
                        BondHandle::Id(BondId(1)),
                        UnpairedElectronsForm {
                            count: NumForm::Undetermined,
                            multiplicity: NumForm::Lit(3),
                        },
                        UnpairedElectronsForm::from((2_u8, 3_u8)),
                    ),
                ],
            },
        ]),
    )]
    fn test_bonds_resolver_plan_resolve(#[case] molecule: Molecule, #[case] expected: Edits) {
        assert_eq!(BondsResolver::new().plan_resolve(&molecule), expected);
    }

    #[rstest]
    #[case::determined(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c+#u2#s1"]]}"#))]
    fn test_bonds_resolver_plan_resolve_identity(#[case] molecule: Molecule) {
        assert_eq!(BondsResolver::new().plan_resolve(&molecule), Edits::new());
    }

    #[rstest]
    #[case::partial_unpaired_electrons(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#s3"]]}"#),
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c0#u2#s3"]]}"#)
    )]
    fn test_bonds_resolver_resolve(#[case] mut molecule: Molecule, #[case] expected: Molecule) {
        assert_eq!(
            BondsResolver::new().resolve(molecule.clone()),
            Ok(Solution::Determined(expected.clone()))
        );
        assert_eq!(
            BondsResolver::new().resolve_into(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    fn test_bonds_resolver_plan_resolve_stale(
        #[values(RollbackMethod::Clone, RollbackMethod::Journal)] rollback: RollbackMethod,
    ) {
        let mut molecule = mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#);
        let edits = BondsResolver::new().plan_resolve(&molecule);
        molecule.bond_mut(BondId(1)).attributes_mut().charge = NumForm::Lit(9);
        let expected = molecule.clone();
        assert_eq!(
            molecule.clone().edit().apply(edits.clone()),
            Err(MoleculeApplyError::Transaction(
                TransactionError::OldStateMismatch
            ))
        );
        assert_eq!(
            molecule.transact([edits], rollback),
            Err(MoleculeApplyError::Transaction(
                TransactionError::OldStateMismatch
            ))
        );
        assert_eq!(molecule, expected);
    }
}
