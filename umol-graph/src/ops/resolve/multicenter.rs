//! Multicenter-bond resolver: fills `MulticenterBondForm` charge and unpaired-electron
//! defaults parallel to `BondsResolver`.

use thiserror::Error;
use umol_graph_ir::ir::{
    AtomConstraintKey, Edits, Lattice, Molecule, MulticenterBondHandle, Normalize, NumForm,
    UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

use crate::ops::validate::{
    DerivedKind, IncidenceConstraintInvariantsContradiction, IncidenceConstraintInvariantsValidator,
};

/// Resolves multicenter-bond defaults from the molecule's evidence.
///
/// # Semantic properties
///
/// Plans apply to their input and preserve molecule integrity. Consuming and in-place
/// operations agree on outcomes and determined molecules for the same input.
#[derive(Clone, Debug, Default)]
pub struct MulticenterBondsResolver;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum MulticenterBondsContradiction {
    #[error(transparent)]
    Constraint(#[from] IncidenceConstraintInvariantsContradiction),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum MulticenterBondsError {}

impl MulticenterBondsResolver {
    pub fn new() -> Self {
        Self
    }

    /// Construct charge and unpaired-electron default edits without mutating `molecule`.
    pub fn plan_resolve(
        &self,
        molecule: &Molecule,
    ) -> Solution<Edits, MulticenterBondsContradiction> {
        for atom in molecule.atoms().ids() {
            match IncidenceConstraintInvariantsValidator
                .validate_molecule_atom_constraint(
                    molecule,
                    atom,
                    AtomConstraintKey::MulticenterValence,
                    DerivedKind::DerivedComplete,
                )
                .expect("atom id came from the molecule atom store")
            {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => {
                    return Solution::Underdetermined(Edits::new());
                }
                Solution::Contradictory(contradiction) => {
                    return Solution::Contradictory(contradiction.into());
                }
            }
        }

        let mut charges = Vec::new();
        let mut unpaired_electrons = Vec::new();
        for bond_id in molecule.multicenter_bonds().ids() {
            let bond = molecule.multicenter_bond(bond_id).attributes();
            let mut selected_unpaired_electrons = bond.unpaired_electrons.clone();
            if matches!(bond.charge, NumForm::Undetermined) {
                charges.push((
                    MulticenterBondHandle::Id(bond_id),
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
                    MulticenterBondHandle::Id(bond_id),
                    bond.unpaired_electrons.clone(),
                    selected_unpaired_electrons,
                ));
            }
        }
        let mut edits = Edits::new();
        if !charges.is_empty() {
            edits.modify_multicenter_bond_charges(charges);
        }
        if !unpaired_electrons.is_empty() {
            edits.modify_multicenter_bond_unpaired_electrons(unpaired_electrons);
        }
        Solution::Determined(edits)
    }

    /// Resolves the molecule, consuming its storage without a recovery journal.
    ///
    /// Uses the same plan as [MulticenterBondsResolver::resolve_into]. A determined result has checked
    /// molecule integrity. Rejection drops the input.
    pub fn resolve(
        &self,
        molecule: Molecule,
    ) -> Result<Solution<Molecule, MulticenterBondsContradiction, ()>, MulticenterBondsError> {
        let edits = match self.plan_resolve(&molecule) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("multicenter-bond resolve plan applies to its input");
        let molecule = editor
            .finish()
            .expect("multicenter-bond resolve plan preserves molecule integrity");
        Ok(Solution::Determined(molecule))
    }

    /// Plans and applies multicenter-bond defaults without a recovery journal.
    ///
    /// # Semantic properties
    ///
    /// Planning finishes before mutation. Only Determined changes the molecule;
    /// every other outcome leaves it unchanged.
    pub fn resolve_into(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), MulticenterBondsContradiction>, MulticenterBondsError> {
        let edits = match self.plan_resolve(molecule) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        if !edits.is_empty() {
            molecule.mutate(edits);
        }
        Ok(Solution::Determined(()))
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{
        AtomConstraintForm, AtomId, Edit, Edits, MoleculeApplyError, MulticenterBondId,
        MulticenterValenceForm, RollbackMethod, TransactionError,
    };
    use umol_graph_ir::mol_dsl;

    use super::*;

    #[rstest]
    #[case::undetermined(
        mol_dsl!(
            r#"{:atoms ["B" "H" "B"]
                       :multicenter-bonds [{:atoms [0 1 2] :attrs "[1, 0, 1]"}]}"#
        ),
        Edits::from_iter([
            Edit::ModifyMulticenterBondCharges {
                changes: vec![(
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    NumForm::Undetermined,
                    NumForm::Lit(0),
                )],
            },
            Edit::ModifyMulticenterBondUnpairedElectrons {
                changes: vec![(
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    UnpairedElectronsForm::default(),
                    UnpairedElectronsForm::closed_shell(),
                )],
            },
        ]),
    )]
    #[case::partial_unpaired_electrons(
        mol_dsl!(
            r#"{:atoms ["B" "H" "B"]
                       :multicenter-bonds [{:atoms [0 1 2] :attrs "[1, 0, 1]#c-#s3"}]}"#
        ),
        Edits::from_iter([Edit::ModifyMulticenterBondUnpairedElectrons {
            changes: vec![(
                MulticenterBondHandle::Id(MulticenterBondId(0)),
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
            r#"{:atoms ["B" "H" "B" "H" "B" "H" "B"]
            :multicenter-bonds [
                {:atoms [0 1 2] :attrs "[1, 0, 1]"}
                {:atoms [2 3 4] :attrs "[1, 0, 1]#c-#s3"}
                {:atoms [4 5 6] :attrs "[1, 0, 1]#u0#s"}]}"#
        ),
        Edits::from_iter([
            Edit::ModifyMulticenterBondCharges {
                changes: vec![
                    (
                        MulticenterBondHandle::Id(MulticenterBondId(0)),
                        NumForm::Undetermined,
                        NumForm::Lit(0),
                    ),
                    (
                        MulticenterBondHandle::Id(MulticenterBondId(2)),
                        NumForm::Undetermined,
                        NumForm::Lit(0),
                    ),
                ],
            },
            Edit::ModifyMulticenterBondUnpairedElectrons {
                changes: vec![
                    (
                        MulticenterBondHandle::Id(MulticenterBondId(0)),
                        UnpairedElectronsForm::default(),
                        UnpairedElectronsForm::closed_shell(),
                    ),
                    (
                        MulticenterBondHandle::Id(MulticenterBondId(1)),
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
    fn test_multicenter_bonds_resolver_plan_resolve(
        #[case] molecule: Molecule,
        #[case] expected: Edits,
    ) {
        assert_eq!(
            MulticenterBondsResolver::new().plan_resolve(&molecule),
            Solution::Determined(expected)
        );
    }

    #[rstest]
    #[case::determined(mol_dsl!(r#"{:atoms ["B" "H" "B"]
        :multicenter-bonds [{:atoms [0 1 2] :attrs "[1, 0, 1]#c-#u2#s1"}]}"#))]
    fn test_multicenter_bonds_resolver_plan_resolve_identity(#[case] molecule: Molecule) {
        assert_eq!(
            MulticenterBondsResolver::new().plan_resolve(&molecule),
            Solution::Determined(Edits::new())
        );
    }

    #[rstest]
    #[case::contradictory(
        mol_dsl!(r#"{:atoms ["C#m1"]}"#),
        Solution::Contradictory(MulticenterBondsContradiction::Constraint(
            IncidenceConstraintInvariantsContradiction::Atom {
                atom: AtomId(0),
                constraint: AtomConstraintForm::multicenter_valence(
                    MulticenterValenceForm::multicenter(1),
                ),
            },
        )),
    )]
    #[case::underdetermined(
        mol_dsl!(r#"{:atoms ["C#m1" "C" "C"]
                       :multicenter-bonds [{:atoms [0 1 2] :attrs "*"}]}"#),
        Solution::Underdetermined(Edits::new()),
    )]
    #[case::vacuous(
        mol_dsl!(r#"{:atoms ["C#m*"]}"#),
        Solution::Determined(Edits::new()),
    )]
    fn test_multicenter_bonds_resolver_plan_resolve_constraints(
        #[case] mut molecule: Molecule,
        #[case] expected: Solution<Edits, MulticenterBondsContradiction>,
    ) {
        let original = molecule.clone();
        let resolver = MulticenterBondsResolver::new();
        assert_eq!(resolver.plan_resolve(&molecule), expected);
        assert_eq!(
            resolver.resolve_into(&mut molecule),
            Ok(expected.map(|_| ()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::partial_unpaired_electrons(
        mol_dsl!(r#"{:atoms ["B" "H" "B"]
                       :multicenter-bonds [{:atoms [0 1 2] :attrs "[1, 0, 1]#s3"}]}"#),
        mol_dsl!(r#"{:atoms ["B" "H" "B"]
                       :multicenter-bonds [{:atoms [0 1 2] :attrs "[1, 0, 1]#c0#u2#s3"}]}"#)
    )]
    fn test_multicenter_bonds_resolver_resolve(
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        assert_eq!(
            MulticenterBondsResolver::new().resolve(molecule.clone()),
            Ok(Solution::Determined(expected.clone()))
        );
        assert_eq!(
            MulticenterBondsResolver::new().resolve_into(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    fn test_multicenter_bonds_resolver_plan_resolve_stale(
        #[values(RollbackMethod::Clone, RollbackMethod::Journal)] rollback: RollbackMethod,
    ) {
        let mut molecule = mol_dsl!(
            r#"{
            :atoms ["B" "H" "B" "H" "B"]
            :multicenter-bonds [
                {:atoms [0 1 2] :attrs "[1, 0, 1]"}
                {:atoms [2 3 4] :attrs "[1, 0, 1]"}
            ]
        }"#
        );
        let Solution::Determined(edits) = MulticenterBondsResolver::new().plan_resolve(&molecule)
        else {
            panic!("fixture must produce a determined edit plan");
        };
        molecule
            .multicenter_bond_mut(MulticenterBondId(1))
            .attributes_mut()
            .charge = NumForm::Lit(9);
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
