//! Multicenter-bond resolver: fills `MulticenterBondForm` charge and unpaired-electron
//! defaults parallel to `BondsResolver`.

use thiserror::Error;
use umol_graph_ir::ir::{
    AtomConstraintKey, Edits, IncidenceConstraintContradiction, IncidenceConstraintValidator,
    Lattice, Molecule, MulticenterBondHandle, MulticenterBondUpdate, NumForm, TransactionError,
    UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

#[derive(Clone, Debug, Default)]
pub struct MulticenterBondsResolver;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum MulticenterBondsContradiction {
    #[error(transparent)]
    Constraint(#[from] IncidenceConstraintContradiction),
}

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
    pub fn plan(&self, ast: &Molecule) -> Solution<Edits, MulticenterBondsContradiction> {
        for atom in ast.atoms().ids() {
            match IncidenceConstraintValidator
                .validate_molecule_atom_constraint(ast, atom, AtomConstraintKey::MulticenterValence)
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

        let mut edits = Edits::new();
        for bond_id in ast.multicenter_bonds().ids() {
            let bond = ast.multicenter_bond(bond_id).attributes;
            let mut selected_unpaired_electrons = bond.unpaired_electrons.clone();
            let mut update = MulticenterBondUpdate::default();
            if matches!(bond.charge, NumForm::Undetermined) {
                update.charge = Some(NumForm::Lit(0));
            }
            if selected_unpaired_electrons.is_undetermined() {
                selected_unpaired_electrons = UnpairedElectronsForm::closed_shell();
            } else {
                selected_unpaired_electrons.high_spin_complete();
            }
            update.unpaired_electrons = bond
                .unpaired_electrons
                .difference_to(&selected_unpaired_electrons);
            edits.update_multicenter_bond(MulticenterBondHandle::Id(bond_id), bond, &update);
        }
        Solution::Determined(edits)
    }

    /// Plan and atomically apply multicenter-bond defaults.
    pub fn resolve(
        &self,
        ast: &mut Molecule,
    ) -> Result<Solution<(), MulticenterBondsContradiction>, MulticenterBondsError> {
        let edits = match self.plan(ast) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let mut editor = ast.edit();
        editor.transact(edits)?;
        *ast = editor.build();
        Ok(Solution::Determined(()))
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{
        AtomConstraintForm, AtomId, Edit, Edits, MulticenterBondFieldChange, MulticenterBondId,
        MulticenterValenceForm,
    };
    use umol_graph_ir::mol_dsl;

    use super::*;

    #[rstest]
    #[case::undetermined(
        mol_dsl!(r#"{:atoms ["B" "H" "B"]
                       :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]"}]}"#),
        Edits::from_iter([
            Edit::ModifyMulticenterBondField {
                id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                change: MulticenterBondFieldChange::Charge {
                    old: NumForm::Undetermined,
                    new: NumForm::Lit(0),
                },
            },
            Edit::ModifyMulticenterBondField {
                id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                change: MulticenterBondFieldChange::UnpairedElectrons {
                    old: UnpairedElectronsForm::default(),
                    new: UnpairedElectronsForm::closed_shell(),
                },
            },
        ])
    )]
    #[case::partial_unpaired_electrons(
        mol_dsl!(r#"{:atoms ["B" "H" "B"]
                       :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]#c-#s3"}]}"#),
        Edits::from_iter([Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            change: MulticenterBondFieldChange::UnpairedElectrons {
                old: UnpairedElectronsForm {
                    count: NumForm::Undetermined,
                    multiplicity: NumForm::Lit(3),
                },
                new: UnpairedElectronsForm::from((2_u8, 3_u8)),
            },
        }])
    )]
    fn test_multicenter_bonds_resolver_plan(#[case] molecule: Molecule, #[case] expected: Edits) {
        assert_eq!(
            MulticenterBondsResolver::new().plan(&molecule),
            Solution::Determined(expected)
        );
    }

    #[rstest]
    #[case::determined(mol_dsl!(r#"{:atoms ["B" "H" "B"]
        :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]#c-#u2#s1"}]}"#))]
    fn test_multicenter_bonds_resolver_plan_identity(#[case] molecule: Molecule) {
        assert_eq!(
            MulticenterBondsResolver::new().plan(&molecule),
            Solution::Determined(Edits::new())
        );
    }

    #[rstest]
    #[case::contradictory(
        mol_dsl!(r#"{:atoms ["C#m1"]}"#),
        Solution::Contradictory(MulticenterBondsContradiction::Constraint(
            IncidenceConstraintContradiction::Atom {
                atom: AtomId(0),
                constraint: AtomConstraintForm::multicenter_valence(
                    MulticenterValenceForm::multicenter(1),
                ),
            },
        )),
    )]
    #[case::underdetermined(
        mol_dsl!(r#"{:atoms ["C#m1" "C" "C"]
                       :multicenter-bonds [{:atoms [0 1 2] :type "*"}]}"#),
        Solution::Underdetermined(Edits::new()),
    )]
    #[case::vacuous(
        mol_dsl!(r#"{:atoms ["C#m*"]}"#),
        Solution::Determined(Edits::new()),
    )]
    fn test_multicenter_bonds_resolver_plan_constraints(
        #[case] molecule: Molecule,
        #[case] expected: Solution<Edits, MulticenterBondsContradiction>,
    ) {
        assert_eq!(MulticenterBondsResolver::new().plan(&molecule), expected);
    }

    #[rstest]
    #[case::partial_unpaired_electrons(
        mol_dsl!(r#"{:atoms ["B" "H" "B"]
                       :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]#s3"}]}"#),
        mol_dsl!(r#"{:atoms ["B" "H" "B"]
                       :multicenter-bonds [{:atoms [0 1 2] :type "[1, 0, 1]#c0#u2#s3"}]}"#)
    )]
    fn test_multicenter_bonds_resolver_resolve(
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
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
        let Solution::Determined(edits) = MulticenterBondsResolver::new().plan(&molecule) else {
            panic!("fixture must produce a determined edit plan");
        };
        molecule
            .multicenter_bond_mut(MulticenterBondId(1))
            .attributes
            .charge = NumForm::Lit(9);
        let expected = molecule.clone();
        let mut editor = molecule.edit();
        assert_eq!(
            editor.transact(edits),
            Err(TransactionError::OldStateMismatch)
        );
        assert_eq!(editor.build(), expected);
    }
}
