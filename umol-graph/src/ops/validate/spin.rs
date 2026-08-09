//! Spin-state invariant validation for entity fields and molecule-level
//! unpaired-electron coupling targets. A complete literal unpaired-electron
//! count and multiplicity must form a physically valid [`SpinState`]; a pair
//! with either component still non-literal is underdetermined.

use thiserror::Error;
use umol_chem::error::SpinStateError;
use umol_chem::spin::SpinState;
use umol_graph_ir::ir::{
    AromaticSystemId, AsLit, AtomForm, AtomId, BondId, Constraint, Lattice, Molecule,
    MoleculeConstraint, MulticenterBondId, UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

#[derive(Clone, Copy, Debug, Default)]
pub struct SpinInvariantsValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SpinInvariantsContradiction {
    #[error("atom has invalid unpaired electrons: {error}")]
    Atom { error: SpinStateError },
    #[error("molecule atom {atom} has invalid unpaired electrons: {error}")]
    MoleculeAtom { atom: AtomId, error: SpinStateError },
    #[error("bond {bond} has invalid unpaired electrons: {error}")]
    Bond { bond: BondId, error: SpinStateError },
    #[error("aromatic system {system} has invalid unpaired electrons: {error}")]
    AromaticSystem {
        system: AromaticSystemId,
        error: SpinStateError,
    },
    #[error("multicenter bond {bond} has invalid unpaired electrons: {error}")]
    MulticenterBond {
        bond: MulticenterBondId,
        error: SpinStateError,
    },
    #[error("unpaired-electron coupling in constraint {constraint_index} is invalid: {error}")]
    UnpairedElectronCoupling {
        constraint_index: usize,
        error: SpinStateError,
    },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SpinInvariantsError {}

impl SpinInvariantsValidator {
    pub fn validate(
        &self,
        ast: &Molecule,
    ) -> Result<Solution<(), SpinInvariantsContradiction>, SpinInvariantsError> {
        let mut any_undetermined = false;

        for atom in ast.atoms().iter() {
            match validate_unpaired_electrons(atom.unpaired_electrons()) {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => any_undetermined = true,
                Solution::Contradictory(error) => {
                    return Ok(Solution::Contradictory(
                        SpinInvariantsContradiction::MoleculeAtom {
                            atom: atom.id,
                            error,
                        },
                    ));
                }
            }
        }
        for bond in ast.bonds().iter() {
            match validate_unpaired_electrons(bond.unpaired_electrons()) {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => any_undetermined = true,
                Solution::Contradictory(error) => {
                    return Ok(Solution::Contradictory(SpinInvariantsContradiction::Bond {
                        bond: bond.id,
                        error,
                    }));
                }
            }
        }
        for system in ast.aromatic_systems().iter() {
            match validate_unpaired_electrons(system.unpaired_electrons()) {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => any_undetermined = true,
                Solution::Contradictory(error) => {
                    return Ok(Solution::Contradictory(
                        SpinInvariantsContradiction::AromaticSystem {
                            system: system.id,
                            error,
                        },
                    ));
                }
            }
        }
        for bond in ast.multicenter_bonds().iter() {
            match validate_unpaired_electrons(bond.unpaired_electrons()) {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => any_undetermined = true,
                Solution::Contradictory(error) => {
                    return Ok(Solution::Contradictory(
                        SpinInvariantsContradiction::MulticenterBond {
                            bond: bond.id,
                            error,
                        },
                    ));
                }
            }
        }
        for (constraint_index, constraint) in ast.constraints().iter().enumerate() {
            match validate_couplings_in_constraint(constraint) {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => any_undetermined = true,
                Solution::Contradictory(error) => {
                    return Ok(Solution::Contradictory(
                        SpinInvariantsContradiction::UnpairedElectronCoupling {
                            constraint_index,
                            error,
                        },
                    ));
                }
            }
        }

        Ok(if any_undetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        })
    }

    pub fn validate_atom(
        &self,
        atom: &AtomForm,
    ) -> Result<Solution<(), SpinInvariantsContradiction>, SpinInvariantsError> {
        Ok(validate_unpaired_electrons(&atom.unpaired_electrons)
            .map_contradiction(|error| SpinInvariantsContradiction::Atom { error }))
    }
}

fn validate_couplings_in_constraint(constraint: &Constraint) -> Solution<(), SpinStateError> {
    match constraint {
        Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
            unpaired_electrons,
            ..
        }) => {
            if unpaired_electrons.is_undetermined() {
                return Solution::Determined(());
            }
            match validate_unpaired_electrons(unpaired_electrons) {
                Solution::Contradictory(error) => Solution::Contradictory(error),
                Solution::Determined(()) | Solution::Underdetermined(()) => {
                    Solution::Underdetermined(())
                }
            }
        }
        Constraint::And(children) | Constraint::Or(children) => {
            let mut any_undetermined = false;
            for child in children {
                match validate_couplings_in_constraint(child) {
                    Solution::Determined(()) => {}
                    Solution::Underdetermined(()) => any_undetermined = true,
                    Solution::Contradictory(error) => {
                        return Solution::Contradictory(error);
                    }
                }
            }
            if any_undetermined {
                Solution::Underdetermined(())
            } else {
                Solution::Determined(())
            }
        }
        Constraint::Not(child) => validate_couplings_in_constraint(child),
        _ => Solution::Determined(()),
    }
}

fn validate_unpaired_electrons(
    unpaired_electrons: &UnpairedElectronsForm,
) -> Solution<(), SpinStateError> {
    let Some(unpaired_electrons) = unpaired_electrons.as_lit() else {
        return Solution::Underdetermined(());
    };
    match SpinState::try_from(unpaired_electrons) {
        Ok(_) => Solution::Determined(()),
        Err(error) => Solution::Contradictory(error),
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use proptest::prelude::*;
    use rstest::rstest;
    use umol_chem::spin::SpinMultiplicity;
    use umol_graph_ir::ir::{
        AromaticSystemForm, BondForm, MoleculeEntries, MulticenterBondForm, NumForm,
    };

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::closed_shell((0_u8, 1_u8).into(), Solution::Determined(()))]
    #[case::doublet((1_u8, 2_u8).into(), Solution::Determined(()))]
    #[case::open_shell_singlet((2_u8, 1_u8).into(), Solution::Determined(()))]
    #[case::triplet((2_u8, 3_u8).into(), Solution::Determined(()))]
    #[case::count_undetermined(
        UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(1) },
        Solution::Underdetermined(()),
    )]
    #[case::multiplicity_undetermined(
        UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined },
        Solution::Underdetermined(()),
    )]
    #[case::count_negative(
        UnpairedElectronsForm { count: NumForm::Lit(-1), multiplicity: NumForm::Lit(1) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::UnpairedElectronsOutOfRange { count: -1 },
        }),
    )]
    #[case::count_above_u8(
        UnpairedElectronsForm { count: NumForm::Lit(256), multiplicity: NumForm::Lit(1) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::UnpairedElectronsOutOfRange { count: 256 },
        }),
    )]
    #[case::multiplicity_zero(
        UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Lit(0) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::MultiplicityOutOfRange { multiplicity: 0 },
        }),
    )]
    #[case::multiplicity_above_u8(
        UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Lit(256) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::MultiplicityOutOfRange { multiplicity: 256 },
        }),
    )]
    #[case::parity(
        UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Lit(2) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::Incompatible {
                unpaired_electrons: 2,
                multiplicity: SpinMultiplicity::DOUBLET,
            },
        }),
    )]
    #[case::above_maximum(
        UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Lit(2) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::Incompatible {
                unpaired_electrons: 0,
                multiplicity: SpinMultiplicity::DOUBLET,
            },
        }),
    )]
    fn test_spin_invariants_validator_validate_atom(
        #[case] unpaired_electrons: UnpairedElectronsForm,
        #[case] expected: Solution<(), SpinInvariantsContradiction>,
    ) {
        assert_eq!(
            SpinInvariantsValidator
                .validate_atom(&AtomForm {
                    unpaired_electrons,
                    ..Default::default()
                })
                .unwrap(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::all_entity_pairs_valid(
        MoleculeEntries {
            atoms: vec![
                AtomForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() },
                AtomForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() },
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() })],
            aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() })],
            multicenter: vec![(vec![AtomId(0), AtomId(1)], MulticenterBondForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() })],
            ..Default::default()
        },
        Solution::Determined(()),
    )]
    #[case::partial_pair(
        MoleculeEntries {
            atoms: vec![AtomForm {
                unpaired_electrons: UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(1) },
                ..Default::default()
            }],
            ..Default::default()
        },
        Solution::Underdetermined(()),
    )]
    #[case::molecule_atom_reports_id(
        MoleculeEntries {
            atoms: vec![
                AtomForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() },
                AtomForm { unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)), ..Default::default() },
            ],
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::MoleculeAtom {
            atom: AtomId(1),
            error: SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::DOUBLET },
        }),
    )]
    #[case::bond_reports_id(
        MoleculeEntries {
            atoms: vec![
                AtomForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() },
                AtomForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() },
                AtomForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() },
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() }),
                (AtomId(1), AtomId(2), BondForm { unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)), ..Default::default() }),
            ],
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::Bond {
            bond: BondId(1),
            error: SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::DOUBLET },
        }),
    )]
    #[case::aromatic_system_reports_id(
        MoleculeEntries {
            atoms: vec![AtomForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() }; 4],
            aromatic: vec![
                (vec![AtomId(0), AtomId(1)], AromaticSystemForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() }),
                (vec![AtomId(2), AtomId(3)], AromaticSystemForm { unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)), ..Default::default() }),
            ],
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::AromaticSystem {
            system: AromaticSystemId(1),
            error: SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::DOUBLET },
        }),
    )]
    #[case::multicenter_bond_reports_id(
        MoleculeEntries {
            atoms: vec![AtomForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() }; 6],
            multicenter: vec![
                (vec![AtomId(0), AtomId(1), AtomId(2)], MulticenterBondForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() }),
                (vec![AtomId(3), AtomId(4), AtomId(5)], MulticenterBondForm { unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)), ..Default::default() }),
            ],
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::MulticenterBond {
            bond: MulticenterBondId(1),
            error: SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::DOUBLET },
        }),
    )]
    #[case::underdetermined_does_not_mask_later_contradiction(
        MoleculeEntries {
            atoms: vec![
                AtomForm {
                    unpaired_electrons: UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(1) },
                    ..Default::default()
                },
                AtomForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() },
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() })],
            aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() })],
            multicenter: vec![(vec![AtomId(0), AtomId(1)], MulticenterBondForm { unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)), ..Default::default() })],
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::MulticenterBond {
            bond: MulticenterBondId(0),
            error: SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::DOUBLET },
        }),
    )]
    fn test_spin_invariants_validator_validate(
        #[case] entries: MoleculeEntries,
        #[case] expected: Solution<(), SpinInvariantsContradiction>,
    ) {
        assert_eq!(
            SpinInvariantsValidator
                .validate(&Molecule::from_entries(entries))
                .unwrap(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::unrelated_constraint(
        MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }).into(),
            ..Default::default()
        },
        Solution::Determined(()),
    )]
    #[case::vacuous_coupling(
        MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm::default(),
            }).into(),
            ..Default::default()
        },
        Solution::Determined(()),
    )]
    #[case::partial_coupling(
        MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm {
                    count: NumForm::Lit(2),
                    multiplicity: NumForm::Undetermined,
                },
            }).into(),
            ..Default::default()
        },
        Solution::Underdetermined(()),
    )]
    #[case::valid_coupling_not_yet_evaluated(
        MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 3_u8)),
            }).into(),
            ..Default::default()
        },
        Solution::Underdetermined(()),
    )]
    #[case::invalid_coupling(
        MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
            }).into(),
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::UnpairedElectronCoupling {
            constraint_index: 0,
            error: SpinStateError::Incompatible {
                unpaired_electrons: 2,
                multiplicity: SpinMultiplicity::DOUBLET,
            },
        }),
    )]
    #[case::invalid_coupling_nested_in_and_or_not(
        MoleculeEntries {
            constraints: vec![
                Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
                Constraint::Or(vec![
                    Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
                    Constraint::And(vec![Constraint::Not(Box::new(Constraint::Molecule(
                        MoleculeConstraint::UnpairedElectronCoupling {
                            atoms: None,
                            unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
                        },
                    )))]),
                ]),
            ].into(),
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::UnpairedElectronCoupling {
            constraint_index: 1,
            error: SpinStateError::Incompatible {
                unpaired_electrons: 2,
                multiplicity: SpinMultiplicity::DOUBLET,
            },
        }),
    )]
    #[case::valid_coupling_nested_in_not(
        MoleculeEntries {
            constraints: Constraint::Not(Box::new(Constraint::Molecule(
                MoleculeConstraint::UnpairedElectronCoupling {
                    atoms: None,
                    unpaired_electrons: UnpairedElectronsForm::from((2_u8, 3_u8)),
                },
            ))).into(),
            ..Default::default()
        },
        Solution::Underdetermined(()),
    )]
    #[case::earlier_underdetermination_does_not_mask_invalid_coupling(
        MoleculeEntries {
            constraints: vec![
                Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                    atoms: None,
                    unpaired_electrons: UnpairedElectronsForm::from((2_u8, 3_u8)),
                }),
                Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
                Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                    atoms: None,
                    unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
                }),
            ].into(),
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::UnpairedElectronCoupling {
            constraint_index: 2,
            error: SpinStateError::Incompatible {
                unpaired_electrons: 2,
                multiplicity: SpinMultiplicity::DOUBLET,
            },
        }),
    )]
    #[case::underdetermined_entity_does_not_mask_invalid_coupling(
        MoleculeEntries {
            atoms: vec![AtomForm {
                unpaired_electrons: UnpairedElectronsForm {
                    count: NumForm::Undetermined,
                    multiplicity: NumForm::Lit(1),
                },
                ..Default::default()
            }],
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
            }).into(),
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::UnpairedElectronCoupling {
            constraint_index: 0,
            error: SpinStateError::Incompatible {
                unpaired_electrons: 2,
                multiplicity: SpinMultiplicity::DOUBLET,
            },
        }),
    )]
    fn test_spin_invariants_validator_validate_constraints(
        #[case] entries: MoleculeEntries,
        #[case] expected: Solution<(), SpinInvariantsContradiction>,
    ) {
        assert_eq!(
            SpinInvariantsValidator
                .validate(&Molecule::from_entries(entries))
                .unwrap(),
            expected,
        );
    }

    proptest! {
        #[test]
        fn test_spin_invariants_validator_entity_aggregation(
            atom_state in 0_u8..3,
            bond_state in 0_u8..3,
            aromatic_state in 0_u8..3,
            multicenter_state in 0_u8..3,
        ) {
            let state_pair = |state| match state {
                0 => UnpairedElectronsForm::closed_shell(),
                1 => UnpairedElectronsForm {
                    count: NumForm::Undetermined,
                    multiplicity: NumForm::Lit(1),
                },
                2 => UnpairedElectronsForm::from((2_u8, 2_u8)),
                _ => unreachable!("strategy only generates states 0..3"),
            };
            let molecule = Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm { unpaired_electrons: state_pair(atom_state), ..Default::default() },
                    AtomForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() },
                ],
                bonds: vec![(
                    AtomId(0),
                    AtomId(1),
                    BondForm { unpaired_electrons: state_pair(bond_state), ..Default::default() },
                )],
                aromatic: vec![(
                    vec![AtomId(0), AtomId(1)],
                    AromaticSystemForm { unpaired_electrons: state_pair(aromatic_state), ..Default::default() },
                )],
                multicenter: vec![(
                    vec![AtomId(0), AtomId(1)],
                    MulticenterBondForm { unpaired_electrons: state_pair(multicenter_state), ..Default::default() },
                )],
                ..Default::default()
            });
            let error = SpinStateError::Incompatible {
                unpaired_electrons: 2,
                multiplicity: SpinMultiplicity::DOUBLET,
            };
            let expected = if atom_state == 2 {
                Solution::Contradictory(SpinInvariantsContradiction::MoleculeAtom {
                    atom: AtomId(0),
                    error,
                })
            } else if bond_state == 2 {
                Solution::Contradictory(SpinInvariantsContradiction::Bond {
                    bond: BondId(0),
                    error,
                })
            } else if aromatic_state == 2 {
                Solution::Contradictory(SpinInvariantsContradiction::AromaticSystem {
                    system: AromaticSystemId(0),
                    error,
                })
            } else if multicenter_state == 2 {
                Solution::Contradictory(SpinInvariantsContradiction::MulticenterBond {
                    bond: MulticenterBondId(0),
                    error,
                })
            } else if [atom_state, bond_state, aromatic_state, multicenter_state].contains(&1) {
                Solution::Underdetermined(())
            } else {
                Solution::Determined(())
            };

            prop_assert_eq!(SpinInvariantsValidator.validate(&molecule).unwrap(), expected);
        }
    }
}
