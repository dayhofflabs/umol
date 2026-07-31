//! Spin-state invariant validation. A complete literal unpaired-electron count
//! and multiplicity must form a physically valid [`SpinState`]; a pair with
//! either component still non-literal is underdetermined.

use thiserror::Error;
use umol_ast::ast::{
    AromaticSystemId, AsLit, AtomAst, AtomId, BondId, MoleculeAst, MulticenterBondId,
    UnpairedElectronsAst,
};
use umol_chem::error::SpinStateError;
use umol_chem::spin::SpinState;
use umol_utils::solution::Solution;

#[derive(Clone, Copy, Debug, Default)]
pub struct SpinInvariantsValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SpinInvariantsContradiction {
    #[error("atom has invalid unpaired electrons: {error}")]
    Atom { error: SpinStateError },
    #[error("molecule atom {id} has invalid unpaired electrons: {error}")]
    MoleculeAtom { id: AtomId, error: SpinStateError },
    #[error("bond {id} has invalid unpaired electrons: {error}")]
    Bond { id: BondId, error: SpinStateError },
    #[error("aromatic system {id} has invalid unpaired electrons: {error}")]
    AromaticSystem {
        id: AromaticSystemId,
        error: SpinStateError,
    },
    #[error("multicenter bond {id} has invalid unpaired electrons: {error}")]
    MulticenterBond {
        id: MulticenterBondId,
        error: SpinStateError,
    },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SpinInvariantsError {}

impl SpinInvariantsValidator {
    pub fn validate(
        &self,
        ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), SpinInvariantsContradiction>, SpinInvariantsError> {
        let ast = ast.as_ref();
        let mut any_undetermined = false;

        for atom in ast.atoms().iter() {
            match validate_unpaired_electrons(atom.unpaired_electrons()) {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => any_undetermined = true,
                Solution::Contradictory(error) => {
                    return Ok(Solution::Contradictory(
                        SpinInvariantsContradiction::MoleculeAtom { id: atom.id, error },
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
                        id: bond.id,
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
                            id: system.id,
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
                        SpinInvariantsContradiction::MulticenterBond { id: bond.id, error },
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
        atom: &AtomAst,
    ) -> Result<Solution<(), SpinInvariantsContradiction>, SpinInvariantsError> {
        Ok(validate_unpaired_electrons(&atom.unpaired_electrons)
            .map_contradiction(|error| SpinInvariantsContradiction::Atom { error }))
    }
}

fn validate_unpaired_electrons(
    unpaired_electrons: &UnpairedElectronsAst,
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
    use umol_ast::ast::{AromaticSystemAst, BondAst, MoleculeParts, MulticenterBondAst, ValueAst};
    use umol_chem::spin::SpinMultiplicity;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::closed_shell((0_u8, 1_u8).into(), Solution::Determined(()))]
    #[case::doublet((1_u8, 2_u8).into(), Solution::Determined(()))]
    #[case::open_shell_singlet((2_u8, 1_u8).into(), Solution::Determined(()))]
    #[case::triplet((2_u8, 3_u8).into(), Solution::Determined(()))]
    #[case::count_undetermined(
        UnpairedElectronsAst { count: ValueAst::Undetermined, multiplicity: ValueAst::Lit(1) },
        Solution::Underdetermined(()),
    )]
    #[case::multiplicity_undetermined(
        UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined },
        Solution::Underdetermined(()),
    )]
    #[case::count_negative(
        UnpairedElectronsAst { count: ValueAst::Lit(-1), multiplicity: ValueAst::Lit(1) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::UnpairedElectronsOutOfRange { count: -1 },
        }),
    )]
    #[case::count_above_u8(
        UnpairedElectronsAst { count: ValueAst::Lit(256), multiplicity: ValueAst::Lit(1) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::UnpairedElectronsOutOfRange { count: 256 },
        }),
    )]
    #[case::multiplicity_zero(
        UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Lit(0) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::MultiplicityOutOfRange { multiplicity: 0 },
        }),
    )]
    #[case::multiplicity_above_u8(
        UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Lit(256) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::MultiplicityOutOfRange { multiplicity: 256 },
        }),
    )]
    #[case::parity(
        UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(2) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::Incompatible {
                unpaired_electrons: 2,
                multiplicity: SpinMultiplicity::DOUBLET,
            },
        }),
    )]
    #[case::above_maximum(
        UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Lit(2) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::Incompatible {
                unpaired_electrons: 0,
                multiplicity: SpinMultiplicity::DOUBLET,
            },
        }),
    )]
    fn test_spin_invariants_validator_validate_atom(
        #[case] unpaired_electrons: UnpairedElectronsAst,
        #[case] expected: Solution<(), SpinInvariantsContradiction>,
    ) {
        assert_eq!(
            SpinInvariantsValidator
                .validate_atom(&AtomAst {
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
        MoleculeParts {
            atoms: vec![
                AtomAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() },
                AtomAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() },
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() })],
            aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() })],
            multicenter: vec![(vec![AtomId(0), AtomId(1)], MulticenterBondAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() })],
            ..Default::default()
        },
        Solution::Determined(()),
    )]
    #[case::partial_pair(
        MoleculeParts {
            atoms: vec![AtomAst {
                unpaired_electrons: UnpairedElectronsAst { count: ValueAst::Undetermined, multiplicity: ValueAst::Lit(1) },
                ..Default::default()
            }],
            ..Default::default()
        },
        Solution::Underdetermined(()),
    )]
    #[case::molecule_atom_reports_id(
        MoleculeParts {
            atoms: vec![
                AtomAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() },
                AtomAst { unpaired_electrons: UnpairedElectronsAst::from((2_u8, 2_u8)), ..Default::default() },
            ],
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::MoleculeAtom {
            id: AtomId(1),
            error: SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::DOUBLET },
        }),
    )]
    #[case::bond_reports_id(
        MoleculeParts {
            atoms: vec![
                AtomAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() },
                AtomAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() },
                AtomAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() },
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() }),
                (AtomId(1), AtomId(2), BondAst { unpaired_electrons: UnpairedElectronsAst::from((2_u8, 2_u8)), ..Default::default() }),
            ],
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::Bond {
            id: BondId(1),
            error: SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::DOUBLET },
        }),
    )]
    #[case::aromatic_system_reports_id(
        MoleculeParts {
            atoms: vec![AtomAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() }; 4],
            aromatic: vec![
                (vec![AtomId(0), AtomId(1)], AromaticSystemAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() }),
                (vec![AtomId(2), AtomId(3)], AromaticSystemAst { unpaired_electrons: UnpairedElectronsAst::from((2_u8, 2_u8)), ..Default::default() }),
            ],
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::AromaticSystem {
            id: AromaticSystemId(1),
            error: SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::DOUBLET },
        }),
    )]
    #[case::multicenter_bond_reports_id(
        MoleculeParts {
            atoms: vec![AtomAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() }; 6],
            multicenter: vec![
                (vec![AtomId(0), AtomId(1), AtomId(2)], MulticenterBondAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() }),
                (vec![AtomId(3), AtomId(4), AtomId(5)], MulticenterBondAst { unpaired_electrons: UnpairedElectronsAst::from((2_u8, 2_u8)), ..Default::default() }),
            ],
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::MulticenterBond {
            id: MulticenterBondId(1),
            error: SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::DOUBLET },
        }),
    )]
    #[case::underdetermined_does_not_mask_later_contradiction(
        MoleculeParts {
            atoms: vec![
                AtomAst {
                    unpaired_electrons: UnpairedElectronsAst { count: ValueAst::Undetermined, multiplicity: ValueAst::Lit(1) },
                    ..Default::default()
                },
                AtomAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() },
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() })],
            aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() })],
            multicenter: vec![(vec![AtomId(0), AtomId(1)], MulticenterBondAst { unpaired_electrons: UnpairedElectronsAst::from((2_u8, 2_u8)), ..Default::default() })],
            ..Default::default()
        },
        Solution::Contradictory(SpinInvariantsContradiction::MulticenterBond {
            id: MulticenterBondId(0),
            error: SpinStateError::Incompatible { unpaired_electrons: 2, multiplicity: SpinMultiplicity::DOUBLET },
        }),
    )]
    fn test_spin_invariants_validator_validate(
        #[case] parts: MoleculeParts,
        #[case] expected: Solution<(), SpinInvariantsContradiction>,
    ) {
        assert_eq!(
            SpinInvariantsValidator
                .validate(MoleculeAst::from_parts(parts))
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
                0 => UnpairedElectronsAst::closed_shell(),
                1 => UnpairedElectronsAst {
                    count: ValueAst::Undetermined,
                    multiplicity: ValueAst::Lit(1),
                },
                2 => UnpairedElectronsAst::from((2_u8, 2_u8)),
                _ => unreachable!("strategy only generates states 0..3"),
            };
            let molecule = MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![
                    AtomAst { unpaired_electrons: state_pair(atom_state), ..Default::default() },
                    AtomAst { unpaired_electrons: UnpairedElectronsAst::closed_shell(), ..Default::default() },
                ],
                bonds: vec![(
                    AtomId(0),
                    AtomId(1),
                    BondAst { unpaired_electrons: state_pair(bond_state), ..Default::default() },
                )],
                aromatic: vec![(
                    vec![AtomId(0), AtomId(1)],
                    AromaticSystemAst { unpaired_electrons: state_pair(aromatic_state), ..Default::default() },
                )],
                multicenter: vec![(
                    vec![AtomId(0), AtomId(1)],
                    MulticenterBondAst { unpaired_electrons: state_pair(multicenter_state), ..Default::default() },
                )],
                ..Default::default()
            });
            let error = SpinStateError::Incompatible {
                unpaired_electrons: 2,
                multiplicity: SpinMultiplicity::DOUBLET,
            };
            let expected = if atom_state == 2 {
                Solution::Contradictory(SpinInvariantsContradiction::MoleculeAtom {
                    id: AtomId(0),
                    error,
                })
            } else if bond_state == 2 {
                Solution::Contradictory(SpinInvariantsContradiction::Bond {
                    id: BondId(0),
                    error,
                })
            } else if aromatic_state == 2 {
                Solution::Contradictory(SpinInvariantsContradiction::AromaticSystem {
                    id: AromaticSystemId(0),
                    error,
                })
            } else if multicenter_state == 2 {
                Solution::Contradictory(SpinInvariantsContradiction::MulticenterBond {
                    id: MulticenterBondId(0),
                    error,
                })
            } else if [atom_state, bond_state, aromatic_state, multicenter_state].contains(&1) {
                Solution::Underdetermined(())
            } else {
                Solution::Determined(())
            };

            prop_assert_eq!(SpinInvariantsValidator.validate(molecule).unwrap(), expected);
        }
    }
}
