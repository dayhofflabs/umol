//! Molecule-scope aggregate and connectivity constraint evaluation.

use std::collections::BTreeSet;

use thiserror::Error;
use umol_graph_core::ConnectedComponentsAlgorithm;
use umol_utils::solution::Solution;

use super::super::super::constraint::MoleculeConstraint;
use super::super::super::entity::Entity;
use super::super::super::id::{AtomId, BondId};
use super::super::super::molecule::MoleculeAst;
use super::super::super::traits::Lattice;
use super::super::super::value::ValueAst;
use super::ConstraintError;

/// Evaluates one molecule-scope aggregate or connectivity constraint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MoleculeConstraintValidator;

#[derive(Clone, Debug, Error, PartialEq, Eq)]
#[error("molecule constraint is not satisfied: {constraint:?}")]
pub struct MoleculeConstraintContradiction {
    pub constraint: MoleculeConstraint,
}

impl MoleculeConstraintValidator {
    pub fn validate(
        &self,
        ast: &MoleculeAst,
        constraint: &MoleculeConstraint,
        connected_components_algorithm: ConnectedComponentsAlgorithm,
    ) -> Result<Solution<(), MoleculeConstraintContradiction>, ConstraintError> {
        let determined = match constraint {
            MoleculeConstraint::ChargeSum { atoms, sum } => {
                let atoms = atom_subset(ast, atoms.as_deref())?;
                let derived = atoms
                    .into_iter()
                    .map(|atom| ast.atom(atom).charge())
                    .fold(ValueAst::Lit(0), |sum, charge| sum + charge);
                return Ok(evaluate(sum, &derived, constraint));
            }
            MoleculeConstraint::UnpairedElectronCoupling {
                atoms,
                unpaired_electrons,
            } => {
                atom_subset(ast, atoms.as_deref())?;
                if unpaired_electrons.is_undetermined() {
                    true
                } else {
                    return Ok(Solution::Underdetermined(()));
                }
            }
            MoleculeConstraint::BondOrderSum { bonds, sum } => {
                let bonds = bond_subset(ast, bonds.as_deref())?;
                let derived = bonds
                    .into_iter()
                    .map(|bond| ast.bond(bond).order())
                    .fold(ValueAst::Lit(0), |sum, order| sum + order);
                return Ok(evaluate(sum, &derived, constraint));
            }
            MoleculeConstraint::Connected { atoms } => {
                let atoms = atom_subset(ast, atoms.as_deref())?;
                connected(ast, &atoms, connected_components_algorithm)
            }
            MoleculeConstraint::SubPattern { .. } => return Ok(Solution::Underdetermined(())),
        };

        Ok(if determined {
            Solution::Determined(())
        } else {
            Solution::Contradictory(MoleculeConstraintContradiction {
                constraint: constraint.clone(),
            })
        })
    }
}

fn evaluate(
    asserted: &ValueAst,
    derived: &ValueAst,
    constraint: &MoleculeConstraint,
) -> Solution<(), MoleculeConstraintContradiction> {
    if asserted.is_undetermined() {
        Solution::Determined(())
    } else if !derived.is_ground() {
        Solution::Underdetermined(())
    } else if asserted.matches(derived) {
        Solution::Determined(())
    } else {
        Solution::Contradictory(MoleculeConstraintContradiction {
            constraint: constraint.clone(),
        })
    }
}

fn atom_subset(
    ast: &MoleculeAst,
    atoms: Option<&[AtomId]>,
) -> Result<Vec<AtomId>, ConstraintError> {
    match atoms {
        Some(atoms) => {
            let mut selected = BTreeSet::new();
            for &atom in atoms {
                if !ast.atoms().contains(atom) {
                    return Err(ConstraintError::InvalidReference {
                        entity: Entity::Atom(atom),
                    });
                }
                selected.insert(atom);
            }
            Ok(selected.into_iter().collect())
        }
        None => Ok(ast.atoms().ids().collect()),
    }
}

fn bond_subset(
    ast: &MoleculeAst,
    bonds: Option<&[BondId]>,
) -> Result<Vec<BondId>, ConstraintError> {
    match bonds {
        Some(bonds) => {
            let mut selected = BTreeSet::new();
            for &bond in bonds {
                if !ast.bonds().contains(bond) {
                    return Err(ConstraintError::InvalidReference {
                        entity: Entity::Bond(bond),
                    });
                }
                selected.insert(bond);
            }
            Ok(selected.into_iter().collect())
        }
        None => Ok(ast.bonds().ids().collect()),
    }
}

/// Whether every selected atom belongs to one localized-bond component. Paths may pass through
/// atoms outside the selected subset; empty and singleton subsets are connected.
fn connected(ast: &MoleculeAst, atoms: &[AtomId], algorithm: ConnectedComponentsAlgorithm) -> bool {
    if atoms.len() < 2 {
        return true;
    }
    let selected: BTreeSet<_> = atoms.iter().copied().collect();
    ast.graph()
        .connected_components(algorithm)
        .into_iter()
        .any(|component| {
            component
                .into_iter()
                .filter(|atom| selected.contains(atom))
                .count()
                == atoms.len()
        })
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};

    use super::*;
    use crate::ast::spin::UnpairedElectronsAst;
    use crate::mol_dsl;

    #[fixture]
    fn aggregate_molecule() -> MoleculeAst {
        mol_dsl!(
            r#"{:atoms ["C#c+" "N#c-" "O#c2" "F#c0" "Cl#c0"]
                :bonds [[0 1 "1"] [1 2 "2"] [3 4 "1"]]}"#
        )
    }

    #[rstest]
    #[case::charge_subset(MoleculeConstraint::ChargeSum {
        atoms: Some(vec![AtomId(0), AtomId(1)]),
        sum: ValueAst::Lit(0),
    })]
    #[case::charge_all(MoleculeConstraint::ChargeSum {
        atoms: None,
        sum: ValueAst::Lit(2),
    })]
    #[case::charge_empty(MoleculeConstraint::ChargeSum {
        atoms: Some(vec![]),
        sum: ValueAst::Lit(0),
    })]
    #[case::bond_subset(MoleculeConstraint::BondOrderSum {
        bonds: Some(vec![BondId(0), BondId(1)]),
        sum: ValueAst::Lit(3),
    })]
    #[case::bond_all(MoleculeConstraint::BondOrderSum {
        bonds: None,
        sum: ValueAst::Lit(4),
    })]
    #[case::bond_empty(MoleculeConstraint::BondOrderSum {
        bonds: Some(vec![]),
        sum: ValueAst::Lit(0),
    })]
    #[case::coupling_vacuous(MoleculeConstraint::UnpairedElectronCoupling {
        atoms: None,
        unpaired_electrons: UnpairedElectronsAst::default(),
    })]
    #[case::connected_subset(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomId(0), AtomId(2)]),
    })]
    #[case::connected_empty(MoleculeConstraint::Connected {
        atoms: Some(vec![]),
    })]
    #[case::connected_singleton(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomId(3)]),
    })]
    fn test_molecule_constraint_validator_validate(
        aggregate_molecule: MoleculeAst,
        #[case] constraint: MoleculeConstraint,
    ) {
        assert_eq!(
            MoleculeConstraintValidator.validate(
                &aggregate_molecule,
                &constraint,
                ConnectedComponentsAlgorithm::Bfs,
            ),
            Ok(Solution::Determined(())),
        );
    }

    #[rstest]
    #[case::charge(MoleculeConstraint::ChargeSum {
        atoms: None,
        sum: ValueAst::Lit(0),
    })]
    #[case::bond_order(MoleculeConstraint::BondOrderSum {
        bonds: None,
        sum: ValueAst::Lit(0),
    })]
    #[case::connected_all(MoleculeConstraint::Connected { atoms: None })]
    #[case::connected_subset(MoleculeConstraint::Connected {
        atoms: Some(vec![]),
    })]
    fn test_molecule_constraint_validator_validate_empty(#[case] constraint: MoleculeConstraint) {
        assert_eq!(
            MoleculeConstraintValidator.validate(
                &MoleculeAst::default(),
                &constraint,
                ConnectedComponentsAlgorithm::Bfs,
            ),
            Ok(Solution::Determined(())),
        );
    }

    #[rstest]
    #[case::charge(
        r#"{:atoms ["C#c+" "C"] :bonds []}"#,
        MoleculeConstraint::ChargeSum {
            atoms: None,
            sum: ValueAst::Lit(1),
        },
    )]
    #[case::bond_order(
        r#"{:atoms ["C" "C"] :bonds [[0 1 "*"]]}"#,
        MoleculeConstraint::BondOrderSum {
            bonds: None,
            sum: ValueAst::Lit(1),
        },
    )]
    #[case::coupling_literal(
        r#"{:atoms ["C#u0#s1"] :bonds []}"#,
        MoleculeConstraint::UnpairedElectronCoupling {
            atoms: None,
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
        },
    )]
    #[case::coupling_partial(
        r#"{:atoms ["C#u0#s1"] :bonds []}"#,
        MoleculeConstraint::UnpairedElectronCoupling {
            atoms: None,
            unpaired_electrons: UnpairedElectronsAst {
                count: ValueAst::Lit(0),
                multiplicity: ValueAst::Undetermined,
            },
        },
    )]
    fn test_molecule_constraint_validator_validate_partial(
        #[case] input: &str,
        #[case] constraint: MoleculeConstraint,
    ) {
        let molecule = mol_dsl!(input);

        assert_eq!(
            MoleculeConstraintValidator.validate(
                &molecule,
                &constraint,
                ConnectedComponentsAlgorithm::Bfs,
            ),
            Ok(Solution::Underdetermined(())),
        );
    }

    #[rstest]
    #[case::charge_subset(MoleculeConstraint::ChargeSum {
        atoms: Some(vec![AtomId(0), AtomId(1)]),
        sum: ValueAst::Lit(1),
    })]
    #[case::charge_all(MoleculeConstraint::ChargeSum {
        atoms: None,
        sum: ValueAst::Lit(0),
    })]
    #[case::bond_subset(MoleculeConstraint::BondOrderSum {
        bonds: Some(vec![BondId(0), BondId(1)]),
        sum: ValueAst::Lit(2),
    })]
    #[case::bond_all(MoleculeConstraint::BondOrderSum {
        bonds: None,
        sum: ValueAst::Lit(3),
    })]
    #[case::connected_subset(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomId(0), AtomId(3)]),
    })]
    #[case::connected_all(MoleculeConstraint::Connected { atoms: None })]
    fn test_molecule_constraint_validator_validate_contradiction(
        aggregate_molecule: MoleculeAst,
        #[case] constraint: MoleculeConstraint,
    ) {
        assert_eq!(
            MoleculeConstraintValidator.validate(
                &aggregate_molecule,
                &constraint,
                ConnectedComponentsAlgorithm::Bfs,
            ),
            Ok(Solution::Contradictory(MoleculeConstraintContradiction {
                constraint: constraint.clone(),
            },)),
        );
    }

    #[rstest]
    #[case::charge(
        MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(99)]),
            sum: ValueAst::Undetermined,
        },
        Entity::Atom(AtomId(99)),
    )]
    #[case::coupling(
        MoleculeConstraint::UnpairedElectronCoupling {
            atoms: Some(vec![AtomId(99)]),
            unpaired_electrons: UnpairedElectronsAst::default(),
        },
        Entity::Atom(AtomId(99)),
    )]
    #[case::bond_order(
        MoleculeConstraint::BondOrderSum {
            bonds: Some(vec![BondId(99)]),
            sum: ValueAst::Undetermined,
        },
        Entity::Bond(BondId(99)),
    )]
    #[case::connected(
        MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(99)]),
        },
        Entity::Atom(AtomId(99)),
    )]
    fn test_molecule_constraint_validator_validate_error(
        aggregate_molecule: MoleculeAst,
        #[case] constraint: MoleculeConstraint,
        #[case] entity: Entity,
    ) {
        assert_eq!(
            MoleculeConstraintValidator.validate(
                &aggregate_molecule,
                &constraint,
                ConnectedComponentsAlgorithm::Bfs,
            ),
            Err(ConstraintError::InvalidReference { entity }),
        );
    }
}
