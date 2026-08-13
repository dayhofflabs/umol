//! Constraints evaluated against the fixed Relevant ring projection through size 22.

use thiserror::Error;
use umol_graph_core::RelevantCycleEnumerationAlgorithm;
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomConstraintsForm, AtomId, BondConstraintForm, BondConstraintsForm,
    BondId, DativeBondConstraintForm, DativeBondId, Entity, Lattice, Molecule, NumForm,
    RingAtomView, RingBondView, RingConfig, RingModel,
};
use umol_utils::solution::Solution;

use super::ConstraintInvariantsError;

/// Evaluates ring constraints with an explicit relevant-cycle algorithm selector.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RingConstraintInvariantsValidator;

impl RingConstraintInvariantsValidator {
    /// Validate every inline ring constraint against the fixed Relevant ring projection.
    pub fn validate(
        &self,
        molecule: &Molecule,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> Result<Solution<(), RingConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        if !uses_ring_constraints(molecule) {
            return Ok(Solution::Determined(()));
        }
        let rings = molecule.rings(
            RingModel::default(),
            RingConfig {
                relevant_cycle_algorithm,
                ..RingConfig::default()
            },
        );
        let mut any_underdetermined = false;

        for id in molecule.atoms().ids() {
            if let Some(contradiction) = observe(
                validate_atom_constraints(
                    &rings.atom(id),
                    id,
                    &molecule.atom(id).attributes.constraints,
                ),
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }

        for id in molecule.bonds().ids() {
            if let Some(contradiction) = observe(
                validate_bond_constraints(
                    &rings.bond(id),
                    id,
                    &molecule.bond(id).attributes.constraints,
                ),
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }

        for bond in molecule.dative_bonds().iter() {
            if bond.constraints().iter().any(|constraint| {
                matches!(
                    constraint,
                    DativeBondConstraintForm::RingMembership(membership)
                        if !membership.is_undetermined()
                )
            }) {
                return Err(
                    ConstraintInvariantsError::DativeBondRingMembershipUnsupported {
                        bond: bond.id,
                    },
                );
            }
        }

        Ok(if any_underdetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        })
    }

    /// Validate all inline ring constraints on one molecule atom.
    pub fn validate_molecule_atom(
        &self,
        molecule: &Molecule,
        atom_id: AtomId,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> Result<Solution<(), RingConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let atom =
            molecule
                .atoms()
                .get(atom_id)
                .ok_or(ConstraintInvariantsError::InvalidReference {
                    entity: Entity::Atom(atom_id),
                })?;
        if !atom.constraints().iter().any(is_atom_ring_constraint) {
            return Ok(Solution::Determined(()));
        }
        let rings = molecule.rings(
            RingModel::default(),
            RingConfig {
                relevant_cycle_algorithm,
                ..RingConfig::default()
            },
        );
        Ok(validate_atom_constraints(
            &rings.atom(atom_id),
            atom_id,
            &atom.attributes.constraints,
        ))
    }

    /// Validate all inline ring constraints on one molecule bond.
    pub fn validate_molecule_bond(
        &self,
        molecule: &Molecule,
        bond_id: BondId,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> Result<Solution<(), RingConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let bond =
            molecule
                .bonds()
                .get(bond_id)
                .ok_or(ConstraintInvariantsError::InvalidReference {
                    entity: Entity::Bond(bond_id),
                })?;
        if !bond.constraints().iter().any(is_bond_ring_constraint) {
            return Ok(Solution::Determined(()));
        }
        let rings = molecule.rings(
            RingModel::default(),
            RingConfig {
                relevant_cycle_algorithm,
                ..RingConfig::default()
            },
        );
        Ok(validate_bond_constraints(
            &rings.bond(bond_id),
            bond_id,
            &bond.attributes.constraints,
        ))
    }

    /// Validate inline ring constraints on one molecule dative bond.
    pub fn validate_molecule_dative_bond(
        &self,
        molecule: &Molecule,
        bond_id: DativeBondId,
    ) -> Result<Solution<(), RingConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let bond = molecule.dative_bonds().get(bond_id).ok_or(
            ConstraintInvariantsError::InvalidReference {
                entity: Entity::DativeBond(bond_id),
            },
        )?;
        if bond.constraints().iter().any(|constraint| {
            matches!(
                constraint,
                DativeBondConstraintForm::RingMembership(membership)
                    if !membership.is_undetermined()
            )
        }) {
            Err(ConstraintInvariantsError::DativeBondRingMembershipUnsupported { bond: bond_id })
        } else {
            Ok(Solution::Determined(()))
        }
    }
}

fn validate_atom_constraints(
    view: &RingAtomView<'_>,
    atom_id: AtomId,
    constraints: &AtomConstraintsForm,
) -> Solution<(), RingConstraintInvariantsContradiction> {
    conjunction(constraints.iter().filter_map(|constraint| {
        let (asserted, derived) = match constraint {
            AtomConstraintForm::RingDegree(asserted) => (asserted, view.ring_degree()),
            AtomConstraintForm::RingValence(asserted) => (asserted, view.ring_valence()),
            AtomConstraintForm::RingMembership(membership) => {
                (&membership.count, view.ring_membership(membership.scope))
            }
            _ => return None,
        };
        Some(evaluate(
            asserted,
            &derived,
            RingConstraintInvariantsContradiction::Atom {
                atom: atom_id,
                constraint: constraint.clone(),
            },
        ))
    }))
}

fn validate_bond_constraints(
    view: &RingBondView<'_>,
    bond_id: BondId,
    constraints: &BondConstraintsForm,
) -> Solution<(), RingConstraintInvariantsContradiction> {
    conjunction(constraints.iter().filter_map(|constraint| {
        let BondConstraintForm::RingMembership(membership) = constraint else {
            return None;
        };
        Some(evaluate(
            &membership.count,
            &view.ring_membership(membership.scope),
            RingConstraintInvariantsContradiction::Bond {
                bond: bond_id,
                constraint: constraint.clone(),
            },
        ))
    }))
}

#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum RingConstraintInvariantsContradiction {
    #[error("atom {atom:?} does not satisfy ring constraint {constraint:?}")]
    Atom {
        atom: AtomId,
        constraint: AtomConstraintForm,
    },
    #[error("bond {bond:?} does not satisfy ring constraint {constraint:?}")]
    Bond {
        bond: BondId,
        constraint: BondConstraintForm,
    },
}

fn evaluate(
    asserted: &NumForm,
    derived: &NumForm,
    contradiction: RingConstraintInvariantsContradiction,
) -> Solution<(), RingConstraintInvariantsContradiction> {
    if asserted.is_undetermined() {
        Solution::Determined(())
    } else if !derived.is_ground() {
        Solution::Underdetermined(())
    } else if asserted.matches(derived) {
        Solution::Determined(())
    } else {
        Solution::Contradictory(contradiction)
    }
}

fn observe(
    outcome: Solution<(), RingConstraintInvariantsContradiction>,
    any_underdetermined: &mut bool,
) -> Option<RingConstraintInvariantsContradiction> {
    match outcome {
        Solution::Determined(()) => None,
        Solution::Underdetermined(()) => {
            *any_underdetermined = true;
            None
        }
        Solution::Contradictory(contradiction) => Some(contradiction),
    }
}

fn conjunction(
    outcomes: impl IntoIterator<Item = Solution<(), RingConstraintInvariantsContradiction>>,
) -> Solution<(), RingConstraintInvariantsContradiction> {
    let mut any_underdetermined = false;
    for outcome in outcomes {
        if let Some(contradiction) = observe(outcome, &mut any_underdetermined) {
            return Solution::Contradictory(contradiction);
        }
    }
    if any_underdetermined {
        Solution::Underdetermined(())
    } else {
        Solution::Determined(())
    }
}

fn uses_ring_constraints(molecule: &Molecule) -> bool {
    molecule
        .atoms()
        .iter()
        .any(|atom| atom.constraints().iter().any(is_atom_ring_constraint))
        || molecule
            .bonds()
            .iter()
            .any(|bond| bond.constraints().iter().any(is_bond_ring_constraint))
        || molecule.dative_bonds().iter().any(|bond| {
            bond.constraints().iter().any(|constraint| {
                matches!(
                    constraint,
                    DativeBondConstraintForm::RingMembership(membership)
                        if !membership.is_undetermined()
                )
            })
        })
}

fn is_atom_ring_constraint(constraint: &AtomConstraintForm) -> bool {
    matches!(
        constraint,
        AtomConstraintForm::RingDegree(value) | AtomConstraintForm::RingValence(value)
            if !value.is_undetermined()
    ) || matches!(
        constraint,
        AtomConstraintForm::RingMembership(membership) if !membership.is_undetermined()
    )
}

fn is_bond_ring_constraint(constraint: &BondConstraintForm) -> bool {
    matches!(
        constraint,
        BondConstraintForm::RingMembership(membership) if !membership.is_undetermined()
    )
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{DativeBondId, RingScope};
    use umol_graph_ir::mol_dsl;

    use super::*;

    #[rstest]
    #[case::acyclic(
        r#"{:atoms ["C#x0#y0#R!" "C"] :bonds [[0 1 "1#R!"]]}"#,
        Ok(Solution::Determined(())),
    )]
    #[case::monocyclic(
        r#"{:atoms ["C#x2#y2#R#R(3)" "C" "C"] :bonds [[0 1 "1#R#R(3)"] [1 2 "1"] [2 0 "1"]]}"#,
        Ok(Solution::Determined(())),
    )]
    #[case::fused(
        r#"{:atoms ["C" "C#x3#R2#R(3)2" "C" "C"] :bonds [[0 1 "1"] [1 2 "1#R2#R(3)2"] [2 0 "1"] [1 3 "1"] [3 2 "1"]]}"#,
        Ok(Solution::Determined(())),
    )]
    #[case::bridged(
        r#"{:atoms ["C#x3#R3#R(4)3" "C" "C" "C" "C"] :bonds [[0 2 "1#R2#R(4)2"] [2 1 "1"] [0 3 "1"] [3 1 "1"] [0 4 "1"] [4 1 "1"]]}"#,
        Ok(Solution::Determined(())),
    )]
    #[case::size_filtered(
        r#"{:atoms ["C#R(5)!#R(6)" "C" "C" "C" "C" "C"] :bonds [[0 1 "1#R(5)!#R(6)"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#,
        Ok(Solution::Determined(())),
    )]
    #[case::vacuous(
        r#"{:atoms ["C#x*#y*#R*" "N"] :bonds [[0 1 "1#R*"]] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1#R*"}]}"#,
        Ok(Solution::Determined(())),
    )]
    #[case::finite_set(
        r#"{:atoms ["C#R{1,2}" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#,
        Ok(Solution::Determined(())),
    )]
    #[case::partial_ring_valence(
        r#"{:atoms ["C#y2" "C" "C"] :bonds [[0 1 "*"] [1 2 "1"] [2 0 "1"]]}"#,
        Ok(Solution::Underdetermined(())),
    )]
    #[case::atom_ring_degree_contradiction(
        r#"{:atoms ["C#x1" "C"] :bonds [[0 1 "1"]]}"#,
        Ok(Solution::Contradictory(RingConstraintInvariantsContradiction::Atom {
            atom: AtomId(0),
            constraint: AtomConstraintForm::ring_degree(1),
        })),
    )]
    #[case::atom_ring_valence_contradiction(
        r#"{:atoms ["C#y3" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#,
        Ok(Solution::Contradictory(RingConstraintInvariantsContradiction::Atom {
            atom: AtomId(0),
            constraint: AtomConstraintForm::ring_valence(3),
        })),
    )]
    #[case::atom_ring_membership_contradiction(
        r#"{:atoms ["C#R(6)" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#,
        Ok(Solution::Contradictory(RingConstraintInvariantsContradiction::Atom {
            atom: AtomId(0),
            constraint: AtomConstraintForm::ring_membership(RingScope::Size(6), 1),
        })),
    )]
    #[case::bond_ring_membership_contradiction(
        r#"{:atoms ["C" "C"] :bonds [[0 1 "1#R"]]}"#,
        Ok(Solution::Contradictory(RingConstraintInvariantsContradiction::Bond {
            bond: BondId(0),
            constraint: BondConstraintForm::ring_membership(RingScope::All, 1),
        })),
    )]
    #[case::dative_ring_membership_error(
        r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1#R"}]}"#,
        Err(ConstraintInvariantsError::DativeBondRingMembershipUnsupported {
            bond: DativeBondId(0),
        }),
    )]
    fn test_ring_constraint_validator_validate(
        #[case] input: &str,
        #[case] expected: Result<
            Solution<(), RingConstraintInvariantsContradiction>,
            ConstraintInvariantsError,
        >,
    ) {
        assert_eq!(
            RingConstraintInvariantsValidator
                .validate(&mol_dsl!(input), RelevantCycleEnumerationAlgorithm::Vismara,),
            expected,
        );
    }
}
