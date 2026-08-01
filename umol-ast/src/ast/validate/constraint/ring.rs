//! Constraints evaluated against the fixed Relevant ring projection through size 22.

use thiserror::Error;
use umol_graph_core::RelevantCycleEnumerationAlgorithm;
use umol_utils::solution::Solution;

use super::super::super::constraint::{
    AtomConstraintAst, BondConstraintAst, DativeBondConstraintAst,
};
use super::super::super::entity::Entity;
use super::super::super::id::{AtomId, BondId, DativeBondId};
use super::super::super::molecule::MoleculeAst;
use super::super::super::ring::{RingConfig, RingModel};
use super::super::super::traits::Lattice;
use super::super::super::value::ValueAst;
use super::super::super::view::{RingAtomView, RingBondView};
use super::ConstraintError;

/// Evaluates ring constraints with an explicit relevant-cycle algorithm selector.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RingConstraintValidator;

impl RingConstraintValidator {
    /// Validate every inline ring constraint against the fixed Relevant ring projection.
    pub fn validate(
        &self,
        ast: &MoleculeAst,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> Result<Solution<(), RingConstraintContradiction>, ConstraintError> {
        if !uses_ring_constraints(ast) {
            return Ok(Solution::Determined(()));
        }
        let rings = ast.rings(
            RingModel::default(),
            RingConfig {
                relevant_cycle_algorithm,
                ..RingConfig::default()
            },
        );
        let mut any_underdetermined = false;

        for id in ast.atoms().ids() {
            if let Some(contradiction) = observe(
                rings
                    .atom(id)
                    .validate_constraints(id, ast.atom(id).constraints()),
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }

        for id in ast.bonds().ids() {
            if let Some(contradiction) = observe(
                rings
                    .bond(id)
                    .validate_constraints(id, ast.bond(id).constraints()),
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }

        for bond in ast.dative_bonds().iter() {
            if bond.constraints().iter().any(|constraint| {
                matches!(
                    constraint,
                    DativeBondConstraintAst::RingMembership(membership)
                        if !membership.is_undetermined()
                )
            }) {
                return Err(ConstraintError::DativeBondRingMembershipUnsupported { bond: bond.id });
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
        ast: &MoleculeAst,
        atom_id: AtomId,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> Result<Solution<(), RingConstraintContradiction>, ConstraintError> {
        let atom = ast
            .atoms()
            .get(atom_id)
            .ok_or(ConstraintError::InvalidReference {
                entity: Entity::Atom(atom_id),
            })?;
        if !atom.constraints().iter().any(is_atom_ring_constraint) {
            return Ok(Solution::Determined(()));
        }
        let rings = ast.rings(
            RingModel::default(),
            RingConfig {
                relevant_cycle_algorithm,
                ..RingConfig::default()
            },
        );
        Ok(rings
            .atom(atom_id)
            .validate_constraints(atom_id, atom.constraints()))
    }

    /// Validate all inline ring constraints on one molecule bond.
    pub fn validate_molecule_bond(
        &self,
        ast: &MoleculeAst,
        bond_id: BondId,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> Result<Solution<(), RingConstraintContradiction>, ConstraintError> {
        let bond = ast
            .bonds()
            .get(bond_id)
            .ok_or(ConstraintError::InvalidReference {
                entity: Entity::Bond(bond_id),
            })?;
        if !bond.constraints().iter().any(is_bond_ring_constraint) {
            return Ok(Solution::Determined(()));
        }
        let rings = ast.rings(
            RingModel::default(),
            RingConfig {
                relevant_cycle_algorithm,
                ..RingConfig::default()
            },
        );
        Ok(rings
            .bond(bond_id)
            .validate_constraints(bond_id, bond.constraints()))
    }

    /// Validate inline ring constraints on one molecule dative bond.
    pub fn validate_molecule_dative_bond(
        &self,
        ast: &MoleculeAst,
        bond_id: DativeBondId,
    ) -> Result<Solution<(), RingConstraintContradiction>, ConstraintError> {
        let bond = ast
            .dative_bonds()
            .get(bond_id)
            .ok_or(ConstraintError::InvalidReference {
                entity: Entity::DativeBond(bond_id),
            })?;
        if bond.constraints().iter().any(|constraint| {
            matches!(
                constraint,
                DativeBondConstraintAst::RingMembership(membership)
                    if !membership.is_undetermined()
            )
        }) {
            Err(ConstraintError::DativeBondRingMembershipUnsupported { bond: bond_id })
        } else {
            Ok(Solution::Determined(()))
        }
    }
}

impl RingAtomView<'_> {
    fn validate_constraints(
        &self,
        atom_id: AtomId,
        constraints: &super::super::super::constraint::AtomConstraintsAst,
    ) -> Solution<(), RingConstraintContradiction> {
        conjunction(constraints.iter().filter_map(|constraint| {
            let (asserted, derived) = match constraint {
                AtomConstraintAst::RingDegree(asserted) => (asserted, self.ring_degree()),
                AtomConstraintAst::RingValence(asserted) => (asserted, self.ring_valence()),
                AtomConstraintAst::RingMembership(membership) => {
                    (&membership.count, self.ring_membership(membership.scope))
                }
                _ => return None,
            };
            Some(evaluate(
                asserted,
                &derived,
                RingConstraintContradiction::Atom {
                    atom: atom_id,
                    constraint: constraint.clone(),
                },
            ))
        }))
    }
}

impl RingBondView<'_> {
    fn validate_constraints(
        &self,
        bond_id: BondId,
        constraints: &super::super::super::constraint::BondConstraintsAst,
    ) -> Solution<(), RingConstraintContradiction> {
        conjunction(constraints.iter().filter_map(|constraint| {
            let BondConstraintAst::RingMembership(membership) = constraint else {
                return None;
            };
            Some(evaluate(
                &membership.count,
                &self.ring_membership(membership.scope),
                RingConstraintContradiction::Bond {
                    bond: bond_id,
                    constraint: constraint.clone(),
                },
            ))
        }))
    }
}

#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum RingConstraintContradiction {
    #[error("atom {atom:?} does not satisfy ring constraint {constraint:?}")]
    Atom {
        atom: AtomId,
        constraint: AtomConstraintAst,
    },
    #[error("bond {bond:?} does not satisfy ring constraint {constraint:?}")]
    Bond {
        bond: BondId,
        constraint: BondConstraintAst,
    },
}

fn evaluate(
    asserted: &ValueAst,
    derived: &ValueAst,
    contradiction: RingConstraintContradiction,
) -> Solution<(), RingConstraintContradiction> {
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
    outcome: Solution<(), RingConstraintContradiction>,
    any_underdetermined: &mut bool,
) -> Option<RingConstraintContradiction> {
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
    outcomes: impl IntoIterator<Item = Solution<(), RingConstraintContradiction>>,
) -> Solution<(), RingConstraintContradiction> {
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

fn uses_ring_constraints(ast: &MoleculeAst) -> bool {
    ast.atoms()
        .iter()
        .any(|atom| atom.constraints().iter().any(is_atom_ring_constraint))
        || ast
            .bonds()
            .iter()
            .any(|bond| bond.constraints().iter().any(is_bond_ring_constraint))
        || ast.dative_bonds().iter().any(|bond| {
            bond.constraints().iter().any(|constraint| {
                matches!(
                    constraint,
                    DativeBondConstraintAst::RingMembership(membership)
                        if !membership.is_undetermined()
                )
            })
        })
}

fn is_atom_ring_constraint(constraint: &AtomConstraintAst) -> bool {
    matches!(
        constraint,
        AtomConstraintAst::RingDegree(value) | AtomConstraintAst::RingValence(value)
            if !value.is_undetermined()
    ) || matches!(
        constraint,
        AtomConstraintAst::RingMembership(membership) if !membership.is_undetermined()
    )
}

fn is_bond_ring_constraint(constraint: &BondConstraintAst) -> bool {
    matches!(
        constraint,
        BondConstraintAst::RingMembership(membership) if !membership.is_undetermined()
    )
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::ast::constraint::RingScope;
    use crate::ast::id::DativeBondId;
    use crate::mol_dsl;

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
        r#"{:atoms ["C#x*#y*#R*" "N"] :bonds [[0 1 "1#R*"]] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R*"}]}"#,
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
        Ok(Solution::Contradictory(RingConstraintContradiction::Atom {
            atom: AtomId(0),
            constraint: AtomConstraintAst::ring_degree(1),
        })),
    )]
    #[case::atom_ring_valence_contradiction(
        r#"{:atoms ["C#y3" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#,
        Ok(Solution::Contradictory(RingConstraintContradiction::Atom {
            atom: AtomId(0),
            constraint: AtomConstraintAst::ring_valence(3),
        })),
    )]
    #[case::atom_ring_membership_contradiction(
        r#"{:atoms ["C#R(6)" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#,
        Ok(Solution::Contradictory(RingConstraintContradiction::Atom {
            atom: AtomId(0),
            constraint: AtomConstraintAst::ring_membership(RingScope::Size(6), 1),
        })),
    )]
    #[case::bond_ring_membership_contradiction(
        r#"{:atoms ["C" "C"] :bonds [[0 1 "1#R"]]}"#,
        Ok(Solution::Contradictory(RingConstraintContradiction::Bond {
            bond: BondId(0),
            constraint: BondConstraintAst::ring_membership(RingScope::All, 1),
        })),
    )]
    #[case::dative_ring_membership_error(
        r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}]}"#,
        Err(ConstraintError::DativeBondRingMembershipUnsupported {
            bond: DativeBondId(0),
        }),
    )]
    fn test_ring_constraint_validator_validate(
        #[case] input: &str,
        #[case] expected: Result<Solution<(), RingConstraintContradiction>, ConstraintError>,
    ) {
        assert_eq!(
            RingConstraintValidator
                .validate(&mol_dsl!(input), RelevantCycleEnumerationAlgorithm::Vismara,),
            expected,
        );
    }
}
