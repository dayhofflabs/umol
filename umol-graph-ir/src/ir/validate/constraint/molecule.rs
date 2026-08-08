//! Molecule-scope aggregate and connectivity constraint evaluation.

use std::collections::BTreeSet;
use std::ops::ControlFlow;

use thiserror::Error;
use umol_graph_core::ConnectedComponentsAlgorithm;
use umol_utils::solution::Solution;

use super::super::super::constraint::{MoleculeConstraint, SubPatternAnchor};
use super::super::super::correspondence::MoleculeCorrespondence;
use super::super::super::entity::Entity;
use super::super::super::id::{AtomId, BondId};
use super::super::super::molecule::MoleculeAst;
use super::super::super::substructure::SubstructureMatchConfig;
use super::super::super::traits::Lattice;
use super::super::super::value::ValueAst;
use super::{ConstraintError, ConstraintValidateConfig};

/// Evaluates one molecule-scope aggregate, connectivity, or subpattern constraint.
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
        config: ConstraintValidateConfig,
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
                connected(ast, &atoms, config.connected_components_algorithm)
            }
            MoleculeConstraint::SubPattern { anchor, pattern } => {
                validate_anchor(ast, pattern, anchor)?;
                pattern
                    .visit_substructure_matches(
                        ast,
                        SubstructureMatchConfig {
                            match_algorithm: config.substructure_match_algorithm,
                            subgraph_isomorphism_algorithm: config.subgraph_isomorphism_algorithm,
                            relevant_cycle_algorithm: config.relevant_cycle_algorithm,
                        },
                        |correspondence| {
                            if anchor_matches(anchor, &correspondence) {
                                ControlFlow::Break(())
                            } else {
                                ControlFlow::Continue(())
                            }
                        },
                    )
                    .is_break()
            }
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

fn validate_anchor(
    host: &MoleculeAst,
    pattern: &MoleculeAst,
    anchor: &SubPatternAnchor,
) -> Result<(), ConstraintError> {
    for &(target, pattern_id) in anchor.atoms() {
        require_anchor_pair(
            host,
            pattern,
            Entity::Atom(target),
            Entity::Atom(pattern_id),
        )?;
    }
    for &(target, pattern_id) in anchor.bonds() {
        require_anchor_pair(
            host,
            pattern,
            Entity::Bond(target),
            Entity::Bond(pattern_id),
        )?;
    }
    for &(target, pattern_id) in anchor.dative_bonds() {
        require_anchor_pair(
            host,
            pattern,
            Entity::DativeBond(target),
            Entity::DativeBond(pattern_id),
        )?;
    }
    for &(target, pattern_id) in anchor.aromatic_systems() {
        require_anchor_pair(
            host,
            pattern,
            Entity::AromaticSystem(target),
            Entity::AromaticSystem(pattern_id),
        )?;
    }
    for &(target, pattern_id) in anchor.multicenter_bonds() {
        require_anchor_pair(
            host,
            pattern,
            Entity::MulticenterBond(target),
            Entity::MulticenterBond(pattern_id),
        )?;
    }
    for &(target, pattern_id) in anchor.noncovalent_bonds() {
        require_anchor_pair(
            host,
            pattern,
            Entity::NoncovalentBond(target),
            Entity::NoncovalentBond(pattern_id),
        )?;
    }
    for &(target, pattern_id) in anchor.stereo_atoms() {
        require_anchor_pair(
            host,
            pattern,
            Entity::StereoAtom(target),
            Entity::StereoAtom(pattern_id),
        )?;
    }
    for &(target, pattern_id) in anchor.stereo_bonds() {
        require_anchor_pair(
            host,
            pattern,
            Entity::StereoBond(target),
            Entity::StereoBond(pattern_id),
        )?;
    }
    Ok(())
}

fn require_anchor_pair(
    host: &MoleculeAst,
    pattern: &MoleculeAst,
    target: Entity,
    pattern_entity: Entity,
) -> Result<(), ConstraintError> {
    if !contains_entity(host, target) {
        return Err(ConstraintError::InvalidReference { entity: target });
    }
    if !contains_entity(pattern, pattern_entity) {
        return Err(ConstraintError::InvalidReference {
            entity: pattern_entity,
        });
    }
    Ok(())
}

fn contains_entity(molecule: &MoleculeAst, entity: Entity) -> bool {
    match entity {
        Entity::Atom(id) => molecule.atoms().contains(id),
        Entity::Bond(id) => molecule.bonds().contains(id),
        Entity::DativeBond(id) => molecule.dative_bonds().contains(id),
        Entity::AromaticSystem(id) => molecule.aromatic_systems().contains(id),
        Entity::MulticenterBond(id) => molecule.multicenter_bonds().contains(id),
        Entity::NoncovalentBond(id) => molecule.noncovalent_bonds().contains(id),
        Entity::StereoAtom(id) => molecule.stereo_atoms().contains(id),
        Entity::StereoBond(id) => molecule.stereo_bonds().contains(id),
    }
}

fn anchor_matches(anchor: &SubPatternAnchor, correspondence: &MoleculeCorrespondence) -> bool {
    anchor.atoms().iter().all(|&(target, pattern)| {
        correspondence.right_of(Entity::Atom(pattern)) == Some(Entity::Atom(target))
    }) && anchor.bonds().iter().all(|&(target, pattern)| {
        correspondence.right_of(Entity::Bond(pattern)) == Some(Entity::Bond(target))
    }) && anchor.dative_bonds().iter().all(|&(target, pattern)| {
        correspondence.right_of(Entity::DativeBond(pattern)) == Some(Entity::DativeBond(target))
    }) && anchor.aromatic_systems().iter().all(|&(target, pattern)| {
        correspondence.right_of(Entity::AromaticSystem(pattern))
            == Some(Entity::AromaticSystem(target))
    }) && anchor.multicenter_bonds().iter().all(|&(target, pattern)| {
        correspondence.right_of(Entity::MulticenterBond(pattern))
            == Some(Entity::MulticenterBond(target))
    }) && anchor.noncovalent_bonds().iter().all(|&(target, pattern)| {
        correspondence.right_of(Entity::NoncovalentBond(pattern))
            == Some(Entity::NoncovalentBond(target))
    }) && anchor.stereo_atoms().iter().all(|&(target, pattern)| {
        correspondence.right_of(Entity::StereoAtom(pattern)) == Some(Entity::StereoAtom(target))
    }) && anchor.stereo_bonds().iter().all(|&(target, pattern)| {
        correspondence.right_of(Entity::StereoBond(pattern)) == Some(Entity::StereoBond(target))
    })
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
        .enumerate_connected_components(algorithm)
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
    use umol_graph_core::SubgraphIsomorphismAlgorithm::{
        ArcMatch, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit,
    };
    use umol_graph_core::{
        RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
        ARCMATCH_DEFAULT_PATH_LENGTH,
    };

    use super::*;
    use crate::ir::id::{
        AromaticSystemId, DativeBondId, MulticenterBondId, NoncovalentBondId, StereoAtomId,
        StereoBondId,
    };
    use crate::ir::spin::UnpairedElectronsAst;
    use crate::ir::substructure::SubstructureMatchAlgorithm;
    use crate::mol_dsl;

    const CONFIG: ConstraintValidateConfig = ConstraintValidateConfig {
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
        connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
        substructure_match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
        subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit,
    };
    const SUBGRAPH_ISOMORPHISM_ALGORITHMS: [SubgraphIsomorphismAlgorithm; 6] = [
        Vf2,
        Ullmann,
        Ri,
        ArcMatch {
            path_length: ARCMATCH_DEFAULT_PATH_LENGTH,
        },
        Vf2Rdkit,
        RayKirsch,
    ];
    const SUBSTRUCTURE_MATCH_ALGORITHMS: [SubstructureMatchAlgorithm; 2] = [
        SubstructureMatchAlgorithm::GraphAndOverlays,
        SubstructureMatchAlgorithm::Incidence,
    ];

    #[fixture]
    fn aggregate_molecule() -> MoleculeAst {
        mol_dsl!(
            r#"{:atoms ["C#c+" "N#c-" "O#c2" "F#c0" "Cl#c0"]
                :bonds [[0 1 "1"] [1 2 "2"] [3 4 "1"]]}"#
        )
    }

    #[fixture]
    fn anchored_molecule() -> MoleculeAst {
        mol_dsl!(
            r#"{
                :atoms ["C" "F" "Cl" "Br" "I"]
                :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :dative-bonds [{:donors [1] :acceptor 0 :type "1"}]
                :aromatic-systems [{:atoms [0 1] :type "*#e2"}]
                :multicenter-bonds [{:atoms [0 1] :type "*#e2"}]
                :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]
                :stereo-bonds [{:site 0 :ligands [2 3] :type "Ct1"}]
            }"#
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
            MoleculeConstraintValidator.validate(&aggregate_molecule, &constraint, CONFIG,),
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
            MoleculeConstraintValidator.validate(&MoleculeAst::default(), &constraint, CONFIG,),
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
            MoleculeConstraintValidator.validate(&molecule, &constraint, CONFIG,),
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
            MoleculeConstraintValidator.validate(&aggregate_molecule, &constraint, CONFIG,),
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
            MoleculeConstraintValidator.validate(&aggregate_molecule, &constraint, CONFIG,),
            Err(ConstraintError::InvalidReference { entity }),
        );
    }

    #[rstest]
    fn test_molecule_constraint_validator_validate_subpattern_unanchored() {
        let host = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#);
        let constraint = MoleculeConstraint::SubPattern {
            anchor: SubPatternAnchor::new(),
            pattern: Box::new(mol_dsl!(r#"{:atoms ["O"] :bonds []}"#)),
        };

        for match_algorithm in SUBSTRUCTURE_MATCH_ALGORITHMS {
            for subgraph_isomorphism_algorithm in SUBGRAPH_ISOMORPHISM_ALGORITHMS {
                assert_eq!(
                    MoleculeConstraintValidator.validate(
                        &host,
                        &constraint,
                        ConstraintValidateConfig {
                            substructure_match_algorithm: match_algorithm,
                            subgraph_isomorphism_algorithm,
                            ..CONFIG
                        },
                    ),
                    Ok(Solution::Determined(())),
                    "{match_algorithm:?}/{subgraph_isomorphism_algorithm:?}",
                );
            }
        }
    }

    #[rstest]
    fn test_molecule_constraint_validator_validate_subpattern_absent() {
        let host = mol_dsl!(r#"{:atoms ["C"] :bonds []}"#);
        let constraint = MoleculeConstraint::SubPattern {
            anchor: SubPatternAnchor::new(),
            pattern: Box::new(mol_dsl!(r#"{:atoms ["O"] :bonds []}"#)),
        };

        for match_algorithm in SUBSTRUCTURE_MATCH_ALGORITHMS {
            for subgraph_isomorphism_algorithm in SUBGRAPH_ISOMORPHISM_ALGORITHMS {
                assert_eq!(
                    MoleculeConstraintValidator.validate(
                        &host,
                        &constraint,
                        ConstraintValidateConfig {
                            substructure_match_algorithm: match_algorithm,
                            subgraph_isomorphism_algorithm,
                            ..CONFIG
                        },
                    ),
                    Ok(Solution::Contradictory(MoleculeConstraintContradiction {
                        constraint: constraint.clone(),
                    },)),
                    "{match_algorithm:?}/{subgraph_isomorphism_algorithm:?}",
                );
            }
        }
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#), { let mut anchor = SubPatternAnchor::new(); anchor.push_atom(AtomId(0), AtomId(0)); anchor })]
    #[case::bond(mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#), { let mut anchor = SubPatternAnchor::new(); anchor.push_bond(BondId(0), BondId(0)); anchor })]
    #[case::dative_bond(mol_dsl!(r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1"}]}"#), { let mut anchor = SubPatternAnchor::new(); anchor.push_dative_bond(DativeBondId(0), DativeBondId(0)); anchor })]
    #[case::aromatic_system(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*#e2"}]}"#), { let mut anchor = SubPatternAnchor::new(); anchor.push_aromatic_system(AromaticSystemId(0), AromaticSystemId(0)); anchor })]
    #[case::multicenter_bond(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1] :type "*#e2"}]}"#), { let mut anchor = SubPatternAnchor::new(); anchor.push_multicenter_bond(MulticenterBondId(0), MulticenterBondId(0)); anchor })]
    #[case::noncovalent_bond(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"#), { let mut anchor = SubPatternAnchor::new(); anchor.push_noncovalent_bond(NoncovalentBondId(0), NoncovalentBondId(0)); anchor })]
    #[case::stereo_atom(mol_dsl!(r#"{:atoms ["C #h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]}"#), { let mut anchor = SubPatternAnchor::new(); anchor.push_stereo_atom(StereoAtomId(0), StereoAtomId(0)); anchor })]
    #[case::stereo_bond(mol_dsl!(r#"{:atoms ["F" "Cl" "C" "N" "Br" "I"] :bonds [[2 3 "2"]] :stereo-bonds [{:site 0 :ligands [0 1 4 5] :type "Ct1"}]}"#), { let mut anchor = SubPatternAnchor::new(); anchor.push_stereo_bond(StereoBondId(0), StereoBondId(0)); anchor })]
    fn test_molecule_constraint_validator_validate_subpattern_anchor(
        #[case] molecule: MoleculeAst,
        #[case] anchor: SubPatternAnchor,
    ) {
        let constraint = MoleculeConstraint::SubPattern {
            anchor,
            pattern: Box::new(molecule.clone()),
        };

        for match_algorithm in SUBSTRUCTURE_MATCH_ALGORITHMS {
            for subgraph_isomorphism_algorithm in SUBGRAPH_ISOMORPHISM_ALGORITHMS {
                assert_eq!(
                    MoleculeConstraintValidator.validate(
                        &molecule,
                        &constraint,
                        ConstraintValidateConfig {
                            substructure_match_algorithm: match_algorithm,
                            subgraph_isomorphism_algorithm,
                            ..CONFIG
                        },
                    ),
                    Ok(Solution::Determined(())),
                    "{match_algorithm:?}/{subgraph_isomorphism_algorithm:?}",
                );
            }
        }
    }

    #[rstest]
    fn test_molecule_constraint_validator_validate_subpattern_anchor_contradiction() {
        let host = mol_dsl!(r#"{:atoms ["C" "O"] :bonds []}"#);
        let pattern = mol_dsl!(r#"{:atoms ["O"] :bonds []}"#);
        let mut anchor = SubPatternAnchor::new();
        anchor.push_atom(AtomId(0), AtomId(0));
        let constraint = MoleculeConstraint::SubPattern {
            anchor,
            pattern: Box::new(pattern),
        };

        for match_algorithm in SUBSTRUCTURE_MATCH_ALGORITHMS {
            for subgraph_isomorphism_algorithm in SUBGRAPH_ISOMORPHISM_ALGORITHMS {
                assert_eq!(
                    MoleculeConstraintValidator.validate(
                        &host,
                        &constraint,
                        ConstraintValidateConfig {
                            substructure_match_algorithm: match_algorithm,
                            subgraph_isomorphism_algorithm,
                            ..CONFIG
                        },
                    ),
                    Ok(Solution::Contradictory(MoleculeConstraintContradiction {
                        constraint: constraint.clone(),
                    },)),
                    "{match_algorithm:?}/{subgraph_isomorphism_algorithm:?}",
                );
            }
        }
    }

    #[rstest]
    fn test_molecule_constraint_validator_validate_subpattern_ring_constraint() {
        let host = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 2 "1"] [1 3 "1"] [2 3 "1"]]}"#
        );
        let constraint = MoleculeConstraint::SubPattern {
            anchor: SubPatternAnchor::new(),
            pattern: Box::new(mol_dsl!(r#"{:atoms ["C#R3"] :bonds []}"#)),
        };

        for match_algorithm in SUBSTRUCTURE_MATCH_ALGORITHMS {
            for subgraph_isomorphism_algorithm in SUBGRAPH_ISOMORPHISM_ALGORITHMS {
                assert_eq!(
                    MoleculeConstraintValidator.validate(
                        &host,
                        &constraint,
                        ConstraintValidateConfig {
                            substructure_match_algorithm: match_algorithm,
                            subgraph_isomorphism_algorithm,
                            ..CONFIG
                        },
                    ),
                    Ok(Solution::Determined(())),
                    "{match_algorithm:?}/{subgraph_isomorphism_algorithm:?}",
                );
            }
        }
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::target_atom({ let mut anchor = SubPatternAnchor::new(); anchor.push_atom(AtomId(99), AtomId(0)); anchor }, Entity::Atom(AtomId(99)))]
    #[case::target_bond({ let mut anchor = SubPatternAnchor::new(); anchor.push_bond(BondId(99), BondId(0)); anchor }, Entity::Bond(BondId(99)))]
    #[case::target_dative_bond({ let mut anchor = SubPatternAnchor::new(); anchor.push_dative_bond(DativeBondId(99), DativeBondId(0)); anchor }, Entity::DativeBond(DativeBondId(99)))]
    #[case::target_aromatic_system({ let mut anchor = SubPatternAnchor::new(); anchor.push_aromatic_system(AromaticSystemId(99), AromaticSystemId(0)); anchor }, Entity::AromaticSystem(AromaticSystemId(99)))]
    #[case::target_multicenter_bond({ let mut anchor = SubPatternAnchor::new(); anchor.push_multicenter_bond(MulticenterBondId(99), MulticenterBondId(0)); anchor }, Entity::MulticenterBond(MulticenterBondId(99)))]
    #[case::target_noncovalent_bond({ let mut anchor = SubPatternAnchor::new(); anchor.push_noncovalent_bond(NoncovalentBondId(99), NoncovalentBondId(0)); anchor }, Entity::NoncovalentBond(NoncovalentBondId(99)))]
    #[case::target_stereo_atom({ let mut anchor = SubPatternAnchor::new(); anchor.push_stereo_atom(StereoAtomId(99), StereoAtomId(0)); anchor }, Entity::StereoAtom(StereoAtomId(99)))]
    #[case::target_stereo_bond({ let mut anchor = SubPatternAnchor::new(); anchor.push_stereo_bond(StereoBondId(99), StereoBondId(0)); anchor }, Entity::StereoBond(StereoBondId(99)))]
    #[case::pattern_atom({ let mut anchor = SubPatternAnchor::new(); anchor.push_atom(AtomId(0), AtomId(99)); anchor }, Entity::Atom(AtomId(99)))]
    fn test_molecule_constraint_validator_validate_subpattern_anchor_error(
        anchored_molecule: MoleculeAst,
        #[case] anchor: SubPatternAnchor,
        #[case] entity: Entity,
    ) {
        let constraint = MoleculeConstraint::SubPattern {
            anchor,
            pattern: Box::new(anchored_molecule.clone()),
        };

        assert_eq!(
            MoleculeConstraintValidator.validate(
                &anchored_molecule,
                &constraint,
                CONFIG,
            ),
            Err(ConstraintError::InvalidReference { entity }),
        );
    }
}
