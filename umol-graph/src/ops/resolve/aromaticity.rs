//! Aromaticity resolver. Perception reads aromatic valence constraints from the
//! materialized valence stage; planning emits complete aromatic systems and
//! localized bond aromatic constraints without mutating the source molecule.

use std::collections::BTreeSet;

use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemHandle, AromaticSystemId, AromaticValenceAst,
    AtomConstraintAst, AtomHandle, AtomId, AtomUpdate, BondConstraintAst, BondHandle, BondUpdate,
    BooleanForm, Edits, MoleculeAst,
};
use umol_utils::solution::Solution;

use crate::ops::aromaticity::{
    AromaticityConfig, AromaticityContradiction, AromaticityError, AromaticityInconsistency,
    AromaticityPerception,
};
use crate::ops::model::AromaticityModel;

/// How aromaticity resolution handles an independently invalid constraint or entity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AromaticityFailurePolicy {
    Error,
    Keep,
}

/// How aromaticity resolution handles a valid aromatic-valence constraint that disagrees with a
/// valid aromatic system.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AromaticityMismatchPolicy {
    Error,
    Keep,
    RemoveConstraint,
    ReplaceEntity,
}

/// How aromaticity resolution handles a valid localized-bond aromatic constraint that disagrees
/// with a valid aromatic system.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AromaticBondConstraintMismatchPolicy {
    Error,
    Keep,
    RemoveConstraint,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AromaticityResolveConfig {
    pub perception: AromaticityConfig,
    pub aromatic_valence_failure: AromaticityFailurePolicy,
    pub aromatic_system_failure: AromaticityFailurePolicy,
    pub aromatic_valence_mismatch: AromaticityMismatchPolicy,
    pub aromatic_bond_constraint_mismatch: AromaticBondConstraintMismatchPolicy,
    pub reset_aromatic_valence: bool,
}

impl Default for AromaticityResolveConfig {
    fn default() -> Self {
        Self {
            perception: AromaticityConfig::default(),
            aromatic_valence_failure: AromaticityFailurePolicy::Error,
            aromatic_system_failure: AromaticityFailurePolicy::Error,
            aromatic_valence_mismatch: AromaticityMismatchPolicy::Error,
            aromatic_bond_constraint_mismatch: AromaticBondConstraintMismatchPolicy::Error,
            reset_aromatic_valence: false,
        }
    }
}

#[derive(Clone, Debug)]
pub struct AromaticityResolver {
    perception: AromaticityPerception,
    config: AromaticityResolveConfig,
}

impl AromaticityResolver {
    pub fn new(model: &AromaticityModel) -> Self {
        Self::with_config(model, AromaticityResolveConfig::default())
    }

    pub fn with_config(model: &AromaticityModel, config: AromaticityResolveConfig) -> Self {
        Self {
            perception: AromaticityPerception::new(model),
            config,
        }
    }

    /// Construct the complete aromaticity edit plan without mutating `ast`.
    pub fn plan(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<Edits, AromaticityContradiction>, AromaticityError> {
        let outcome = self.perception.derive(ast, self.config.perception)?;

        match outcome {
            Solution::Determined(derivation) => {
                for &inconsistency in &derivation.inconsistencies {
                    let error = match inconsistency {
                        AromaticityInconsistency::AromaticValenceFailure { .. } => {
                            self.config.aromatic_valence_failure == AromaticityFailurePolicy::Error
                        }
                        AromaticityInconsistency::AromaticSystemFailure { .. } => {
                            self.config.aromatic_system_failure == AromaticityFailurePolicy::Error
                        }
                        AromaticityInconsistency::AromaticValenceMismatch { .. } => {
                            self.config.aromatic_valence_mismatch
                                == AromaticityMismatchPolicy::Error
                        }
                        AromaticityInconsistency::AromaticBondConstraintMismatch { .. } => {
                            self.config.aromatic_bond_constraint_mismatch
                                == AromaticBondConstraintMismatchPolicy::Error
                        }
                    };
                    if error {
                        return Ok(Solution::Contradictory(inconsistency.into()));
                    }
                }

                let existing: BTreeSet<Vec<AtomId>> = ast
                    .aromatic_systems()
                    .iter()
                    .map(|system| {
                        let mut atoms: Vec<AtomId> = system.atom_ids().collect();
                        atoms.sort_unstable();
                        atoms
                    })
                    .collect();

                let mut edits = Edits::new();
                let mut remove_constraints = BTreeSet::new();
                let mut remove_bond_constraints = BTreeSet::new();
                let mut replacements = BTreeSet::new();
                let mut suppressed = BTreeSet::new();

                for inconsistency in derivation.inconsistencies {
                    match inconsistency {
                        AromaticityInconsistency::AromaticValenceMismatch { atom, system } => {
                            let existing_members: BTreeSet<AtomId> =
                                ast.aromatic_system(system).atom_ids().collect();
                            let candidate = derivation.systems.iter().position(|(atoms, _)| {
                                atoms.iter().copied().collect::<BTreeSet<_>>() == existing_members
                            });
                            match self.config.aromatic_valence_mismatch {
                                AromaticityMismatchPolicy::Error => unreachable!(),
                                AromaticityMismatchPolicy::Keep => {
                                    if let Some(candidate) = candidate {
                                        suppressed.insert(candidate);
                                    }
                                }
                                AromaticityMismatchPolicy::RemoveConstraint => {
                                    remove_constraints.insert(atom);
                                    if let Some(candidate) = candidate {
                                        suppressed.insert(candidate);
                                    }
                                }
                                AromaticityMismatchPolicy::ReplaceEntity => {
                                    if let Some(candidate) = candidate {
                                        replacements.insert((system, candidate));
                                    }
                                }
                            }
                        }
                        AromaticityInconsistency::AromaticBondConstraintMismatch {
                            bond, ..
                        } => match self.config.aromatic_bond_constraint_mismatch {
                            AromaticBondConstraintMismatchPolicy::Error => unreachable!(),
                            AromaticBondConstraintMismatchPolicy::Keep => {}
                            AromaticBondConstraintMismatchPolicy::RemoveConstraint => {
                                remove_bond_constraints.insert(bond);
                            }
                        },
                        AromaticityInconsistency::AromaticValenceFailure { .. }
                        | AromaticityInconsistency::AromaticSystemFailure { .. } => {}
                    }
                }

                if !replacements.is_empty() {
                    let removes = replacements
                        .iter()
                        .map(|&(system, _)| {
                            let view = ast.aromatic_system(system);
                            (
                                AromaticSystemHandle::Id(system),
                                view.atom_ids().map(AtomHandle::Id).collect(),
                                view.ast.clone(),
                            )
                        })
                        .collect();
                    edits.remove_aromatic_systems(removes);
                }

                for atom in remove_constraints {
                    let mut update = AtomUpdate::default();
                    update.constraints.set(AtomConstraintAst::AromaticValence(
                        AromaticValenceAst::Undetermined,
                    ));
                    edits.update_atom(AtomHandle::Id(atom), ast.atom(atom).ast, &update);
                }
                for bond in remove_bond_constraints {
                    let mut update = BondUpdate::default();
                    update
                        .constraints
                        .set(BondConstraintAst::Aromatic(BooleanForm::Undetermined));
                    edits.update_bond(BondHandle::Id(bond), ast.bond(bond).ast, &update);
                }

                let replaced_candidates: BTreeSet<usize> = replacements
                    .iter()
                    .map(|&(_, candidate)| candidate)
                    .collect();
                let replaced_entities: BTreeSet<AromaticSystemId> =
                    replacements.iter().map(|&(system, _)| system).collect();
                let retained_existing: BTreeSet<Vec<AtomId>> = ast
                    .aromatic_systems()
                    .iter()
                    .filter(|system| !replaced_entities.contains(&system.id))
                    .map(|system| {
                        let mut atoms: Vec<AtomId> = system.atom_ids().collect();
                        atoms.sort_unstable();
                        atoms
                    })
                    .collect();

                for (candidate, (atoms, system)) in derivation.systems.into_iter().enumerate() {
                    let mut key = atoms.clone();
                    key.sort_unstable();
                    if replaced_candidates.contains(&candidate)
                        || (!suppressed.contains(&candidate)
                            && !existing.contains(&key)
                            && !retained_existing.contains(&key))
                    {
                        for edit in self.plan_system(ast, atoms, system) {
                            edits.push(edit);
                        }
                    }
                }
                Ok(Solution::Determined(edits))
            }
            Solution::Underdetermined(_) => Ok(Solution::Underdetermined(Edits::new())),
            Solution::Contradictory(contradiction) => Ok(Solution::Contradictory(contradiction)),
        }
    }

    /// Plan and atomically apply aromaticity resolution.
    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), AromaticityContradiction>, AromaticityError> {
        let edits = match self.plan(ast)? {
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

    fn plan_system(
        &self,
        ast: &MoleculeAst,
        atoms: Vec<AtomId>,
        system: AromaticSystemForm,
    ) -> Edits {
        let mut atom_updates = Vec::new();
        if self.config.reset_aromatic_valence {
            for &atom_id in &atoms {
                let mut update = AtomUpdate::default();
                update.constraints.set(AtomConstraintAst::AromaticValence(
                    AromaticValenceAst::Undetermined,
                ));
                atom_updates.push((atom_id, update));
            }
        }

        let mut edits = Edits::new();
        edits.add_aromatic_system(atoms.iter().copied().map(AtomHandle::Id).collect(), system);
        for (atom_id, update) in atom_updates {
            edits.update_atom(AtomHandle::Id(atom_id), ast.atom(atom_id).ast, &update);
        }

        let members: BTreeSet<AtomId> = atoms.iter().copied().collect();
        let mut bond_ids = BTreeSet::new();
        for &atom_id in &atoms {
            for neighbor in ast.atom(atom_id).neighbors() {
                if members.contains(&neighbor.atom_id()) {
                    bond_ids.insert(neighbor.bond_id());
                }
            }
        }
        for bond_id in bond_ids {
            if matches!(
                ast.bond(bond_id).ast.constraints.aromatic(),
                BooleanForm::Lit(_)
            ) {
                continue;
            }
            let mut update = BondUpdate::default();
            update
                .constraints
                .set(BondConstraintAst::Aromatic(BooleanForm::Lit(true)));
            edits.update_bond(BondHandle::Id(bond_id), ast.bond(bond_id).ast, &update);
        }
        edits
    }
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};
    use umol_graph_core::{
        ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm,
        RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm,
    };
    use umol_graph_ir::ir::{
        AromaticSystemId, BondConstraintKey, BondId, Edit, Edits, NumForm, RingConfig,
        UnpairedElectronsForm,
    };
    use umol_graph_ir::{mol_dsl, mol_dsl_ground};

    use super::*;
    use crate::ops::model::{ElementScope, RingLimits};

    #[fixture]
    fn aromaticity_model() -> AromaticityModel {
        AromaticityModel::HueckelRule {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        }
    }

    #[fixture]
    fn benzene() -> MoleculeAst {
        mol_dsl!(
            r#"{
            :atoms ["C#i=#c0#h#n0#u0#s#v2#a" "C#i=#c0#h#n0#u0#s#v2#a"
                    "C#i=#c0#h#n0#u0#s#v2#a" "C#i=#c0#h#n0#u0#s#v2#a"
                    "C#i=#c0#h#n0#u0#s#v2#a" "C#i=#c0#h#n0#u0#s#v2#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#
        )
    }

    #[fixture]
    fn aromatic_valence_mismatch() -> MoleculeAst {
        mol_dsl!(
            r#"{
            :atoms ["C#a2" "C#a0" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
        }"#
        )
    }

    #[fixture]
    fn aromatic_bond_constraint_mismatch() -> MoleculeAst {
        mol_dsl!(
            r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a!"] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
        }"#
        )
    }

    #[rstest]
    fn test_aromaticity_resolve_config_default() {
        assert_eq!(
            AromaticityResolveConfig::default(),
            AromaticityResolveConfig {
                perception: AromaticityConfig::default(),
                aromatic_valence_failure: AromaticityFailurePolicy::Error,
                aromatic_system_failure: AromaticityFailurePolicy::Error,
                aromatic_valence_mismatch: AromaticityMismatchPolicy::Error,
                aromatic_bond_constraint_mismatch: AromaticBondConstraintMismatchPolicy::Error,
                reset_aromatic_valence: false,
            }
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_plan(aromaticity_model: AromaticityModel, benzene: MoleculeAst) {
        assert_eq!(
            AromaticityResolver::with_config(
                &aromaticity_model,
                AromaticityResolveConfig {
                    perception: AromaticityConfig {
                        ring_config: RingConfig {
                            simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm::ReadTarjan,
                            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                        },
                        connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                        maximum_independent_set_algorithm:
                            MaximumIndependentSetAlgorithm::BranchAndBound,
                    },
                    reset_aromatic_valence: false,
                    ..AromaticityResolveConfig::default()
                },
            )
            .plan(&benzene),
            Ok(Solution::Determined(Edits::from_iter([
                Edit::AddAromaticSystem {
                    atoms: (0..6).map(|id| AtomHandle::Id(AtomId(id))).collect(),
                    ast: AromaticSystemForm::from_electrons(vec![1; 6])
                        .with_charge(0)
                        .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(0)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanForm::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(1)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanForm::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(2)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanForm::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(3)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanForm::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(4)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanForm::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(5)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanForm::Lit(true))),
                },
            ])))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_plan_partial(aromaticity_model: AromaticityModel) {
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C#a+" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"]]
        }"#
        );

        assert_eq!(
            AromaticityResolver::new(&aromaticity_model).plan(&molecule),
            Ok(Solution::Underdetermined(Edits::new()))
        );
    }

    #[rstest]
    #[case::error(
        AromaticityMismatchPolicy::Error,
        Solution::Contradictory(AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticValenceMismatch {
                atom: AtomId(0),
                system: AromaticSystemId(0),
            }
        ))
    )]
    #[case::keep(AromaticityMismatchPolicy::Keep, Solution::Determined(Edits::new()))]
    #[case::remove_constraint(
        AromaticityMismatchPolicy::RemoveConstraint,
        Solution::Determined(Edits::from_iter([
            Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(0)),
                old: Some(AtomConstraintAst::AromaticValence(
                    AromaticValenceAst::Aromatic(NumForm::Lit(2)),
                )),
                new: None,
            },
            Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(1)),
                old: Some(AtomConstraintAst::AromaticValence(
                    AromaticValenceAst::Aromatic(NumForm::Lit(0)),
                )),
                new: None,
            },
        ]))
    )]
    #[case::replace_entity(
        AromaticityMismatchPolicy::ReplaceEntity,
        Solution::Determined(Edits::from_iter([
            Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    (0..6).map(|id| AtomHandle::Id(AtomId(id))).collect(),
                    AromaticSystemForm::from_electrons(vec![1; 6]),
                )],
            },
            Edit::AddAromaticSystem {
                atoms: (0..6).map(|id| AtomHandle::Id(AtomId(id))).collect(),
                ast: AromaticSystemForm::from_electrons(vec![2, 0, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            },
        ]))
    )]
    fn test_aromaticity_resolver_plan_aromatic_valence_mismatch(
        aromaticity_model: AromaticityModel,
        aromatic_valence_mismatch: MoleculeAst,
        #[case] policy: AromaticityMismatchPolicy,
        #[case] expected: Solution<Edits, AromaticityContradiction>,
    ) {
        let resolver = AromaticityResolver::with_config(
            &aromaticity_model,
            AromaticityResolveConfig {
                aromatic_valence_mismatch: policy,
                ..AromaticityResolveConfig::default()
            },
        );

        assert_eq!(resolver.plan(&aromatic_valence_mismatch), Ok(expected));
    }

    #[rstest]
    fn test_aromaticity_resolver_resolve_aromatic_valence_mismatch_reset(
        aromaticity_model: AromaticityModel,
        mut aromatic_valence_mismatch: MoleculeAst,
    ) {
        let resolver = AromaticityResolver::with_config(
            &aromaticity_model,
            AromaticityResolveConfig {
                aromatic_valence_mismatch: AromaticityMismatchPolicy::ReplaceEntity,
                reset_aromatic_valence: true,
                ..AromaticityResolveConfig::default()
            },
        );
        let expected = mol_dsl!(
            r#"{
            :atoms ["C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{
                :atoms [0 1 2 3 4 5]
                :type "[2,0,1,1,1,1]#c0#u0#s"
            }]
        }"#
        );

        assert_eq!(
            resolver.resolve(&mut aromatic_valence_mismatch),
            Ok(Solution::Determined(()))
        );
        assert_eq!(aromatic_valence_mismatch, expected);
    }

    #[rstest]
    #[case::error(
        AromaticBondConstraintMismatchPolicy::Error,
        Solution::Contradictory(AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticBondConstraintMismatch {
                bond: BondId(0),
                system: AromaticSystemId(0),
            }
        ))
    )]
    #[case::keep(
        AromaticBondConstraintMismatchPolicy::Keep,
        Solution::Determined(Edits::new())
    )]
    #[case::remove_constraint(
        AromaticBondConstraintMismatchPolicy::RemoveConstraint,
        Solution::Determined(Edits::from_iter([Edit::ModifyBondConstraint {
            id: BondHandle::Id(BondId(0)),
            old: Some(BondConstraintAst::Aromatic(BooleanForm::Lit(false))),
            new: None,
        }]))
    )]
    fn test_aromaticity_resolver_plan_aromatic_bond_constraint_mismatch(
        aromaticity_model: AromaticityModel,
        aromatic_bond_constraint_mismatch: MoleculeAst,
        #[case] policy: AromaticBondConstraintMismatchPolicy,
        #[case] expected: Solution<Edits, AromaticityContradiction>,
    ) {
        let resolver = AromaticityResolver::with_config(
            &aromaticity_model,
            AromaticityResolveConfig {
                aromatic_bond_constraint_mismatch: policy,
                ..AromaticityResolveConfig::default()
            },
        );

        assert_eq!(
            resolver.plan(&aromatic_bond_constraint_mismatch),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::conformant(
        AromaticityModel::daylight(),
        AromaticityResolveConfig::default(),
        mol_dsl!(r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
        }"#)
    )]
    #[case::rejected_projections(
        AromaticityModel::mdl(),
        AromaticityResolveConfig {
            aromatic_valence_failure: AromaticityFailurePolicy::Keep,
            ..AromaticityResolveConfig::default()
        },
        mol_dsl!(r#"{
            :atoms ["O#n1#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 0 :aromatic]]
        }"#)
    )]
    #[case::rejected_existing_system(
        AromaticityModel::mdl(),
        AromaticityResolveConfig {
            aromatic_valence_failure: AromaticityFailurePolicy::Keep,
            aromatic_system_failure: AromaticityFailurePolicy::Keep,
            ..AromaticityResolveConfig::default()
        },
        mol_dsl!(r#"{
            :atoms ["O#n1#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :type "[2,1,1,1,1]"}]
        }"#)
    )]
    fn test_aromaticity_resolver_plan_identity(
        #[case] model: AromaticityModel,
        #[case] config: AromaticityResolveConfig,
        #[case] molecule: MoleculeAst,
    ) {
        assert_eq!(
            AromaticityResolver::with_config(&model, config).plan(&molecule),
            Ok(Solution::Determined(Edits::new()))
        );
    }

    #[rstest]
    #[case::homogeneous_localized(
        AromaticityResolveConfig::default(),
        mol_dsl_ground!(r#"{:atoms ["C #h #a" "C #h #a" "C #c+ #h #a0"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#),
        NumForm::Lit(0),
        vec![NumForm::Lit(0), NumForm::Lit(0), NumForm::Lit(1)],
        vec![
            Some(AromaticValenceAst::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceAst::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceAst::Aromatic(NumForm::Lit(0))),
        ]
    )]
    #[case::heterogeneous_localized(
        AromaticityResolveConfig::default(),
        mol_dsl_ground!(r#"{:atoms ["N #c+ #h #a" "C #h #a" "C #h #a"
                                      "C #h #a" "C #h #a" "C #h #a"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                                      [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        NumForm::Lit(0),
        vec![
            NumForm::Lit(1), NumForm::Lit(0), NumForm::Lit(0),
            NumForm::Lit(0), NumForm::Lit(0), NumForm::Lit(0),
        ],
        vec![Some(AromaticValenceAst::Aromatic(NumForm::Lit(1))); 6]
    )]
    #[case::accepted_system_with_rejected_projections(
        AromaticityResolveConfig {
            aromatic_valence_failure: AromaticityFailurePolicy::Keep,
            ..AromaticityResolveConfig::default()
        },
        mol_dsl_ground!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"
                    "C#h3#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#),
        NumForm::Lit(0),
        vec![NumForm::Lit(0); 7],
        vec![
            Some(AromaticValenceAst::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceAst::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceAst::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceAst::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceAst::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceAst::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceAst::Aromatic(NumForm::Lit(1))),
        ]
    )]
    #[case::reset_source_constraints(
        AromaticityResolveConfig {
            reset_aromatic_valence: true,
            ..AromaticityResolveConfig::default()
        },
        benzene(),
        NumForm::Lit(0),
        vec![NumForm::Lit(0); 6],
        vec![None; 6]
    )]
    fn test_aromaticity_resolver_resolve(
        aromaticity_model: AromaticityModel,
        #[case] config: AromaticityResolveConfig,
        #[case] mut molecule: MoleculeAst,
        #[case] expected_system_charge: NumForm,
        #[case] expected_atom_charges: Vec<NumForm>,
        #[case] expected_aromatic_valences: Vec<Option<AromaticValenceAst>>,
    ) {
        assert_eq!(
            AromaticityResolver::with_config(&aromaticity_model, config).resolve(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule.aromatic_systems().count(), 1);
        assert_eq!(
            molecule.aromatic_system(AromaticSystemId(0)).ast.charge,
            expected_system_charge
        );
        assert_eq!(
            molecule
                .atoms()
                .iter()
                .map(|atom| atom.ast.charge.clone())
                .collect::<Vec<_>>(),
            expected_atom_charges
        );
        assert_eq!(
            molecule
                .atoms()
                .iter()
                .map(|atom| atom.ast.constraints.aromatic_valence().cloned())
                .collect::<Vec<_>>(),
            expected_aromatic_valences
        );
        assert!(molecule.bonds().iter().all(|bond| matches!(
            bond.ast.constraints.get(BondConstraintKey::Aromatic),
            Some(BondConstraintAst::Aromatic(BooleanForm::Lit(true)))
        )));
    }

    #[rstest]
    fn test_aromaticity_resolver_resolve_identity(
        aromaticity_model: AromaticityModel,
        mut benzene: MoleculeAst,
    ) {
        let resolver = AromaticityResolver::new(&aromaticity_model);
        assert_eq!(resolver.resolve(&mut benzene), Ok(Solution::Determined(())));
        let expected = benzene.clone();

        assert_eq!(resolver.resolve(&mut benzene), Ok(Solution::Determined(())));
        assert_eq!(benzene, expected);
    }

    #[rstest]
    #[case::clar_heterocycle(
        AromaticityModel::Clar {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        },
        mol_dsl_ground!(r#"{:atoms ["N #h #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        AromaticityContradiction::ClarNonBenzenoid(
            "Clar model requires benzenoid input but non-carbon aromatic atoms are present".to_string()
        )
    )]
    #[case::aromatic_valence_failure(
        AromaticityModel::mdl(),
        mol_dsl_ground!(r#"{:atoms ["O #n1 #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) }
        )
    )]
    #[case::aromatic_system_failure(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :type "[1,1,1,1,1]"}]
        }"#),
        AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticSystemFailure {
                system: AromaticSystemId(0)
            }
        )
    )]
    fn test_aromaticity_resolver_resolve_contradiction(
        #[case] model: AromaticityModel,
        #[case] mut molecule: MoleculeAst,
        #[case] expected: AromaticityContradiction,
    ) {
        let original = molecule.clone();
        assert_eq!(
            AromaticityResolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Contradictory(expected))
        );
        assert_eq!(molecule, original);
    }
}
