//! Aromaticity resolver. Perception reads aromatic valence constraints from the
//! materialized valence stage; planning emits complete aromatic systems and
//! localized bond aromatic constraints without mutating the source molecule.

use std::collections::BTreeSet;

use umol_ast::ast::{
    AromaticSystemAst, AromaticValenceAst, AtomConstraintAst, AtomHandle, AtomId, AtomUpdate,
    BondConstraintAst, BondHandle, BondUpdate, BooleanAst, Edit, MoleculeAst,
};
use umol_utils::solution::Solution;

use crate::ops::aromaticity::{
    AromaticityConfig, AromaticityContradiction, AromaticityError, AromaticityPerception,
};
use crate::ops::model::AromaticityModel;

/// How aromaticity resolution handles assertions that disagree with perception.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AromaticityInconsistencyPolicy {
    Keep,
    Error,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AromaticityResolveConfig {
    pub perception: AromaticityConfig,
    pub inconsistency: AromaticityInconsistencyPolicy,
    pub reset_aromatic_valence: bool,
}

impl Default for AromaticityResolveConfig {
    fn default() -> Self {
        Self {
            perception: AromaticityConfig::default(),
            inconsistency: AromaticityInconsistencyPolicy::Error,
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
    ) -> Result<Solution<Vec<Edit>, AromaticityContradiction>, AromaticityError> {
        let outcome = self.perception.derive(ast, self.config.perception)?;

        match outcome {
            Solution::Determined(derivation) => {
                if self.config.inconsistency == AromaticityInconsistencyPolicy::Error {
                    if let Some(&mismatch) = derivation.mismatches.first() {
                        return Ok(Solution::Contradictory(mismatch.into()));
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
                let mut edits = Vec::new();
                for (atoms, system) in derivation.systems {
                    let mut key = atoms.clone();
                    key.sort_unstable();
                    if !existing.contains(&key) {
                        edits.extend(self.plan_system(ast, atoms, system));
                    }
                }
                Ok(Solution::Determined(edits))
            }
            Solution::Underdetermined(_) => Ok(Solution::Underdetermined(Vec::new())),
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
        system: AromaticSystemAst,
    ) -> Vec<Edit> {
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

        let mut edits = vec![Edit::AddAromaticSystem {
            atoms: atoms.iter().copied().map(AtomHandle::Id).collect(),
            ast: system,
        }];
        for (atom_id, update) in atom_updates {
            edits.extend(Edit::for_atom_update(
                AtomHandle::Id(atom_id),
                ast.atom(atom_id).ast,
                &update,
            ));
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
            let mut update = BondUpdate::default();
            update
                .constraints
                .set(BondConstraintAst::Aromatic(BooleanAst::Lit(true)));
            edits.extend(Edit::for_bond_update(
                BondHandle::Id(bond_id),
                ast.bond(bond_id).ast,
                &update,
            ));
        }
        edits
    }
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};
    use umol_ast::ast::{
        AromaticSystemId, BondConstraintKey, BondId, RingConfig, UnpairedElectronsAst, ValueAst,
    };
    use umol_ast::{mol_dsl, mol_dsl_ground};
    use umol_graph_core::{
        ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm,
        RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm,
    };

    use super::*;
    use crate::ops::aromaticity::AromaticityMismatch;
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

    #[rstest]
    fn test_aromaticity_resolve_config_default() {
        assert_eq!(
            AromaticityResolveConfig::default(),
            AromaticityResolveConfig {
                perception: AromaticityConfig::default(),
                inconsistency: AromaticityInconsistencyPolicy::Error,
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
                    inconsistency: AromaticityInconsistencyPolicy::Error,
                    reset_aromatic_valence: false,
                },
            )
            .plan(&benzene),
            Ok(Solution::Determined(vec![
                Edit::AddAromaticSystem {
                    atoms: (0..6).map(|id| AtomHandle::Id(AtomId(id))).collect(),
                    ast: AromaticSystemAst::from_electrons(vec![1; 6])
                        .with_charge(0)
                        .with_unpaired_electrons(UnpairedElectronsAst::closed_shell()),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(0)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(1)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(2)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(3)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(4)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(5)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
                },
            ]))
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
            perception: AromaticityConfig::default(),
            inconsistency: AromaticityInconsistencyPolicy::Keep,
            reset_aromatic_valence: false,
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
            perception: AromaticityConfig::default(),
            inconsistency: AromaticityInconsistencyPolicy::Keep,
            reset_aromatic_valence: false,
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
            Ok(Solution::Determined(Vec::new()))
        );
    }

    #[rstest]
    #[case::homogeneous_localized(
        AromaticityResolveConfig::default(),
        mol_dsl_ground!(r#"{:atoms ["C #h #a" "C #h #a" "C #c+ #h #a0"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#),
        ValueAst::Lit(0),
        vec![ValueAst::Lit(0), ValueAst::Lit(0), ValueAst::Lit(1)],
        vec![
            Some(AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
            Some(AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
            Some(AromaticValenceAst::Aromatic(ValueAst::Lit(0))),
        ]
    )]
    #[case::heterogeneous_localized(
        AromaticityResolveConfig::default(),
        mol_dsl_ground!(r#"{:atoms ["N #c+ #h #a" "C #h #a" "C #h #a"
                                      "C #h #a" "C #h #a" "C #h #a"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                                      [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        ValueAst::Lit(0),
        vec![
            ValueAst::Lit(1), ValueAst::Lit(0), ValueAst::Lit(0),
            ValueAst::Lit(0), ValueAst::Lit(0), ValueAst::Lit(0),
        ],
        vec![Some(AromaticValenceAst::Aromatic(ValueAst::Lit(1))); 6]
    )]
    #[case::accepted_system_with_rejected_projections(
        AromaticityResolveConfig {
            perception: AromaticityConfig::default(),
            inconsistency: AromaticityInconsistencyPolicy::Keep,
            reset_aromatic_valence: false,
        },
        mol_dsl_ground!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"
                    "C#h3#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#),
        ValueAst::Lit(0),
        vec![ValueAst::Lit(0); 7],
        vec![
            Some(AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
            Some(AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
            Some(AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
            Some(AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
            Some(AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
            Some(AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
            Some(AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
        ]
    )]
    #[case::reset_source_constraints(
        AromaticityResolveConfig {
            perception: AromaticityConfig::default(),
            inconsistency: AromaticityInconsistencyPolicy::Error,
            reset_aromatic_valence: true,
        },
        benzene(),
        ValueAst::Lit(0),
        vec![ValueAst::Lit(0); 6],
        vec![None; 6]
    )]
    fn test_aromaticity_resolver_resolve(
        aromaticity_model: AromaticityModel,
        #[case] config: AromaticityResolveConfig,
        #[case] mut molecule: MoleculeAst,
        #[case] expected_system_charge: ValueAst,
        #[case] expected_atom_charges: Vec<ValueAst>,
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
            Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true)))
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
    #[case::projection_mismatch(
        AromaticityModel::mdl(),
        mol_dsl_ground!(r#"{:atoms ["O #n1 #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        AromaticityContradiction::Mismatch(AromaticityMismatch::AtomProjection {
            atom: AtomId(0),
        })
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
