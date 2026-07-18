//! Aromaticity resolver. Perception reads aromatic valence constraints from the
//! materialized valence stage; planning emits complete aromatic systems,
//! optional charge delocalization atom updates, and localized bond aromatic
//! constraints without mutating the source molecule.

use std::collections::{BTreeMap, BTreeSet};

use umol_ast::ast::{
    AromaticSystemAst, AromaticValenceAst, AtomConstraintAst, AtomHandle, AtomId, AtomUpdate,
    BondConstraintAst, BondHandle, BondUpdate, BooleanAst, Edit, MoleculeAst, ValueAst,
};
use umol_utils::solution::Solution;

use crate::ops::aromaticity::{
    derive_charge_equalization, AromaticityContradiction, AromaticityError, AromaticityPerception,
};
use crate::ops::model::AromaticityModel;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AromaticityResolverConfig {
    pub delocalize_charge: bool,
    pub reset_aromatic_valence: bool,
}

impl Default for AromaticityResolverConfig {
    fn default() -> Self {
        Self {
            delocalize_charge: true,
            reset_aromatic_valence: false,
        }
    }
}

#[derive(Clone, Debug)]
pub struct AromaticityResolver {
    perception: AromaticityPerception,
    config: AromaticityResolverConfig,
}

impl AromaticityResolver {
    pub fn new(model: &AromaticityModel) -> Self {
        Self::with_config(model, AromaticityResolverConfig::default())
    }

    pub fn with_config(model: &AromaticityModel, config: AromaticityResolverConfig) -> Self {
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
        let outcome = self.perception.find_systems(ast, |view| {
            if view.is_in_aromatic_system() {
                return None;
            }
            match view
                .ast
                .constraints
                .aromatic_valence()
                .unwrap_or(&AromaticValenceAst::Undetermined)
            {
                AromaticValenceAst::Aromatic(ValueAst::Lit(n)) if *n >= 0 => Some(*n as u8),
                _ => None,
            }
        })?;

        match outcome {
            Solution::Determined(systems) => {
                let mut edits = Vec::new();
                for (atoms, system) in systems {
                    edits.extend(self.plan_system(ast, atoms, system)?);
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
    ) -> Result<Vec<Edit>, AromaticityError> {
        let (system, charge_updates) = if self.config.delocalize_charge {
            derive_charge_equalization(ast, &atoms, &system)?
        } else {
            (system, Vec::new())
        };
        let mut atom_updates: BTreeMap<AtomId, AtomUpdate> = charge_updates.into_iter().collect();
        if self.config.reset_aromatic_valence {
            for &atom_id in &atoms {
                atom_updates.entry(atom_id).or_default().constraints.set(
                    AtomConstraintAst::AromaticValence(AromaticValenceAst::Undetermined),
                );
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
        Ok(edits)
    }
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};
    use umol_ast::ast::{AromaticSystemId, BondConstraintKey, BondId, SpinStateAst};
    use umol_ast::{mol_dsl, mol_dsl_ground};

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

    #[rstest]
    fn test_aromaticity_resolver_config_default() {
        assert_eq!(
            AromaticityResolverConfig::default(),
            AromaticityResolverConfig {
                delocalize_charge: true,
                reset_aromatic_valence: false,
            }
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_plan(aromaticity_model: AromaticityModel, benzene: MoleculeAst) {
        assert_eq!(
            AromaticityResolver::new(&aromaticity_model).plan(&benzene),
            Ok(Solution::Determined(vec![
                Edit::AddAromaticSystem {
                    atoms: (0..6).map(|id| AtomHandle::Id(AtomId(id))).collect(),
                    ast: AromaticSystemAst::from_electrons(vec![1; 6])
                        .with_charge(0)
                        .with_spin(SpinStateAst::closed_shell()),
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
    fn test_aromaticity_resolver_plan_identity(aromaticity_model: AromaticityModel) {
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
        }"#
        );
        assert_eq!(
            AromaticityResolver::new(&aromaticity_model).plan(&molecule),
            Ok(Solution::Determined(Vec::new()))
        );
    }

    #[rstest]
    #[case::homogeneous_delocalized(
        AromaticityResolverConfig::default(),
        mol_dsl_ground!(r#"{:atoms ["C #h #a" "C #h #a" "C #c+ #h #a0"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#),
        ValueAst::Lit(1),
        vec![ValueAst::Lit(0); 3],
        vec![Some(AromaticValenceAst::Aromatic(ValueAst::Lit(1))); 3]
    )]
    #[case::homogeneous_localized(
        AromaticityResolverConfig {
            delocalize_charge: false,
            reset_aromatic_valence: false,
        },
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
        AromaticityResolverConfig::default(),
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
    #[case::reset_source_constraints(
        AromaticityResolverConfig {
            delocalize_charge: true,
            reset_aromatic_valence: true,
        },
        benzene(),
        ValueAst::Lit(0),
        vec![ValueAst::Lit(0); 6],
        vec![None; 6]
    )]
    fn test_aromaticity_resolver_resolve(
        aromaticity_model: AromaticityModel,
        #[case] config: AromaticityResolverConfig,
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

    #[rstest]
    fn test_aromaticity_resolver_resolve_error(
        aromaticity_model: AromaticityModel,
        mut benzene: MoleculeAst,
    ) {
        benzene.atom_mut(AtomId(0)).ast.charge = ValueAst::Undetermined;
        let original = benzene.clone();
        assert_eq!(
            AromaticityResolver::new(&aromaticity_model).resolve(&mut benzene),
            Err(AromaticityError::NonGroundAtom(AtomId(0)))
        );
        assert_eq!(benzene, original);
    }
}
