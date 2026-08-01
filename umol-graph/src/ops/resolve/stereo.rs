//! Structural stereo resolver. Planning reads `#T` / `#C` assertions from the
//! materialized aromaticity state and emits stereo-element additions plus
//! optional source-constraint removals without mutating the molecule.

use std::collections::{BTreeMap, BTreeSet};

use thiserror::Error;
use umol_ast::ast::{
    AtomConstraintAst, AtomHandle, AtomId, AtomUpdate, BondConstraintAst, BondHandle, BondId,
    BondUpdate, CisTransStereoAst, Edit, MoleculeAst, TetrahedralStereoAst, TransactionError,
};
use umol_utils::solution::Solution;

use crate::ops::model::StereoModel;
use crate::ops::stereo::{StereoMismatch, StereoPerception};

/// How stereo resolution handles a `#T`/`#C` assertion it cannot fully
/// realize. The assertion is retained, removed, or reported as a contradiction;
/// it is never silently dropped.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StereoInconsistencyPolicy {
    Keep,
    Strip,
    Error,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StereoResolveConfig {
    pub reset_stereo_constraints: bool,
    pub inconsistency: StereoInconsistencyPolicy,
}

impl Default for StereoResolveConfig {
    fn default() -> Self {
        Self {
            reset_stereo_constraints: false,
            inconsistency: StereoInconsistencyPolicy::Error,
        }
    }
}

#[derive(Clone, Debug)]
pub struct StereoResolver {
    perception: StereoPerception,
    config: StereoResolveConfig,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoContradiction {
    #[error("tetrahedral stereo assertion at atom {0:?} cannot be realized")]
    UnrealizableAtom(AtomId),
    #[error("cis-trans stereo assertion at bond {0:?} cannot be realized")]
    UnrealizableBond(BondId),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoError {
    #[error(transparent)]
    Transaction(#[from] TransactionError),
}

impl StereoResolver {
    pub fn new(model: &StereoModel) -> Self {
        Self::with_config(model, StereoResolveConfig::default())
    }

    pub fn with_config(model: &StereoModel, config: StereoResolveConfig) -> Self {
        Self {
            perception: StereoPerception::new(model),
            config,
        }
    }

    /// Construct the complete stereo edit plan without mutating `ast`.
    pub fn plan(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<Vec<Edit>, StereoContradiction>, StereoError> {
        let derivation = self.perception.derive(ast);
        let mut atoms: BTreeMap<_, _> = derivation
            .atoms
            .into_iter()
            .map(|(id, ligands, stereo)| (id, (ligands, stereo)))
            .collect();
        let mut bonds: BTreeMap<_, _> = derivation
            .bonds
            .into_iter()
            .map(|(id, ligands, stereo)| (id, (ligands, stereo)))
            .collect();
        let unrealizable_atoms: BTreeSet<_> = derivation
            .mismatches
            .iter()
            .filter_map(|mismatch| match mismatch {
                StereoMismatch::UnrealizableAtom { atom } => Some(*atom),
                _ => None,
            })
            .collect();
        let unrealizable_bonds: BTreeSet<_> = derivation
            .mismatches
            .iter()
            .filter_map(|mismatch| match mismatch {
                StereoMismatch::UnrealizableBond { bond } => Some(*bond),
                _ => None,
            })
            .collect();

        let mut edits = Vec::new();
        for id in ast.atoms().ids() {
            if ast.stereo_atoms().is_at(id) {
                continue;
            }
            let Some((ligands, stereo)) = atoms.remove(&id) else {
                if !unrealizable_atoms.contains(&id) {
                    continue;
                }
                match self.config.inconsistency {
                    StereoInconsistencyPolicy::Keep => {}
                    StereoInconsistencyPolicy::Strip => {
                        let mut update = AtomUpdate::default();
                        update.constraints.set(AtomConstraintAst::TetrahedralStereo(
                            TetrahedralStereoAst::Undetermined,
                        ));
                        edits.extend(Edit::for_atom_update(
                            AtomHandle::Id(id),
                            ast.atom(id).ast,
                            &update,
                        ));
                    }
                    StereoInconsistencyPolicy::Error => {
                        return Ok(Solution::Contradictory(
                            StereoContradiction::UnrealizableAtom(id),
                        ));
                    }
                }
                continue;
            };
            edits.push(Edit::AddStereoAtom {
                site: AtomHandle::Id(id),
                ligands: ligands
                    .into_iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                ast: stereo,
            });
            if self.config.reset_stereo_constraints {
                let mut update = AtomUpdate::default();
                update.constraints.set(AtomConstraintAst::TetrahedralStereo(
                    TetrahedralStereoAst::Undetermined,
                ));
                edits.extend(Edit::for_atom_update(
                    AtomHandle::Id(id),
                    ast.atom(id).ast,
                    &update,
                ));
            }
        }
        for id in ast.bonds().ids() {
            if ast.stereo_bonds().is_at(id) {
                continue;
            }
            let Some((ligands, stereo)) = bonds.remove(&id) else {
                if !unrealizable_bonds.contains(&id) {
                    continue;
                }
                match self.config.inconsistency {
                    StereoInconsistencyPolicy::Keep => {}
                    StereoInconsistencyPolicy::Strip => {
                        let mut update = BondUpdate::default();
                        update.constraints.set(BondConstraintAst::CisTransStereo(
                            CisTransStereoAst::Undetermined,
                        ));
                        edits.extend(Edit::for_bond_update(
                            BondHandle::Id(id),
                            ast.bond(id).ast,
                            &update,
                        ));
                    }
                    StereoInconsistencyPolicy::Error => {
                        return Ok(Solution::Contradictory(
                            StereoContradiction::UnrealizableBond(id),
                        ));
                    }
                }
                continue;
            };
            edits.push(Edit::AddStereoBond {
                site: BondHandle::Id(id),
                ligands: ligands
                    .into_iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                ast: stereo,
            });
            if self.config.reset_stereo_constraints {
                let mut update = BondUpdate::default();
                update.constraints.set(BondConstraintAst::CisTransStereo(
                    CisTransStereoAst::Undetermined,
                ));
                edits.extend(Edit::for_bond_update(
                    BondHandle::Id(id),
                    ast.bond(id).ast,
                    &update,
                ));
            }
        }
        Ok(Solution::Determined(edits))
    }

    /// Plan and atomically apply structural stereo resolution.
    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), StereoContradiction>, StereoError> {
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
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};
    use umol_ast::ast::{
        AtomId, BondId, StereoAtomAst, StereoBondAst, StereoCoset, StereoKind, StereoLigandKind,
    };
    use umol_ast::mol_dsl_ground;

    use super::*;

    #[fixture]
    fn stereo_model() -> StereoModel {
        StereoModel::default()
    }

    #[rstest]
    fn test_stereo_resolve_config_default() {
        assert_eq!(
            StereoResolveConfig::default(),
            StereoResolveConfig {
                reset_stereo_constraints: false,
                inconsistency: StereoInconsistencyPolicy::Error,
            }
        );
    }

    #[rstest]
    #[case::tetrahedral(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
                             :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#),
        vec![Edit::AddStereoAtom {
            site: AtomHandle::Id(AtomId(1)),
            ligands: vec![
                (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(1)), StereoLigandKind::ImplicitHydrogen),
            ],
            ast: StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        }]
    )]
    #[case::cis_trans(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
                             :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]}"#),
        vec![Edit::AddStereoBond {
            site: BondHandle::Id(BondId(1)),
            ligands: vec![
                (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(1)), StereoLigandKind::ImplicitHydrogen),
                (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(2)), StereoLigandKind::ImplicitHydrogen),
            ],
            ast: StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        }]
    )]
    fn test_stereo_resolver_plan(
        stereo_model: StereoModel,
        #[case] molecule: MoleculeAst,
        #[case] expected: Vec<Edit>,
    ) {
        assert_eq!(
            StereoResolver::new(&stereo_model).plan(&molecule),
            Ok(Solution::Determined(expected))
        );
    }

    #[rstest]
    #[case::no_assertion(mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h3"] :bonds [[0 1 "1"]]}"#))]
    #[case::existing_atom(mol_dsl_ground!(r#"{
        :atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :type "Th1"}]
    }"#))]
    #[case::existing_atom_mismatch(mol_dsl_ground!(r#"{
        :atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :type "Th2"}]
    }"#))]
    #[case::existing_bond_mismatch(mol_dsl_ground!(r#"{
        :atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
        :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
        :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :type "Ct2"}]
    }"#))]
    fn test_stereo_resolver_plan_identity(
        stereo_model: StereoModel,
        #[case] molecule: MoleculeAst,
    ) {
        assert_eq!(
            StereoResolver::new(&stereo_model).plan(&molecule),
            Ok(Solution::Determined(Vec::new()))
        );
    }

    #[rstest]
    #[case::atom_keep(
        StereoInconsistencyPolicy::Keep,
        mol_dsl_ground!(
            r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ),
        Solution::Determined(Vec::new())
    )]
    #[case::atom_strip(
        StereoInconsistencyPolicy::Strip,
        mol_dsl_ground!(
            r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ),
        Solution::Determined(vec![Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(1)),
            old: Some(AtomConstraintAst::TetrahedralStereo(
                TetrahedralStereoAst::Stereo(StereoCoset::Lit(1)),
            )),
            new: None,
        }])
    )]
    #[case::atom_error(
        StereoInconsistencyPolicy::Error,
        mol_dsl_ground!(
            r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ),
        Solution::Contradictory(StereoContradiction::UnrealizableAtom(AtomId(1)))
    )]
    #[case::bond_keep(
        StereoInconsistencyPolicy::Keep,
        mol_dsl_ground!(
            r#"{:atoms ["C #h3" "C #h2" "C #h1"] :bonds [[0 1 "1"] [1 2 "2#C1"]]}"#
        ),
        Solution::Determined(Vec::new())
    )]
    #[case::bond_strip(
        StereoInconsistencyPolicy::Strip,
        mol_dsl_ground!(
            r#"{:atoms ["C #h3" "C #h2" "C #h1"] :bonds [[0 1 "1"] [1 2 "2#C1"]]}"#
        ),
        Solution::Determined(vec![Edit::ModifyBondConstraint {
            id: BondHandle::Id(BondId(1)),
            old: Some(BondConstraintAst::CisTransStereo(
                CisTransStereoAst::Stereo(StereoCoset::Lit(1)),
            )),
            new: None,
        }])
    )]
    #[case::bond_error(
        StereoInconsistencyPolicy::Error,
        mol_dsl_ground!(
            r#"{:atoms ["C #h3" "C #h2" "C #h1"] :bonds [[0 1 "1"] [1 2 "2#C1"]]}"#
        ),
        Solution::Contradictory(StereoContradiction::UnrealizableBond(BondId(1)))
    )]
    fn test_stereo_resolver_plan_inconsistency(
        stereo_model: StereoModel,
        #[case] policy: StereoInconsistencyPolicy,
        #[case] molecule: MoleculeAst,
        #[case] expected: Solution<Vec<Edit>, StereoContradiction>,
    ) {
        assert_eq!(
            StereoResolver::with_config(
                &stereo_model,
                StereoResolveConfig {
                    reset_stereo_constraints: false,
                    inconsistency: policy,
                },
            )
            .plan(&molecule),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::tetrahedral(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
                             :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#),
        mol_dsl_ground!(r#"{
            :atoms ["C #h3" "C #h1" "N #h2" "O #h1"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :type "Th1"}]
        }"#)
    )]
    #[case::cis_trans(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
                             :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]}"#),
        mol_dsl_ground!(r#"{
            :atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :type "Ct1"}]
        }"#)
    )]
    fn test_stereo_resolver_resolve(
        stereo_model: StereoModel,
        #[case] mut molecule: MoleculeAst,
        #[case] expected: MoleculeAst,
    ) {
        let resolver = StereoResolver::with_config(
            &stereo_model,
            StereoResolveConfig {
                reset_stereo_constraints: true,
                inconsistency: StereoInconsistencyPolicy::Error,
            },
        );
        assert_eq!(
            resolver.resolve(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::atom(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"]
                             :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        StereoContradiction::UnrealizableAtom(AtomId(1))
    )]
    #[case::bond(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h2" "C #h1"]
                             :bonds [[0 1 "1"] [1 2 "2#C1"]]}"#),
        StereoContradiction::UnrealizableBond(BondId(1))
    )]
    fn test_stereo_resolver_resolve_error(
        stereo_model: StereoModel,
        #[case] mut molecule: MoleculeAst,
        #[case] expected: StereoContradiction,
    ) {
        let original = molecule.clone();
        assert_eq!(
            StereoResolver::new(&stereo_model).resolve(&mut molecule),
            Ok(Solution::Contradictory(expected))
        );
        assert_eq!(molecule, original);
    }
}
