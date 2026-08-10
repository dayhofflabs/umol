//! Structural stereo resolver. Planning reads `#T` / `#C` assertions from the
//! materialized aromaticity state and emits stereo-element additions plus
//! optional source-constraint removals without mutating the molecule.

use std::collections::{BTreeMap, BTreeSet};

use thiserror::Error;
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomHandle, AtomUpdate, BondConstraintForm, BondHandle, BondUpdate,
    CisTransStereoForm, Edits, Lattice, Molecule, StereoAtomHandle, StereoBondHandle,
    TetrahedralStereoForm, TransactionError,
};
use umol_utils::solution::Solution;

use crate::ops::model::StereoModel;
use crate::ops::stereo::{StereoInconsistency, StereoPerception};

/// How stereo resolution handles an independently invalid constraint or entity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StereoFailurePolicy {
    Error,
    Keep,
    Remove,
}

/// How stereo resolution handles an independently valid constraint and entity that disagree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StereoMismatchPolicy {
    Error,
    Keep,
    RemoveConstraint,
    ReplaceEntity,
    RemoveBoth,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StereoResolveConfig {
    pub tetrahedral_stereo_failure: StereoFailurePolicy,
    pub stereo_atom_failure: StereoFailurePolicy,
    pub tetrahedral_stereo_mismatch: StereoMismatchPolicy,
    pub cis_trans_stereo_failure: StereoFailurePolicy,
    pub stereo_bond_failure: StereoFailurePolicy,
    pub cis_trans_stereo_mismatch: StereoMismatchPolicy,
    pub reset_stereo_constraints: bool,
}

impl Default for StereoResolveConfig {
    fn default() -> Self {
        Self {
            tetrahedral_stereo_failure: StereoFailurePolicy::Error,
            stereo_atom_failure: StereoFailurePolicy::Error,
            tetrahedral_stereo_mismatch: StereoMismatchPolicy::Error,
            cis_trans_stereo_failure: StereoFailurePolicy::Error,
            stereo_bond_failure: StereoFailurePolicy::Error,
            cis_trans_stereo_mismatch: StereoMismatchPolicy::Error,
            reset_stereo_constraints: false,
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
    #[error("stereo inconsistency: {0}")]
    Inconsistency(#[from] StereoInconsistency),
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
        ast: &Molecule,
    ) -> Result<Solution<Edits, StereoContradiction>, StereoError> {
        let partial_atom_constraint = ast.atoms().iter().any(|atom| {
            atom.constraints()
                .tetrahedral_stereo()
                .is_some_and(|constraint| !constraint.is_undetermined() && !constraint.is_ground())
        });
        let partial_bond_constraint = ast.bonds().iter().any(|bond| {
            bond.constraints()
                .cis_trans_stereo()
                .is_some_and(|constraint| !constraint.is_undetermined() && !constraint.is_ground())
        });
        if partial_atom_constraint || partial_bond_constraint {
            return Ok(Solution::Underdetermined(Edits::new()));
        }

        let derivation = self.perception.derive(ast);
        for &inconsistency in &derivation.inconsistencies {
            let error = match inconsistency {
                StereoInconsistency::TetrahedralStereoFailure { .. } => {
                    self.config.tetrahedral_stereo_failure == StereoFailurePolicy::Error
                }
                StereoInconsistency::StereoAtomFailure { .. } => {
                    self.config.stereo_atom_failure == StereoFailurePolicy::Error
                }
                StereoInconsistency::TetrahedralStereoMismatch { .. } => {
                    self.config.tetrahedral_stereo_mismatch == StereoMismatchPolicy::Error
                }
                StereoInconsistency::CisTransStereoFailure { .. } => {
                    self.config.cis_trans_stereo_failure == StereoFailurePolicy::Error
                }
                StereoInconsistency::StereoBondFailure { .. } => {
                    self.config.stereo_bond_failure == StereoFailurePolicy::Error
                }
                StereoInconsistency::CisTransStereoMismatch { .. } => {
                    self.config.cis_trans_stereo_mismatch == StereoMismatchPolicy::Error
                }
            };
            if error {
                return Ok(Solution::Contradictory(inconsistency.into()));
            }
        }

        let atoms: BTreeMap<_, _> = derivation
            .atoms
            .into_iter()
            .map(|(id, ligands, stereo)| (id, (ligands, stereo)))
            .collect();
        let bonds: BTreeMap<_, _> = derivation
            .bonds
            .into_iter()
            .map(|(id, ligands, stereo)| (id, (ligands, stereo)))
            .collect();

        let mut edits = Edits::new();
        let mut remove_atom_constraints = BTreeSet::new();
        let mut remove_bond_constraints = BTreeSet::new();
        let mut remove_stereo_atoms = BTreeSet::new();
        let mut remove_stereo_bonds = BTreeSet::new();
        let mut suppressed_atoms = BTreeSet::new();
        let mut suppressed_bonds = BTreeSet::new();

        for inconsistency in derivation.inconsistencies {
            match inconsistency {
                StereoInconsistency::TetrahedralStereoFailure { atom } => {
                    match self.config.tetrahedral_stereo_failure {
                        StereoFailurePolicy::Error => unreachable!(),
                        StereoFailurePolicy::Keep => {}
                        StereoFailurePolicy::Remove => {
                            remove_atom_constraints.insert(atom);
                        }
                    }
                }
                StereoInconsistency::StereoAtomFailure { stereo_atom } => {
                    let site = ast.stereo_atom(stereo_atom).site_id();
                    match self.config.stereo_atom_failure {
                        StereoFailurePolicy::Error => unreachable!(),
                        StereoFailurePolicy::Keep => {
                            suppressed_atoms.insert(site);
                        }
                        StereoFailurePolicy::Remove => {
                            remove_stereo_atoms.insert(stereo_atom);
                        }
                    }
                }
                StereoInconsistency::TetrahedralStereoMismatch { atom, stereo_atom } => {
                    match self.config.tetrahedral_stereo_mismatch {
                        StereoMismatchPolicy::Error => unreachable!(),
                        StereoMismatchPolicy::Keep => {
                            suppressed_atoms.insert(atom);
                        }
                        StereoMismatchPolicy::RemoveConstraint => {
                            remove_atom_constraints.insert(atom);
                            suppressed_atoms.insert(atom);
                        }
                        StereoMismatchPolicy::ReplaceEntity => {
                            remove_stereo_atoms.insert(stereo_atom);
                        }
                        StereoMismatchPolicy::RemoveBoth => {
                            remove_atom_constraints.insert(atom);
                            remove_stereo_atoms.insert(stereo_atom);
                            suppressed_atoms.insert(atom);
                        }
                    }
                }
                StereoInconsistency::CisTransStereoFailure { bond } => {
                    match self.config.cis_trans_stereo_failure {
                        StereoFailurePolicy::Error => unreachable!(),
                        StereoFailurePolicy::Keep => {}
                        StereoFailurePolicy::Remove => {
                            remove_bond_constraints.insert(bond);
                        }
                    }
                }
                StereoInconsistency::StereoBondFailure { stereo_bond } => {
                    let site = ast.stereo_bond(stereo_bond).site_id();
                    match self.config.stereo_bond_failure {
                        StereoFailurePolicy::Error => unreachable!(),
                        StereoFailurePolicy::Keep => {
                            suppressed_bonds.insert(site);
                        }
                        StereoFailurePolicy::Remove => {
                            remove_stereo_bonds.insert(stereo_bond);
                        }
                    }
                }
                StereoInconsistency::CisTransStereoMismatch { bond, stereo_bond } => {
                    match self.config.cis_trans_stereo_mismatch {
                        StereoMismatchPolicy::Error => unreachable!(),
                        StereoMismatchPolicy::Keep => {
                            suppressed_bonds.insert(bond);
                        }
                        StereoMismatchPolicy::RemoveConstraint => {
                            remove_bond_constraints.insert(bond);
                            suppressed_bonds.insert(bond);
                        }
                        StereoMismatchPolicy::ReplaceEntity => {
                            remove_stereo_bonds.insert(stereo_bond);
                        }
                        StereoMismatchPolicy::RemoveBoth => {
                            remove_bond_constraints.insert(bond);
                            remove_stereo_bonds.insert(stereo_bond);
                            suppressed_bonds.insert(bond);
                        }
                    }
                }
            }
        }

        if !remove_stereo_atoms.is_empty() {
            edits.remove_stereo_atoms(
                remove_stereo_atoms
                    .iter()
                    .map(|&id| {
                        let view = ast.stereo_atom(id);
                        (
                            StereoAtomHandle::Id(id),
                            AtomHandle::Id(view.site_id()),
                            view.ligands()
                                .map(|ligand| (AtomHandle::Id(ligand.atom_id()), ligand.kind()))
                                .collect(),
                            view.attributes.clone(),
                        )
                    })
                    .collect(),
            );
        }
        if !remove_stereo_bonds.is_empty() {
            edits.remove_stereo_bonds(
                remove_stereo_bonds
                    .iter()
                    .map(|&id| {
                        let view = ast.stereo_bond(id);
                        (
                            StereoBondHandle::Id(id),
                            BondHandle::Id(view.site_id()),
                            view.ligands()
                                .map(|ligand| (AtomHandle::Id(ligand.atom_id()), ligand.kind()))
                                .collect(),
                            view.attributes.clone(),
                        )
                    })
                    .collect(),
            );
        }

        let retained_atom_sites: BTreeSet<_> = ast
            .stereo_atoms()
            .iter()
            .filter(|view| !remove_stereo_atoms.contains(&view.id))
            .map(|view| view.site_id())
            .collect();
        let retained_bond_sites: BTreeSet<_> = ast
            .stereo_bonds()
            .iter()
            .filter(|view| !remove_stereo_bonds.contains(&view.id))
            .map(|view| view.site_id())
            .collect();

        for (id, (ligands, stereo)) in atoms {
            if suppressed_atoms.contains(&id) || retained_atom_sites.contains(&id) {
                continue;
            }
            edits.add_stereo_atom(
                AtomHandle::Id(id),
                ligands
                    .into_iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                stereo,
            );
            if self.config.reset_stereo_constraints {
                remove_atom_constraints.insert(id);
            }
        }
        for (id, (ligands, stereo)) in bonds {
            if suppressed_bonds.contains(&id) || retained_bond_sites.contains(&id) {
                continue;
            }
            edits.add_stereo_bond(
                BondHandle::Id(id),
                ligands
                    .into_iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                stereo,
            );
            if self.config.reset_stereo_constraints {
                remove_bond_constraints.insert(id);
            }
        }

        for atom in remove_atom_constraints {
            let mut update = AtomUpdate::default();
            update
                .constraints
                .set(AtomConstraintForm::TetrahedralStereo(
                    TetrahedralStereoForm::Undetermined,
                ));
            edits.update_atom(AtomHandle::Id(atom), ast.atom(atom).attributes, &update);
        }
        for bond in remove_bond_constraints {
            let mut update = BondUpdate::default();
            update.constraints.set(BondConstraintForm::CisTransStereo(
                CisTransStereoForm::Undetermined,
            ));
            edits.update_bond(BondHandle::Id(bond), ast.bond(bond).attributes, &update);
        }
        Ok(Solution::Determined(edits))
    }

    /// Plan and atomically apply structural stereo resolution.
    pub fn resolve(
        &self,
        ast: &mut Molecule,
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
    use umol_graph_ir::ir::{
        AtomId, BondId, Edit, Edits, StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId,
        StereoCoset, StereoKind, StereoLigandKind,
    };
    use umol_graph_ir::mol_dsl_ground;

    use super::*;

    #[fixture]
    fn stereo_model() -> StereoModel {
        StereoModel::default()
    }

    #[fixture]
    fn tetrahedral_entity_failure_molecule() -> Molecule {
        mol_dsl_ground!(
            r#"{
            :atoms ["C#h3" "C#h" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{:site 1 :ligands [0] :attrs "Th1"}]
        }"#
        )
    }

    #[fixture]
    fn cis_trans_entity_failure_molecule() -> Molecule {
        mol_dsl_ground!(
            r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0] :attrs "Ct1"}]
        }"#
        )
    }

    #[fixture]
    fn tetrahedral_mismatch_molecule() -> Molecule {
        mol_dsl_ground!(
            r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Th0"}]
        }"#
        )
    }

    #[fixture]
    fn cis_trans_mismatch_molecule() -> Molecule {
        mol_dsl_ground!(
            r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct0"}]
        }"#
        )
    }

    #[rstest]
    fn test_stereo_resolve_config_default() {
        assert_eq!(
            StereoResolveConfig::default(),
            StereoResolveConfig {
                reset_stereo_constraints: false,
                tetrahedral_stereo_failure: StereoFailurePolicy::Error,
                stereo_atom_failure: StereoFailurePolicy::Error,
                tetrahedral_stereo_mismatch: StereoMismatchPolicy::Error,
                cis_trans_stereo_failure: StereoFailurePolicy::Error,
                stereo_bond_failure: StereoFailurePolicy::Error,
                cis_trans_stereo_mismatch: StereoMismatchPolicy::Error,
            }
        );
    }

    #[rstest]
    #[case::tetrahedral(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
                             :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#),
        Edits::from_iter([Edit::AddStereoAtom {
            site: AtomHandle::Id(AtomId(1)),
            ligands: vec![
                (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(1)), StereoLigandKind::ImplicitHydrogen),
            ],
            attributes: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        }])
    )]
    #[case::cis_trans(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
                             :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]}"#),
        Edits::from_iter([Edit::AddStereoBond {
            site: BondHandle::Id(BondId(1)),
            ligands: vec![
                (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(1)), StereoLigandKind::ImplicitHydrogen),
                (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(2)), StereoLigandKind::ImplicitHydrogen),
            ],
            attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        }])
    )]
    fn test_stereo_resolver_plan(
        stereo_model: StereoModel,
        #[case] molecule: Molecule,
        #[case] expected: Edits,
    ) {
        assert_eq!(
            StereoResolver::new(&stereo_model).plan(&molecule),
            Ok(Solution::Determined(expected))
        );
    }

    #[rstest]
    #[case::tetrahedral(mol_dsl_ground!(r#"{
        :atoms ["C #h3" "C #h1 #T+" "N #h2" "O #h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
    }"#))]
    #[case::cis_trans(mol_dsl_ground!(r#"{
        :atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
        :bonds [[0 1 "1"] [1 2 "2#C+"] [2 3 "1"]]
    }"#))]
    fn test_stereo_resolver_plan_partial(stereo_model: StereoModel, #[case] molecule: Molecule) {
        assert_eq!(
            StereoResolver::new(&stereo_model).plan(&molecule),
            Ok(Solution::Underdetermined(Edits::new()))
        );
    }

    #[rstest]
    #[case::no_assertion(mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h3"] :bonds [[0 1 "1"]]}"#))]
    #[case::vacuous(mol_dsl_ground!(r#"{:atoms ["C #h4 #T*"]}"#))]
    #[case::existing_atom(mol_dsl_ground!(r#"{
        :atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Th1"}]
    }"#))]
    fn test_stereo_resolver_plan_identity(stereo_model: StereoModel, #[case] molecule: Molecule) {
        assert_eq!(
            StereoResolver::new(&stereo_model).plan(&molecule),
            Ok(Solution::Determined(Edits::new()))
        );
    }

    #[rstest]
    #[case::tetrahedral_keep(
        StereoFailurePolicy::Keep,
        mol_dsl_ground!(
            r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ),
        Solution::Determined(Edits::new())
    )]
    #[case::tetrahedral_remove(
        StereoFailurePolicy::Remove,
        mol_dsl_ground!(
            r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ),
        Solution::Determined(Edits::from_iter([Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(1)),
            old: Some(AtomConstraintForm::TetrahedralStereo(
                TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)),
            )),
            new: None,
        }]))
    )]
    #[case::tetrahedral_error(
        StereoFailurePolicy::Error,
        mol_dsl_ground!(
            r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ),
        Solution::Contradictory(StereoContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) }
        ))
    )]
    #[case::cis_trans_keep(
        StereoFailurePolicy::Keep,
        mol_dsl_ground!(
            r#"{:atoms ["C #h3" "C #h2" "C #h1"] :bonds [[0 1 "1"] [1 2 "2#C1"]]}"#
        ),
        Solution::Determined(Edits::new())
    )]
    #[case::cis_trans_remove(
        StereoFailurePolicy::Remove,
        mol_dsl_ground!(
            r#"{:atoms ["C #h3" "C #h2" "C #h1"] :bonds [[0 1 "1"] [1 2 "2#C1"]]}"#
        ),
        Solution::Determined(Edits::from_iter([Edit::ModifyBondConstraint {
            id: BondHandle::Id(BondId(1)),
            old: Some(BondConstraintForm::CisTransStereo(
                CisTransStereoForm::Stereo(StereoCoset::Lit(1)),
            )),
            new: None,
        }]))
    )]
    #[case::cis_trans_error(
        StereoFailurePolicy::Error,
        mol_dsl_ground!(
            r#"{:atoms ["C #h3" "C #h2" "C #h1"] :bonds [[0 1 "1"] [1 2 "2#C1"]]}"#
        ),
        Solution::Contradictory(StereoContradiction::Inconsistency(
            StereoInconsistency::CisTransStereoFailure { bond: BondId(1) }
        ))
    )]
    fn test_stereo_resolver_plan_constraint_failure(
        stereo_model: StereoModel,
        #[case] policy: StereoFailurePolicy,
        #[case] molecule: Molecule,
        #[case] expected: Solution<Edits, StereoContradiction>,
    ) {
        assert_eq!(
            StereoResolver::with_config(
                &stereo_model,
                StereoResolveConfig {
                    tetrahedral_stereo_failure: policy,
                    cis_trans_stereo_failure: policy,
                    ..StereoResolveConfig::default()
                },
            )
            .plan(&molecule),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::error(
        StereoFailurePolicy::Error,
        Solution::Contradictory(StereoContradiction::Inconsistency(
            StereoInconsistency::StereoAtomFailure {
                stereo_atom: StereoAtomId(0),
            }
        ))
    )]
    #[case::keep(StereoFailurePolicy::Keep, Solution::Determined(Edits::new()))]
    #[case::remove(
        StereoFailurePolicy::Remove,
        Solution::Determined(Edits::from_iter([Edit::RemoveStereoAtoms {
            removes: vec![(
                StereoAtomHandle::Id(StereoAtomId(0)),
                AtomHandle::Id(AtomId(1)),
                vec![(AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
        }]))
    )]
    fn test_stereo_resolver_plan_stereo_atom_failure(
        stereo_model: StereoModel,
        tetrahedral_entity_failure_molecule: Molecule,
        #[case] policy: StereoFailurePolicy,
        #[case] expected: Solution<Edits, StereoContradiction>,
    ) {
        assert_eq!(
            StereoResolver::with_config(
                &stereo_model,
                StereoResolveConfig {
                    stereo_atom_failure: policy,
                    ..StereoResolveConfig::default()
                },
            )
            .plan(&tetrahedral_entity_failure_molecule),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::error(
        StereoFailurePolicy::Error,
        Solution::Contradictory(StereoContradiction::Inconsistency(
            StereoInconsistency::StereoBondFailure {
                stereo_bond: StereoBondId(0),
            }
        ))
    )]
    #[case::keep(StereoFailurePolicy::Keep, Solution::Determined(Edits::new()))]
    #[case::remove(
        StereoFailurePolicy::Remove,
        Solution::Determined(Edits::from_iter([Edit::RemoveStereoBonds {
            removes: vec![(
                StereoBondHandle::Id(StereoBondId(0)),
                BondHandle::Id(BondId(1)),
                vec![(AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom)],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
        }]))
    )]
    fn test_stereo_resolver_plan_stereo_bond_failure(
        stereo_model: StereoModel,
        cis_trans_entity_failure_molecule: Molecule,
        #[case] policy: StereoFailurePolicy,
        #[case] expected: Solution<Edits, StereoContradiction>,
    ) {
        assert_eq!(
            StereoResolver::with_config(
                &stereo_model,
                StereoResolveConfig {
                    stereo_bond_failure: policy,
                    ..StereoResolveConfig::default()
                },
            )
            .plan(&cis_trans_entity_failure_molecule),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::error(
        StereoMismatchPolicy::Error,
        Solution::Contradictory(StereoContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoMismatch {
                atom: AtomId(1),
                stereo_atom: StereoAtomId(0),
            }
        ))
    )]
    #[case::keep(StereoMismatchPolicy::Keep, Solution::Determined(Edits::new()))]
    #[case::remove_constraint(
        StereoMismatchPolicy::RemoveConstraint,
        Solution::Determined(Edits::from_iter([Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(1)),
            old: Some(AtomConstraintForm::TetrahedralStereo(
                TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)),
            )),
            new: None,
        }]))
    )]
    #[case::replace_entity(
        StereoMismatchPolicy::ReplaceEntity,
        Solution::Determined(Edits::from_iter([
            Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(1)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(1)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                )],
            },
            Edit::AddStereoAtom {
                site: AtomHandle::Id(AtomId(1)),
                ligands: vec![
                    (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(1)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                ],
                attributes: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            },
        ]))
    )]
    #[case::remove_both(
        StereoMismatchPolicy::RemoveBoth,
        Solution::Determined(Edits::from_iter([
            Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(1)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(1)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                )],
            },
            Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(1)),
                old: Some(AtomConstraintForm::TetrahedralStereo(
                    TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)),
                )),
                new: None,
            },
        ]))
    )]
    fn test_stereo_resolver_plan_tetrahedral_mismatch(
        stereo_model: StereoModel,
        tetrahedral_mismatch_molecule: Molecule,
        #[case] policy: StereoMismatchPolicy,
        #[case] expected: Solution<Edits, StereoContradiction>,
    ) {
        assert_eq!(
            StereoResolver::with_config(
                &stereo_model,
                StereoResolveConfig {
                    tetrahedral_stereo_mismatch: policy,
                    ..StereoResolveConfig::default()
                },
            )
            .plan(&tetrahedral_mismatch_molecule),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::error(
        StereoMismatchPolicy::Error,
        Solution::Contradictory(StereoContradiction::Inconsistency(
            StereoInconsistency::CisTransStereoMismatch {
                bond: BondId(1),
                stereo_bond: StereoBondId(0),
            }
        ))
    )]
    #[case::keep(StereoMismatchPolicy::Keep, Solution::Determined(Edits::new()))]
    #[case::remove_constraint(
        StereoMismatchPolicy::RemoveConstraint,
        Solution::Determined(Edits::from_iter([Edit::ModifyBondConstraint {
            id: BondHandle::Id(BondId(1)),
            old: Some(BondConstraintForm::CisTransStereo(
                CisTransStereoForm::Stereo(StereoCoset::Lit(1)),
            )),
            new: None,
        }]))
    )]
    #[case::replace_entity(
        StereoMismatchPolicy::ReplaceEntity,
        Solution::Determined(Edits::from_iter([
            Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(1)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(1)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(2)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                )],
            },
            Edit::AddStereoBond {
                site: BondHandle::Id(BondId(1)),
                ligands: vec![
                    (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(1)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                    (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(2)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                ],
                attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            },
        ]))
    )]
    #[case::remove_both(
        StereoMismatchPolicy::RemoveBoth,
        Solution::Determined(Edits::from_iter([
            Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(1)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(1)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(2)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                )],
            },
            Edit::ModifyBondConstraint {
                id: BondHandle::Id(BondId(1)),
                old: Some(BondConstraintForm::CisTransStereo(
                    CisTransStereoForm::Stereo(StereoCoset::Lit(1)),
                )),
                new: None,
            },
        ]))
    )]
    fn test_stereo_resolver_plan_cis_trans_mismatch(
        stereo_model: StereoModel,
        cis_trans_mismatch_molecule: Molecule,
        #[case] policy: StereoMismatchPolicy,
        #[case] expected: Solution<Edits, StereoContradiction>,
    ) {
        assert_eq!(
            StereoResolver::with_config(
                &stereo_model,
                StereoResolveConfig {
                    cis_trans_stereo_mismatch: policy,
                    ..StereoResolveConfig::default()
                },
            )
            .plan(&cis_trans_mismatch_molecule),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::tetrahedral_not_stereo(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h#T!" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Th1"}]
        }"#),
        Edits::from_iter([Edit::RemoveStereoAtoms {
            removes: vec![(
                StereoAtomHandle::Id(StereoAtomId(0)),
                AtomHandle::Id(AtomId(1)),
                vec![
                    (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(1)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
        }])
    )]
    #[case::cis_trans_not_stereo(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C!"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]
        }"#),
        Edits::from_iter([Edit::RemoveStereoBonds {
            removes: vec![(
                StereoBondHandle::Id(StereoBondId(0)),
                BondHandle::Id(BondId(1)),
                vec![
                    (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(1)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                    (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(2)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
        }])
    )]
    fn test_stereo_resolver_plan_not_stereo_mismatch(
        stereo_model: StereoModel,
        #[case] molecule: Molecule,
        #[case] expected: Edits,
    ) {
        assert_eq!(
            StereoResolver::with_config(
                &stereo_model,
                StereoResolveConfig {
                    tetrahedral_stereo_mismatch: StereoMismatchPolicy::ReplaceEntity,
                    cis_trans_stereo_mismatch: StereoMismatchPolicy::ReplaceEntity,
                    ..StereoResolveConfig::default()
                },
            )
            .plan(&molecule),
            Ok(Solution::Determined(expected))
        );
    }

    #[rstest]
    #[case::tetrahedral(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
                             :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#),
        mol_dsl_ground!(r#"{
            :atoms ["C #h3" "C #h1" "N #h2" "O #h1"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Th1"}]
        }"#)
    )]
    #[case::cis_trans(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
                             :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]}"#),
        mol_dsl_ground!(r#"{
            :atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]
        }"#)
    )]
    fn test_stereo_resolver_resolve(
        stereo_model: StereoModel,
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        let resolver = StereoResolver::with_config(
            &stereo_model,
            StereoResolveConfig {
                reset_stereo_constraints: true,
                ..StereoResolveConfig::default()
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
        StereoContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) }
        )
    )]
    #[case::bond(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h2" "C #h1"]
                             :bonds [[0 1 "1"] [1 2 "2#C1"]]}"#),
        StereoContradiction::Inconsistency(
            StereoInconsistency::CisTransStereoFailure { bond: BondId(1) }
        )
    )]
    fn test_stereo_resolver_resolve_error(
        stereo_model: StereoModel,
        #[case] mut molecule: Molecule,
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
