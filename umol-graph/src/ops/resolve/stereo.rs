//! Structural stereo resolver. Planning reads `#T` / `#C` assertions from the
//! materialized aromaticity state and emits stereo-element additions plus
//! optional source-constraint removals without mutating the molecule.

use std::collections::{BTreeMap, BTreeSet};

use thiserror::Error;
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomHandle, AtomId, AtomUpdate, BondConstraintForm, BondHandle, BondId,
    BondUpdate, CisTransStereoForm, Edits, Lattice, Molecule, StereoAtomHandle, StereoAtomId,
    StereoBondHandle, StereoBondId, StereoCoset, StereoKind, TetrahedralStereoForm,
    TransactionError,
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

/// Failures to express stereo entities as fixed-frame atom and bond assertions.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoProjectError {
    #[error("stereo atom {stereo_atom:?} does not have tetrahedral kind")]
    UnsupportedStereoAtom { stereo_atom: StereoAtomId },
    #[error("stereo bond {stereo_bond:?} does not have cis-trans kind")]
    UnsupportedStereoBond { stereo_bond: StereoBondId },
    #[error("stereo atom {stereo_atom:?} cannot be expressed in the tetrahedral assertion frame")]
    StereoAtomFrame { stereo_atom: StereoAtomId },
    #[error("stereo bond {stereo_bond:?} cannot be expressed in the cis-trans assertion frame")]
    StereoBondFrame { stereo_bond: StereoBondId },
    #[error("atom {atom:?} has an assertion incompatible with its stereo configuration")]
    AtomAssertion { atom: AtomId },
    #[error("bond {bond:?} has an assertion incompatible with its stereo configuration")]
    BondAssertion { bond: BondId },
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

    /// Construct the complete stereo edit plan without mutating `molecule`.
    pub fn plan(
        &self,
        molecule: &Molecule,
    ) -> Result<Solution<Edits, StereoContradiction>, StereoError> {
        let partial_atom_constraint = molecule.atoms().iter().any(|atom| {
            atom.constraints()
                .tetrahedral_stereo()
                .is_some_and(|constraint| !constraint.is_undetermined() && !constraint.is_ground())
        });
        let skipped = self.perception.skipped_stereo_bonds(molecule);
        let partial_bond_constraint = molecule.bonds().iter().any(|bond| {
            !skipped.contains(&bond.id)
                && bond
                    .constraints()
                    .cis_trans_stereo()
                    .is_some_and(|constraint| {
                        !constraint.is_undetermined() && !constraint.is_ground()
                    })
        });
        if partial_atom_constraint || partial_bond_constraint {
            return Ok(Solution::Underdetermined(Edits::new()));
        }

        let derivation = self.perception.derive(molecule);
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

        let stereo_atoms: BTreeMap<_, _> = derivation
            .stereo_atoms
            .into_iter()
            .map(|(id, ligands, stereo)| (id, (ligands, stereo)))
            .collect();
        let stereo_bonds: BTreeMap<_, _> = derivation
            .stereo_bonds
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
                    let site = molecule.stereo_atom(stereo_atom).site_id();
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
                    let site = molecule.stereo_bond(stereo_bond).site_id();
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

        for site in derivation.skipped_stereo_bonds {
            remove_bond_constraints.insert(site);
        }

        if !remove_stereo_atoms.is_empty() {
            edits.remove_stereo_atoms(
                remove_stereo_atoms
                    .iter()
                    .map(|&id| {
                        let view = molecule.stereo_atom(id);
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
                        let view = molecule.stereo_bond(id);
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

        let retained_atom_sites: BTreeSet<_> = molecule
            .stereo_atoms()
            .iter()
            .filter(|view| !remove_stereo_atoms.contains(&view.id))
            .map(|view| view.site_id())
            .collect();
        let retained_bond_sites: BTreeSet<_> = molecule
            .stereo_bonds()
            .iter()
            .filter(|view| !remove_stereo_bonds.contains(&view.id))
            .map(|view| view.site_id())
            .collect();

        for (id, (ligands, stereo)) in stereo_atoms {
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
        for (id, (ligands, stereo)) in stereo_bonds {
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
            edits.update_atom(
                AtomHandle::Id(atom),
                molecule.atom(atom).attributes,
                &update,
            );
        }
        for bond in remove_bond_constraints {
            let mut update = BondUpdate::default();
            update.constraints.set(BondConstraintForm::CisTransStereo(
                CisTransStereoForm::Undetermined,
            ));
            edits.update_bond(
                BondHandle::Id(bond),
                molecule.bond(bond).attributes,
                &update,
            );
        }
        Ok(Solution::Determined(edits))
    }

    /// Plan and atomically apply structural stereo resolution.
    pub fn resolve(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), StereoContradiction>, StereoError> {
        let edits = match self.plan(molecule)? {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let mut editor = molecule.edit();
        editor.transact(edits)?;
        *molecule = editor.build();
        Ok(Solution::Determined(()))
    }

    /// Replaces tetrahedral and cis-trans entities with their fixed-frame #T and #C assertions.
    ///
    /// Uses the model's per-site reference frames and transports each stored configuration
    /// into that frame. Actual ligands, implicit hydrogens, and lone pairs remain distinct.
    /// Transports open configurations without choosing a configuration, combining them with
    /// existing assertions by meet. It does not invoke whole-molecule perception.
    /// Resolver failure, mismatch, and constraint-reset policies
    /// do not discard information during projection. Entity removal uses the editor's ordinary
    /// constraint-compaction semantics.
    ///
    /// # Semantic properties
    ///
    /// Projection is idempotent. Transporting an entity's configuration together with its
    /// ligand frame preserves the resulting assertion. Other atom and bond fields are unchanged.
    /// Success publishes all edits together; every error preserves the input exactly.
    ///
    /// # Errors
    ///
    /// Rejects unsupported or undetermined kinds, unavailable model reference frames, ligand
    /// frames that cannot transport to those references, and conflicting existing assertions.
    pub fn project(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), StereoContradiction>, StereoProjectError> {
        if !molecule.has_stereo_atoms() && !molecule.has_stereo_bonds() {
            return Ok(Solution::Determined(()));
        }
        let mut editor = molecule.edit();
        for stereo in molecule.stereo_atoms().iter() {
            let stereo_atom = stereo.id;
            if stereo.attributes.configuration.kind() != Some(StereoKind::Tetrahedral) {
                return Err(StereoProjectError::UnsupportedStereoAtom { stereo_atom });
            }
            let atom = stereo.site_id();
            let (ligands, _) = self
                .perception
                .derive_stereo_atom(molecule, atom, &StereoCoset::Undetermined)
                .ok_or(StereoProjectError::StereoAtomFrame { stereo_atom })?;
            let coset = stereo
                .coset_for(ligands)
                .ok_or(StereoProjectError::StereoAtomFrame { stereo_atom })?;
            let projected = TetrahedralStereoForm::Stereo(coset);
            let attributes = editor.atom_mut(atom).attributes;
            let assertion = match attributes.constraints.tetrahedral_stereo() {
                Some(existing) => existing
                    .meet(&projected)
                    .ok_or(StereoProjectError::AtomAssertion { atom })?,
                None => projected,
            };
            attributes
                .constraints
                .set(AtomConstraintForm::TetrahedralStereo(assertion));
        }
        for stereo in molecule.stereo_bonds().iter() {
            let stereo_bond = stereo.id;
            if stereo.attributes.configuration.kind() != Some(StereoKind::CisTrans) {
                return Err(StereoProjectError::UnsupportedStereoBond { stereo_bond });
            }
            let bond = stereo.site_id();
            let (ligands, _) = self
                .perception
                .derive_stereo_bond(molecule, bond, &StereoCoset::Undetermined)
                .ok_or(StereoProjectError::StereoBondFrame { stereo_bond })?;
            let coset = stereo
                .coset_for(ligands)
                .ok_or(StereoProjectError::StereoBondFrame { stereo_bond })?;
            let projected = CisTransStereoForm::Stereo(coset);
            let attributes = editor.bond_mut(bond).attributes;
            let assertion = match attributes.constraints.cis_trans_stereo() {
                Some(existing) => existing
                    .meet(&projected)
                    .ok_or(StereoProjectError::BondAssertion { bond })?,
                None => projected,
            };
            attributes
                .constraints
                .set(BondConstraintForm::CisTransStereo(assertion));
        }
        if molecule.has_stereo_atoms() {
            editor.remove_stereo_atoms(&molecule.stereo_atoms().ids().collect::<Vec<_>>());
        }
        if molecule.has_stereo_bonds() {
            editor.remove_stereo_bonds(&molecule.stereo_bonds().ids().collect::<Vec<_>>());
        }
        *molecule = editor.build();
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
    use umol_graph_ir::mol_dsl_concrete;
    use umol_io::smiles::SmilesIoConfig;

    use super::*;
    use crate::ingest::ingest_smiles_with;
    use crate::ops::model::{ChemistryModel, ValenceModel};
    use crate::ops::resolve::ResolveConfig;

    #[fixture]
    fn stereo_model() -> StereoModel {
        StereoModel::default()
    }

    #[fixture]
    fn tetrahedral_entity_failure_molecule() -> Molecule {
        mol_dsl_concrete!(
            r#"{
            :atoms ["C#h3" "C#h" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Sp1"}]
        }"#
        )
    }

    #[fixture]
    fn cis_trans_entity_failure_molecule() -> Molecule {
        mol_dsl_concrete!(
            r#"{
            :atoms ["C#h3" "C#h0" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]
        }"#
        )
    }

    #[fixture]
    fn tetrahedral_mismatch_molecule() -> Molecule {
        mol_dsl_concrete!(
            r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Th0"}]
        }"#
        )
    }

    #[fixture]
    fn cis_trans_mismatch_molecule() -> Molecule {
        mol_dsl_concrete!(
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
        mol_dsl_concrete!(r#"{:atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
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
        mol_dsl_concrete!(r#"{:atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
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
    #[case::cyclooctene_asserted(
        mol_dsl_concrete!(r#"{:atoms ["C #h2" "C #h1" "C #h1" "C #h2" "C #h2" "C #h2" "C #h2" "C #h2"]
                             :bonds [[0 7 "1"] [0 1 "1"] [1 2 "2#C1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 7 "1"]]}"#),
        Edits::from_iter([Edit::AddStereoBond {
            site: BondHandle::Id(BondId(2)),
            ligands: vec![
                (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(1)), StereoLigandKind::ImplicitHydrogen),
                (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(2)), StereoLigandKind::ImplicitHydrogen),
            ],
            attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        }])
    )]
    #[case::cyclohexene_asserted(
        mol_dsl_concrete!(r#"{:atoms ["C #h1" "C #h1" "C #h2" "C #h2" "C #h2" "C #h2"]
                             :bonds [[0 5 "1"] [0 1 "2#C1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]}"#),
        Edits::from_iter([Edit::ModifyBondConstraint {
            id: BondHandle::Id(BondId(1)),
            old: Some(BondConstraintForm::CisTransStereo(
                CisTransStereoForm::Stereo(StereoCoset::Lit(1)),
            )),
            new: None,
        }])
    )]
    #[case::cyclohexene_asserted_undetermined(
        mol_dsl_concrete!(r#"{:atoms ["C #h1" "C #h1" "C #h2" "C #h2" "C #h2" "C #h2"]
                             :bonds [[0 5 "1"] [0 1 "2#C+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]}"#),
        Edits::from_iter([Edit::ModifyBondConstraint {
            id: BondHandle::Id(BondId(1)),
            old: Some(BondConstraintForm::CisTransStereo(
                CisTransStereoForm::Stereo(StereoCoset::Undetermined),
            )),
            new: None,
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
    #[case::cyclohexene_realized_at_zero(
        0,
        mol_dsl_concrete!(r#"{:atoms ["C #h1" "C #h1" "C #h2" "C #h2" "C #h2" "C #h2"]
                             :bonds [[0 5 "1"] [0 1 "2#C1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]}"#),
        Edits::from_iter([Edit::AddStereoBond {
            site: BondHandle::Id(BondId(1)),
            ligands: vec![
                (AtomHandle::Id(AtomId(5)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(0)), StereoLigandKind::ImplicitHydrogen),
                (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(1)), StereoLigandKind::ImplicitHydrogen),
            ],
            attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        }])
    )]
    #[case::cyclooctene_skipped_at_nine(
        9,
        mol_dsl_concrete!(r#"{:atoms ["C #h2" "C #h1" "C #h1" "C #h2" "C #h2" "C #h2" "C #h2" "C #h2"]
                             :bonds [[0 7 "1"] [0 1 "1"] [1 2 "2#C1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 7 "1"]]}"#),
        Edits::from_iter([Edit::ModifyBondConstraint {
            id: BondHandle::Id(BondId(2)),
            old: Some(BondConstraintForm::CisTransStereo(
                CisTransStereoForm::Stereo(StereoCoset::Lit(1)),
            )),
            new: None,
        }])
    )]
    fn test_stereo_resolver_plan_stereo_bond_minimum_ring_size(
        #[case] stereo_bond_minimum_ring_size: u32,
        #[case] molecule: Molecule,
        #[case] expected: Edits,
    ) {
        let model = StereoModel {
            stereo_bond_minimum_ring_size,
            ..StereoModel::default()
        };
        assert_eq!(
            StereoResolver::new(&model).plan(&molecule),
            Ok(Solution::Determined(expected))
        );
    }

    #[rstest]
    #[case::tetrahedral(mol_dsl_concrete!(r#"{
        :atoms ["C #h3" "C #h1 #T+" "N #h2" "O #h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
    }"#))]
    #[case::cis_trans(mol_dsl_concrete!(r#"{
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
    #[case::no_assertion(mol_dsl_concrete!(r#"{:atoms ["C #h3" "C #h3"] :bonds [[0 1 "1"]]}"#))]
    #[case::vacuous(mol_dsl_concrete!(r#"{:atoms ["C #h4 #T*"]}"#))]
    #[case::existing_atom(mol_dsl_concrete!(r#"{
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
        mol_dsl_concrete!(
            r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ),
        Solution::Determined(Edits::new())
    )]
    #[case::tetrahedral_remove(
        StereoFailurePolicy::Remove,
        mol_dsl_concrete!(
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
        mol_dsl_concrete!(
            r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ),
        Solution::Contradictory(StereoContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) }
        ))
    )]
    #[case::cis_trans_keep(
        StereoFailurePolicy::Keep,
        mol_dsl_concrete!(
            r#"{:atoms ["C #h3" "C #h2" "C #h1"] :bonds [[0 1 "1"] [1 2 "2#C1"]]}"#
        ),
        Solution::Determined(Edits::new())
    )]
    #[case::cis_trans_remove(
        StereoFailurePolicy::Remove,
        mol_dsl_concrete!(
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
        mol_dsl_concrete!(
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
                vec![
                    (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(1)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                ],
                StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(1)),
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
        mol_dsl_concrete!(r#"{
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
        mol_dsl_concrete!(r#"{
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
        mol_dsl_concrete!(r#"{:atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
                             :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C #h3" "C #h1" "N #h2" "O #h1"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Th1"}]
        }"#)
    )]
    #[case::cis_trans(
        mol_dsl_concrete!(r#"{:atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
                             :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]}"#),
        mol_dsl_concrete!(r#"{
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
        mol_dsl_concrete!(r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"]
                             :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        StereoContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) }
        )
    )]
    #[case::bond(
        mol_dsl_concrete!(r#"{:atoms ["C #h3" "C #h2" "C #h1"]
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

    #[rstest]
    #[case::actual_ligands(
        mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#T1" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#))]
    #[case::implicit_hydrogen(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs "Th0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#))]
    #[case::explicit_hydrogen(
        mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [4 1 2 3] :attrs "Th0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#T1" "F" "Cl" "Br" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#))]
    #[case::lone_pair(
        mol_dsl_concrete!(r#"{:atoms ["N#n1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:lp 0] 1 2 3] :attrs "Th0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["N#n1#T1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#))]
    #[case::tetrahedral_open(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs "Th*"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T+" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#))]
    #[case::existing_assertion(
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T0" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs "Th*"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T0" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#))]
    #[case::bond_side_swap(
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [3 2 4 5] :attrs "Ct0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2#C1"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#))]
    #[case::bond_endpoint_swap(
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[1 0 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [4 5 2 3] :attrs "Ct0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2#C0"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#))]
    #[case::bond_hydrogens(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "C#h1" "F" "Cl"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:h 0] 2 3 [:h 1]] :attrs "Ct0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "C#h1" "F" "Cl"] :bonds [[0 1 "2#C1"] [0 2 "1"] [1 3 "1"]]}"#))]
    #[case::bond_lone_pair(
        mol_dsl_concrete!(r#"{:atoms ["N#n1" "C#h1" "C#h3" "F"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:lp 0] 2 3 [:h 1]] :attrs "Ct0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["N#n1" "C#h1" "C#h3" "F"] :bonds [[0 1 "2#C1"] [0 2 "1"] [1 3 "1"]]}"#))]
    #[case::bond_either(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "C#h1" "F" "Cl"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:h 0] 2 3 [:h 1]] :attrs "Ct*"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "C#h1" "F" "Cl"] :bonds [[0 1 "2#C+"] [0 2 "1"] [1 3 "1"]]}"#))]
    fn test_stereo_resolver_project(
        stereo_model: StereoModel,
        #[values(false, true)] permissive: bool,
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        let config = if permissive {
            StereoResolveConfig {
                tetrahedral_stereo_failure: StereoFailurePolicy::Remove,
                stereo_atom_failure: StereoFailurePolicy::Remove,
                tetrahedral_stereo_mismatch: StereoMismatchPolicy::RemoveBoth,
                cis_trans_stereo_failure: StereoFailurePolicy::Remove,
                stereo_bond_failure: StereoFailurePolicy::Remove,
                cis_trans_stereo_mismatch: StereoMismatchPolicy::RemoveBoth,
                reset_stereo_constraints: true,
            }
        } else {
            StereoResolveConfig::default()
        };
        let resolver = StereoResolver::with_config(&stereo_model, config);
        assert_eq!(
            resolver.project(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
        assert_eq!(
            resolver.project(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::empty(Molecule::new())]
    #[case::assertions(mol_dsl_concrete!(r#"{:atoms ["C#h1#T+" "C#h1#T!"] :bonds [[0 1 "2#C+"]]}"#))]
    fn test_stereo_resolver_project_identity(
        stereo_model: StereoModel,
        #[case] mut molecule: Molecule,
    ) {
        let original = molecule.clone();
        assert_eq!(
            StereoResolver::new(&stereo_model).project(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::atom_kind(
        mol_dsl_concrete!(r#"{:atoms ["Pt" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]}"#),
        StereoProjectError::UnsupportedStereoAtom { stereo_atom: StereoAtomId(0) })]
    #[case::bond_kind(
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ax0"}]}"#),
        StereoProjectError::UnsupportedStereoBond { stereo_bond: StereoBondId(0) })]
    #[case::atom_kind_undetermined(
        mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "*"}]}"#),
        StereoProjectError::UnsupportedStereoAtom { stereo_atom: StereoAtomId(0) })]
    #[case::bond_kind_undetermined(
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "*"}]}"#),
        StereoProjectError::UnsupportedStereoBond { stereo_bond: StereoBondId(0) })]
    #[case::atom_frame(
        mol_dsl_concrete!(r#"{:atoms ["C#h0#n0" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#),
        StereoProjectError::StereoAtomFrame { stereo_atom: StereoAtomId(0) })]
    #[case::virtual_kind_mismatch(
        mol_dsl_concrete!(r#"{:atoms ["N#h1#n1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:lp 0]] :attrs "Th0"}]}"#),
        StereoProjectError::StereoAtomFrame { stereo_atom: StereoAtomId(0) })]
    #[case::bond_frame(
        mol_dsl_concrete!(r#"{:atoms ["C#h0#n0" "C#h1" "F" "Cl"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 [:h 0] 3 [:h 1]] :attrs "Ct0"}]}"#),
        StereoProjectError::StereoBondFrame { stereo_bond: StereoBondId(0) })]
    #[case::atom_assertion(
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#),
        StereoProjectError::AtomAssertion { atom: AtomId(0) })]
    #[case::late_bond_assertion(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "F" "Cl" "Br" "C#h1" "C#h1" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [4 5 "2#C0"] [4 6 "1"] [5 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]
            :stereo-bonds [{:site 3 :ligands [6 [:h 4] 7 [:h 5]] :attrs "Ct1"}]}"#),
        StereoProjectError::BondAssertion { bond: BondId(3) })]
    fn test_stereo_resolver_project_error(
        stereo_model: StereoModel,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        policy: StereoFailurePolicy,
        #[values(
            StereoMismatchPolicy::Error,
            StereoMismatchPolicy::Keep,
            StereoMismatchPolicy::RemoveBoth
        )]
        mismatch: StereoMismatchPolicy,
        #[case] mut molecule: Molecule,
        #[case] expected: StereoProjectError,
    ) {
        let resolver = StereoResolver::with_config(
            &stereo_model,
            StereoResolveConfig {
                tetrahedral_stereo_failure: policy,
                stereo_atom_failure: policy,
                cis_trans_stereo_failure: policy,
                stereo_bond_failure: policy,
                tetrahedral_stereo_mismatch: mismatch,
                cis_trans_stereo_mismatch: mismatch,
                reset_stereo_constraints: true,
            },
        );
        let original = molecule.clone();
        assert_eq!(resolver.project(&mut molecule), Err(expected));
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::atom(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#), StereoKind::Tetrahedral,
        StereoProjectError::StereoAtomFrame { stereo_atom: StereoAtomId(0) })]
    #[case::bond(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "C#h1" "F" "Cl"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 [:h 0] 3 [:h 1]] :attrs "Ct0"}]}"#), StereoKind::CisTrans,
        StereoProjectError::StereoBondFrame { stereo_bond: StereoBondId(0) })]
    fn test_stereo_resolver_project_model(
        #[case] mut molecule: Molecule,
        #[case] kind: StereoKind,
        #[case] expected: StereoProjectError,
    ) {
        let mut model = StereoModel::default();
        model.kind_models[kind as usize] = None;
        let original = molecule.clone();
        assert_eq!(
            StereoResolver::new(&model).project(&mut molecule),
            Err(expected)
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::tetrahedral("[NH2][C@H](F)Cl",
        mol_dsl_concrete!(r#"{:atoms ["N#h2#n1" "C#h1#T0" "F#n3" "Cl#n3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#))]
    #[case::root_hydrogen("[C@H](F)(Cl)Br",
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T1" "F#n3" "Cl#n3" "Br#n3"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#))]
    #[case::explicit_hydrogen("[H][C@](F)(Cl)Br",
        mol_dsl_concrete!(r#"{:atoms ["H" "C#T0" "F#n3" "Cl#n3" "Br#n3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"] [1 4 "1"]]}"#))]
    #[case::lone_pair("[CH3][S@](=[O])[CH2][CH3]",
        mol_dsl_concrete!(r#"{:atoms ["C#h3" "S#n1#T0" "O#n2" "C#h2" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2"] [1 3 "1"] [3 4 "1"]]}"#))]
    #[case::cis_trans("[CH3]/[CH]=[CH]/[CH3]",
        mol_dsl_concrete!(r#"{:atoms ["C#h3" "C#h1" "C#h1" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]}"#))]
    fn test_stereo_resolver_project_input(
        #[values(ValenceModel::smiles(), ValenceModel::default())] valence: ValenceModel,
        #[case] input: &str,
        #[case] expected: Molecule,
    ) {
        let model = ChemistryModel {
            valence,
            ..Default::default()
        };
        let mut molecule = ingest_smiles_with(
            input,
            &SmilesIoConfig::opensmiles(),
            &model,
            &ResolveConfig::default(),
        )
        .unwrap();
        assert_eq!(
            StereoResolver::new(&model.stereo).project(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }
}
