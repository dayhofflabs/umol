//! Structural stereo resolver. Planning reads complete `#T` / `#C` evidence from the
//! materialized aromaticity state and emits stereo-element additions, virtual-ligand
//! identification, and the discharge of source assertions without mutating the molecule.

use std::collections::{BTreeMap, BTreeSet};

use thiserror::Error;
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomConstraintKey, AtomHandle, AtomId, BondConstraintForm,
    BondConstraintKey, BondHandle, BondId, CisTransStereoForm, Edit, Edits, Lattice, Molecule,
    Normalize, RingConfig, StereoAtomHandle, StereoAtomId, StereoBondHandle, StereoBondId,
    StereoKind, TetrahedralStereoForm,
};
use umol_utils::solution::Solution;

use crate::ops::model::StereoModel;
use crate::ops::stereo::{StereoContradiction, StereoInconsistency, StereoPerceiver};

/// How stereo resolution handles an unrealizable constraint or entity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StereoFailurePolicy {
    Error,
    Keep,
    Remove,
}

/// Ring-detection algorithms and recovery policies for structural stereo resolution.
///
/// Ambiguity remains underdetermined under every policy. Mismatches are always contradictory.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StereoResolveConfig {
    /// Policy for unrealizable tetrahedral constraints and stored stereo atoms.
    pub stereo_atom_failure: StereoFailurePolicy,
    /// Policy for unrealizable cis/trans constraints and stored stereo bonds.
    pub stereo_bond_failure: StereoFailurePolicy,
    pub ring_config: RingConfig,
}

impl Default for StereoResolveConfig {
    fn default() -> Self {
        Self {
            stereo_atom_failure: StereoFailurePolicy::Error,
            stereo_bond_failure: StereoFailurePolicy::Error,
            ring_config: RingConfig::default(),
        }
    }
}

/// Resolves structural stereo from the molecule's evidence.
///
/// # Semantic properties
///
/// Plans apply to their input and preserve molecule integrity. Consuming and in-place
/// operations agree on outcomes and determined molecules for the same input and configuration.
#[derive(Clone, Debug)]
pub struct StereoResolver {
    perception: StereoPerceiver,
    config: StereoResolveConfig,
}

/// Failures to express stereo entities as atom and bond kind constraints.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoProjectError {
    #[error("stereo atom {stereo_atom:?} does not have tetrahedral kind")]
    UnsupportedStereoAtom { stereo_atom: StereoAtomId },
    #[error("stereo bond {stereo_bond:?} does not have cis-trans kind")]
    UnsupportedStereoBond { stereo_bond: StereoBondId },
    #[error("atom {atom:?} has an assertion incompatible with its stereo kind")]
    AtomAssertion { atom: AtomId },
    #[error("bond {bond:?} has an assertion incompatible with its stereo kind")]
    BondAssertion { bond: BondId },
}

impl StereoResolver {
    pub fn new(model: &StereoModel) -> Self {
        Self::with_config(model, StereoResolveConfig::default())
    }

    pub fn with_config(model: &StereoModel, config: StereoResolveConfig) -> Self {
        Self {
            perception: StereoPerceiver::new(model),
            config,
        }
    }

    /// Construct the complete stereo edit plan without mutating `molecule`.
    /// New and retained entities discharge independently entailed stereo constraints.
    /// New entities have undetermined cosets; retained configurations are preserved.
    /// Discharge preserves other constraint keys.
    /// Requests below the creation cutoff are retained under every failure policy.
    /// Identified virtual ligand kinds replace stored kinds in place, preserving frame order
    /// and all entity attributes. Keep retains a failed value; Remove clears its assertion or entity.
    /// Mismatched evidence is contradictory regardless of the failure policies.
    /// Established findings are checked before underdetermination. An underdetermined
    /// derivation produces an empty plan and publishes no partial changes.
    pub fn plan_resolve(&self, molecule: &Molecule) -> Solution<Edits, StereoContradiction> {
        let (derivation, underdetermined) =
            match self.perception.derive(molecule, self.config.ring_config) {
                Solution::Determined(derivation) => (derivation, false),
                Solution::Underdetermined(derivation) => (derivation, true),
                Solution::Contradictory(contradiction) => {
                    return Solution::Contradictory(contradiction);
                }
            };
        for &inconsistency in &derivation.inconsistencies {
            let error = match inconsistency {
                StereoInconsistency::TetrahedralStereoFailure { .. } => {
                    self.config.stereo_atom_failure == StereoFailurePolicy::Error
                }
                StereoInconsistency::StereoAtomFailure { .. } => {
                    self.config.stereo_atom_failure == StereoFailurePolicy::Error
                }
                StereoInconsistency::CisTransStereoFailure { .. } => {
                    self.config.stereo_bond_failure == StereoFailurePolicy::Error
                }
                StereoInconsistency::StereoBondFailure { .. } => {
                    self.config.stereo_bond_failure == StereoFailurePolicy::Error
                }
                StereoInconsistency::TetrahedralStereoMismatch { .. }
                | StereoInconsistency::CisTransStereoMismatch { .. } => true,
            };
            if error {
                return Solution::Contradictory(inconsistency.into());
            }
        }

        if underdetermined {
            return Solution::Underdetermined(Edits::new());
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
                    match self.config.stereo_atom_failure {
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
                StereoInconsistency::CisTransStereoFailure { bond } => {
                    match self.config.stereo_bond_failure {
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
                StereoInconsistency::TetrahedralStereoMismatch { .. }
                | StereoInconsistency::CisTransStereoMismatch { .. } => {
                    unreachable!("mismatches were rejected before planning")
                }
            }
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
                            view.attributes().clone(),
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
                            view.attributes().clone(),
                        )
                    })
                    .collect(),
            );
        }

        for (id, (ligands, stereo)) in stereo_atoms {
            if suppressed_atoms.contains(&id) {
                continue;
            }
            if let Some(existing) = molecule
                .stereo_atoms()
                .at(id)
                .filter(|view| !remove_stereo_atoms.contains(&view.id()))
            {
                if existing.ligand_ids() != ligands {
                    edits.push(Edit::ReplaceStereoAtomLigands {
                        id: StereoAtomHandle::Id(existing.id()),
                        old: existing
                            .ligand_ids()
                            .iter()
                            .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                            .collect(),
                        new: ligands
                            .into_iter()
                            .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                            .collect(),
                    });
                }
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
            remove_atom_constraints.insert(id);
        }
        for (id, (ligands, stereo)) in stereo_bonds {
            if suppressed_bonds.contains(&id) {
                continue;
            }
            if let Some(existing) = molecule
                .stereo_bonds()
                .at(id)
                .filter(|view| !remove_stereo_bonds.contains(&view.id()))
            {
                if existing.ligand_ids() != ligands {
                    edits.push(Edit::ReplaceStereoBondLigands {
                        id: StereoBondHandle::Id(existing.id()),
                        old: existing
                            .ligand_ids()
                            .iter()
                            .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                            .collect(),
                        new: ligands
                            .into_iter()
                            .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                            .collect(),
                    });
                }
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
            remove_bond_constraints.insert(id);
        }

        for atom in molecule.atoms().iter() {
            if molecule
                .stereo_atoms()
                .at(atom.id())
                .is_some_and(|stereo| !remove_stereo_atoms.contains(&stereo.id()))
                && atom
                    .attributes()
                    .constraints
                    .get(AtomConstraintKey::TetrahedralStereo)
                    .is_some_and(|asserted| {
                        atom.constraints()
                            .derived(AtomConstraintKey::TetrahedralStereo)
                            .is_some_and(|derived| asserted.matches(&derived))
                    })
            {
                remove_atom_constraints.insert(atom.id());
            }
        }
        for bond in molecule.bonds().iter() {
            if molecule
                .stereo_bonds()
                .at(bond.id())
                .is_some_and(|stereo| !remove_stereo_bonds.contains(&stereo.id()))
                && bond
                    .attributes()
                    .constraints
                    .get(BondConstraintKey::CisTransStereo)
                    .is_some_and(|asserted| {
                        bond.constraints()
                            .derived(BondConstraintKey::CisTransStereo)
                            .is_some_and(|derived| asserted.matches(&derived))
                    })
            {
                remove_bond_constraints.insert(bond.id());
            }
        }

        let mut atom_constraints = Vec::new();
        for atom in remove_atom_constraints {
            if let Some(asserted) = molecule
                .atom(atom)
                .attributes()
                .constraints
                .get(AtomConstraintKey::TetrahedralStereo)
            {
                atom_constraints.push((AtomHandle::Id(atom), Some(asserted.clone()), None));
            }
        }
        if !atom_constraints.is_empty() {
            edits.modify_atom_constraints(AtomConstraintKey::TetrahedralStereo, atom_constraints);
        }
        let mut bond_constraints = Vec::new();
        for bond in remove_bond_constraints {
            if let Some(asserted) = molecule
                .bond(bond)
                .attributes()
                .constraints
                .get(BondConstraintKey::CisTransStereo)
            {
                bond_constraints.push((BondHandle::Id(bond), Some(asserted.clone()), None));
            }
        }
        if !bond_constraints.is_empty() {
            edits.modify_bond_constraints(BondConstraintKey::CisTransStereo, bond_constraints);
        }
        Solution::Determined(edits)
    }

    /// Resolves the molecule, consuming its storage without a recovery journal.
    ///
    /// Uses the same plan as [StereoResolver::resolve_into]. A determined result has checked
    /// molecule integrity; its stereo configurations may remain open.
    /// Underdetermined and Contradictory drop the input.
    pub fn resolve(&self, molecule: Molecule) -> Solution<Molecule, StereoContradiction, ()> {
        let edits = match self.plan_resolve(&molecule) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Solution::Underdetermined(()),
            Solution::Contradictory(contradiction) => {
                return Solution::Contradictory(contradiction);
            }
        };
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("stereo resolve plan applies to its input");
        let molecule = editor
            .finish()
            .expect("stereo resolve plan preserves molecule integrity");
        Solution::Determined(molecule)
    }

    /// Plans and applies structural stereo resolution without a recovery journal.
    ///
    /// # Semantic properties
    ///
    /// Planning finishes before mutation. Only Determined changes the molecule;
    /// every other outcome leaves it unchanged.
    /// With consistent, realizable stereo evidence, resolution is idempotent under `==`.
    /// Stereo configurations may remain open in a determined result.
    pub fn resolve_into(&self, molecule: &mut Molecule) -> Solution<(), StereoContradiction> {
        let edits = match self.plan_resolve(molecule) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Solution::Underdetermined(()),
            Solution::Contradictory(contradiction) => {
                return Solution::Contradictory(contradiction);
            }
        };
        if !edits.is_empty() {
            molecule.mutate(edits);
        }
        Solution::Determined(())
    }

    /// Replaces supported stereo entities with kind constraints, discarding their configurations.
    ///
    /// # Errors
    ///
    /// Rejects unsupported kinds and conflicting scalar constraints before mutation.
    pub fn project(molecule: Molecule) -> Result<Molecule, StereoProjectError> {
        let edits = Self::plan_project(&molecule)?;
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("stereo project plan applies to its input");
        Ok(editor
            .finish()
            .expect("stereo project plan preserves molecule integrity"))
    }

    /// Replaces supported stereo entities with positive scalar constraints in place.
    ///
    /// Ligand frames, configurations, and entity-local constraints are discarded.
    /// Other atom and bond data are preserved, with ordinary reference compaction.
    ///
    /// # Semantic properties
    ///
    /// Projection is idempotent and independent of ligand order and configuration.
    /// Every rejection leaves the molecule unchanged.
    ///
    /// # Errors
    ///
    /// Rejects unsupported or undetermined kinds and conflicting scalar constraints.
    pub fn project_into(molecule: &mut Molecule) -> Result<(), StereoProjectError> {
        let edits = Self::plan_project(molecule)?;
        if !edits.is_empty() {
            molecule.mutate(edits);
        }
        Ok(())
    }

    /// Plans replacement of supported stereo entities with positive scalar constraints.
    ///
    /// The plan discards frames, configurations, and entity-local constraints without
    /// interpreting them. Applying it to the unchanged input preserves integrity.
    ///
    /// # Errors
    ///
    /// Rejects unsupported or undetermined kinds and conflicting scalar constraints.
    pub fn plan_project(molecule: &Molecule) -> Result<Edits, StereoProjectError> {
        let mut atom_constraints = Vec::new();
        for stereo in molecule.stereo_atoms().iter() {
            if stereo.attributes().configuration.kind() != Some(StereoKind::Tetrahedral) {
                return Err(StereoProjectError::UnsupportedStereoAtom {
                    stereo_atom: stereo.id(),
                });
            }
            let atom = stereo.site_id();
            let attributes = molecule.atom(atom).attributes();
            let assertion = attributes
                .constraints
                .tetrahedral_stereo()
                .meet(&TetrahedralStereoForm::Stereo)
                .ok_or(StereoProjectError::AtomAssertion { atom })?;
            let new = AtomConstraintForm::TetrahedralStereo(assertion);
            let old = attributes
                .constraints
                .get(AtomConstraintKey::TetrahedralStereo);
            if old.is_none_or(|old| !old.normalized_eq(&new)) {
                atom_constraints.push((AtomHandle::Id(atom), old.cloned(), Some(new)));
            }
        }
        let mut bond_constraints = Vec::new();
        for stereo in molecule.stereo_bonds().iter() {
            if stereo.attributes().configuration.kind() != Some(StereoKind::CisTrans) {
                return Err(StereoProjectError::UnsupportedStereoBond {
                    stereo_bond: stereo.id(),
                });
            }
            let bond = stereo.site_id();
            let attributes = molecule.bond(bond).attributes();
            let assertion = attributes
                .constraints
                .cis_trans_stereo()
                .meet(&CisTransStereoForm::Stereo)
                .ok_or(StereoProjectError::BondAssertion { bond })?;
            let new = BondConstraintForm::CisTransStereo(assertion);
            let old = attributes
                .constraints
                .get(BondConstraintKey::CisTransStereo);
            if old.is_none_or(|old| !old.normalized_eq(&new)) {
                bond_constraints.push((BondHandle::Id(bond), old.cloned(), Some(new)));
            }
        }
        let mut edits = Edits::new();
        if !atom_constraints.is_empty() {
            edits.modify_atom_constraints(AtomConstraintKey::TetrahedralStereo, atom_constraints);
        }
        if !bond_constraints.is_empty() {
            edits.modify_bond_constraints(BondConstraintKey::CisTransStereo, bond_constraints);
        }
        if molecule.has_stereo_atoms() {
            edits.remove_stereo_atoms(
                molecule
                    .stereo_atoms()
                    .iter()
                    .map(|stereo| {
                        (
                            StereoAtomHandle::Id(stereo.id()),
                            AtomHandle::Id(stereo.site_id()),
                            stereo
                                .ligands()
                                .map(|ligand| (AtomHandle::Id(ligand.atom_id()), ligand.kind()))
                                .collect(),
                            stereo.attributes().clone(),
                        )
                    })
                    .collect(),
            );
        }
        if molecule.has_stereo_bonds() {
            edits.remove_stereo_bonds(
                molecule
                    .stereo_bonds()
                    .iter()
                    .map(|stereo| {
                        (
                            StereoBondHandle::Id(stereo.id()),
                            BondHandle::Id(stereo.site_id()),
                            stereo
                                .ligands()
                                .map(|ligand| (AtomHandle::Id(ligand.atom_id()), ligand.kind()))
                                .collect(),
                            stereo.attributes().clone(),
                        )
                    })
                    .collect(),
            );
        }
        Ok(edits)
    }
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{
        AddStereoAtom, AddStereoBond, AtomForm, AtomId, BondForm, BondId, Edit, Edits,
        MoleculeEntries, NumForm, StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId,
        StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    };
    use umol_graph_ir::{mol_dsl, mol_dsl_concrete};
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
            :atoms ["C#h3" "C#h#T!" "N#h2" "O#h"]
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
            :bonds [[0 1 "1"] [1 2 "2#C!"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct0"}]
        }"#
        )
    }

    #[fixture]
    fn underdetermined_molecule() -> Molecule {
        let mut molecule = mol_dsl_concrete!(
            r#"{
            :atoms ["C" "F" "Cl" "Br" "I" "C#h1" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
                    [5 6 "1"] [5 7 "1"] [5 8 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}
                           {:site 5 :ligands [[:virt 5] 8 7 6] :attrs "Th0"}]
        }"#
        );
        molecule
            .atom_mut(AtomId(0))
            .attributes_mut()
            .implicit_hydrogens = NumForm::var("hydrogens");
        molecule
    }

    #[rstest]
    fn test_stereo_resolve_config_default() {
        assert_eq!(
            StereoResolveConfig::default(),
            StereoResolveConfig {
                stereo_atom_failure: StereoFailurePolicy::Error,
                stereo_bond_failure: StereoFailurePolicy::Error,
                ring_config: Default::default(),
            }
        );
    }

    #[rstest]
    #[case::tetrahedral(
        mol_dsl_concrete!(
            r#"{:atoms ["C #h3" "C #h1 #T+" "N #h2" "O #h1"]
                             :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#
        ),
        Edits::from_iter([
            Edit::AddStereoAtom(Box::new(AddStereoAtom {
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
                attributes: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
            })),
            Edit::ModifyAtomConstraints {
                key: AtomConstraintKey::TetrahedralStereo,
                changes: vec![(
                    AtomHandle::Id(AtomId(1)),
                    Some(AtomConstraintForm::TetrahedralStereo(
                        TetrahedralStereoForm::Stereo,
                    )),
                    None,
                )],
            },
        ]),
    )]
    #[case::cis_trans(
        mol_dsl_concrete!(
            r#"{:atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
                             :bonds [[0 1 "1"] [1 2 "2#C+"] [2 3 "1"]]}"#
        ),
        Edits::from_iter([
            Edit::AddStereoBond(Box::new(AddStereoBond {
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
                attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
            })),
            Edit::ModifyBondConstraints {
                key: BondConstraintKey::CisTransStereo,
                changes: vec![(
                    BondHandle::Id(BondId(1)),
                    Some(BondConstraintForm::CisTransStereo(
                        CisTransStereoForm::Stereo,
                    )),
                    None,
                )],
            },
        ]),
    )]
    #[case::cyclooctene_asserted(
        mol_dsl_concrete!(
            r#"{:atoms ["C #h2" "C #h1" "C #h1" "C #h2" "C #h2" "C #h2" "C #h2" "C #h2"]
                             :bonds [[0 7 "1"] [0 1 "1"] [1 2 "2#C+"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 7 "1"]]}"#
        ),
        Edits::from_iter([
            Edit::AddStereoBond(Box::new(AddStereoBond {
                site: BondHandle::Id(BondId(2)),
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
                attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
            })),
            Edit::ModifyBondConstraints {
                key: BondConstraintKey::CisTransStereo,
                changes: vec![(
                    BondHandle::Id(BondId(2)),
                    Some(BondConstraintForm::CisTransStereo(
                        CisTransStereoForm::Stereo,
                    )),
                    None,
                )],
            },
        ]),
    )]
    #[case::cyclohexene_asserted(
        mol_dsl_concrete!(
            r#"{:atoms ["C #h1" "C #h1" "C #h2" "C #h2" "C #h2" "C #h2"]
                             :bonds [[0 5 "1"] [0 1 "2#C+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]}"#
        ),
        Edits::new(),
    )]
    #[case::cyclohexene_asserted_undetermined(
        mol_dsl_concrete!(
            r#"{:atoms ["C #h1" "C #h1" "C #h2" "C #h2" "C #h2" "C #h2"]
                             :bonds [[0 5 "1"] [0 1 "2#C+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]}"#
        ),
        Edits::new(),
    )]
    fn test_stereo_resolver_plan_resolve(
        stereo_model: StereoModel,
        #[case] molecule: Molecule,
        #[case] expected: Edits,
    ) {
        assert_eq!(
            StereoResolver::new(&stereo_model).plan_resolve(&molecule),
            Solution::Determined(expected)
        );
    }

    #[rstest]
    #[case::cyclohexene_realized_at_zero(
        0,
        mol_dsl_concrete!(
            r#"{:atoms ["C #h1" "C #h1" "C #h2" "C #h2" "C #h2" "C #h2"]
                             :bonds [[0 5 "1"] [0 1 "2#C+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]}"#
        ),
        Edits::from_iter([
            Edit::AddStereoBond(Box::new(AddStereoBond {
                site: BondHandle::Id(BondId(1)),
                ligands: vec![
                    (AtomHandle::Id(AtomId(5)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(0)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                    (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(1)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                ],
                attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
            })),
            Edit::ModifyBondConstraints {
                key: BondConstraintKey::CisTransStereo,
                changes: vec![(
                    BondHandle::Id(BondId(1)),
                    Some(BondConstraintForm::CisTransStereo(
                        CisTransStereoForm::Stereo,
                    )),
                    None,
                )],
            },
        ]),
    )]
    #[case::cyclooctene_skipped_at_nine(
        9,
        mol_dsl_concrete!(
            r#"{:atoms ["C #h2" "C #h1" "C #h1" "C #h2" "C #h2" "C #h2" "C #h2" "C #h2"]
                             :bonds [[0 7 "1"] [0 1 "1"] [1 2 "2#C+"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 7 "1"]]}"#
        ),
        Edits::new(),
    )]
    fn test_stereo_resolver_plan_resolve_ring(
        #[case] stereo_bond_min_ring_size: u32,
        #[case] molecule: Molecule,
        #[case] expected: Edits,
    ) {
        let resolver = StereoResolver::with_config(
            &StereoModel {
                stereo_bond_min_ring_size,
                ..Default::default()
            },
            StereoResolveConfig::default(),
        );
        assert_eq!(
            resolver.plan_resolve(&molecule),
            Solution::Determined(expected)
        );
    }

    #[rstest]
    #[case::boundary(Some(8), 8)]
    #[case::disabled(Some(6), 0)]
    #[case::large_boundary(Some(24), 24)]
    #[case::acyclic(None, 64)]
    fn test_stereo_resolver_resolve_ring(
        #[case] ring_size: Option<u32>,
        #[case] stereo_bond_min_ring_size: u32,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        failure: StereoFailurePolicy,
    ) {
        let size = ring_size.unwrap_or(4);
        let mut entries = MoleculeEntries {
            atoms: (0..size)
                .map(|id| {
                    AtomForm::from_element(Element::C)
                        .with_implicit_hydrogens(if id < 2 {
                            1
                        } else if ring_size.is_some() {
                            2
                        } else {
                            3
                        })
                        .into_concrete()
                })
                .collect(),
            bonds: if ring_size.is_some() {
                (0..size)
                    .map(|id| {
                        (
                            AtomId(id),
                            AtomId((id + 1) % size),
                            BondForm::from_order(if id == 0 { 2 } else { 1 }).into_concrete(),
                        )
                    })
                    .collect()
            } else {
                vec![
                    (
                        AtomId(0),
                        AtomId(1),
                        BondForm::from_order(2).into_concrete(),
                    ),
                    (
                        AtomId(0),
                        AtomId(2),
                        BondForm::from_order(1).into_concrete(),
                    ),
                    (
                        AtomId(1),
                        AtomId(3),
                        BondForm::from_order(1).into_concrete(),
                    ),
                ]
            },
            ..Default::default()
        };
        let unasserted_bond = entries.bonds[0].2.clone();
        entries.bonds[0].2 =
            unasserted_bond
                .clone()
                .with_constraint(BondConstraintForm::cis_trans_stereo(
                    CisTransStereoForm::Stereo,
                ));
        let molecule = Molecule::from_entries(entries.clone());
        entries.bonds[0].2 = unasserted_bond;
        entries.stereo_bonds = vec![(
            BondId(0),
            vec![
                StereoLigand::new(
                    AtomId(if ring_size.is_some() { size - 1 } else { 2 }),
                    StereoLigandKind::Atom,
                ),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(
                    AtomId(if ring_size.is_some() { 2 } else { 3 }),
                    StereoLigandKind::Atom,
                ),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
        )];
        let expected = Molecule::from_entries(entries);
        let resolver = StereoResolver::with_config(
            &StereoModel {
                stereo_bond_min_ring_size,
                ..Default::default()
            },
            StereoResolveConfig {
                stereo_atom_failure: failure,
                stereo_bond_failure: failure,
                ring_config: RingConfig::default(),
            },
        );
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Solution::Determined(expected.clone())
        );
        let mut in_place = molecule;
        assert_eq!(
            resolver.resolve_into(&mut in_place),
            Solution::Determined(())
        );
        assert_eq!(in_place, expected);
    }

    #[rstest]
    #[case::small(6, 8)]
    #[case::above(8, 9)]
    #[case::large(24, 25)]
    fn test_stereo_resolver_resolve_ring_identity(
        #[case] size: u32,
        #[case] stereo_bond_min_ring_size: u32,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        failure: StereoFailurePolicy,
    ) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: (0..size)
                .map(|id| {
                    AtomForm::from_element(Element::C)
                        .with_implicit_hydrogens(if id < 2 { 1 } else { 2 })
                        .into_concrete()
                })
                .collect(),
            bonds: (0..size)
                .map(|id| {
                    let form = if id == 0 {
                        BondForm::from_order(2).with_constraint(
                            BondConstraintForm::cis_trans_stereo(CisTransStereoForm::Stereo),
                        )
                    } else {
                        BondForm::from_order(1)
                    };
                    (AtomId(id), AtomId((id + 1) % size), form.into_concrete())
                })
                .collect(),
            ..Default::default()
        });
        let resolver = StereoResolver::with_config(
            &StereoModel {
                stereo_bond_min_ring_size,
                ..Default::default()
            },
            StereoResolveConfig {
                stereo_atom_failure: failure,
                stereo_bond_failure: failure,
                ring_config: RingConfig::default(),
            },
        );
        assert_eq!(
            resolver.plan_resolve(&molecule),
            Solution::Determined(Edits::new())
        );
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Solution::Determined(molecule.clone())
        );
        let mut in_place = molecule.clone();
        assert_eq!(
            resolver.resolve_into(&mut in_place),
            Solution::Determined(())
        );
        assert_eq!(in_place, molecule);
    }

    #[rstest]
    #[case::no_assertion(mol_dsl_concrete!(r#"{:atoms ["C #h3" "C #h3"] :bonds [[0 1 "1"]]}"#))]
    #[case::vacuous(mol_dsl_concrete!(r#"{:atoms ["C #h4 #T*"]}"#))]
    #[case::vacuous_bond(mol_dsl_concrete!(r#"{
        :atoms ["C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "2#C*"] [0 2 "1"] [1 3 "1"]]
    }"#))]
    #[case::existing_atom(mol_dsl_concrete!(r#"{
        :atoms ["C #h3" "C #h1 " "N #h2" "O #h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Th1"}]
    }"#))]
    #[case::refined_atom(mol_dsl_concrete!(r#"{
        :atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]
    }"#))]
    #[case::refined_bond(mol_dsl_concrete!(r#"{
        :atoms ["C" "C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct1"}]
    }"#))]
    #[case::reframed_atom(mol_dsl_concrete!(r#"{
        :atoms ["C#h" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs "Th0"}]
    }"#))]
    #[case::reframed_bond(mol_dsl_concrete!(r#"{
        :atoms ["C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
        :stereo-bonds [{:site 0 :ligands [[:h 1] 3 2 [:h 0]] :attrs "Ct0"}]
    }"#))]
    #[case::open_atom(mol_dsl_concrete!(r#"{
        :atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]
    }"#))]
    #[case::open_bond(mol_dsl_concrete!(r#"{
        :atoms ["C" "C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct*"}]
    }"#))]
    #[case::small_ring_bond(mol_dsl_concrete!(r#"{
        :atoms ["C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h2"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :stereo-bonds [{:site 0 :ligands [[:h 1] 2 5 [:h 0]] :attrs "Ct1#g/#o(0,1)/#p(0,1)!#f(0,1)"}]
    }"#))]
    fn test_stereo_resolver_plan_resolve_identity(
        stereo_model: StereoModel,
        #[case] molecule: Molecule,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        failure: StereoFailurePolicy,
        #[values(0, 8, 64)] stereo_bond_min_ring_size: u32,
    ) {
        let resolver = StereoResolver::with_config(
            &StereoModel {
                stereo_bond_min_ring_size,
                ..stereo_model
            },
            StereoResolveConfig {
                stereo_atom_failure: failure,
                stereo_bond_failure: failure,
                ring_config: RingConfig::default(),
            },
        );
        assert_eq!(
            resolver.plan_resolve(&molecule),
            Solution::Determined(Edits::new())
        );
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Solution::Determined(molecule)
        );
    }

    #[rstest]
    #[case::hydrogen(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 1 "2#C+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:virt 1] 2 5 [:virt 0]] :attrs "Ct1#g/#o(0,1)/#p(0,1)!#f(0,1)"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:h 1] 2 5 [:h 0]] :attrs "Ct1#g/#o(0,1)/#p(0,1)!#f(0,1)"}]
        }"#),
    )]
    #[case::lone_pair(
        mol_dsl_concrete!(r#"{
            :atoms ["N#n" "N#n" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 1 "2#C+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:virt 1] 2 5 [:virt 0]] :attrs "Ct0"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["N#n" "N#n" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:lp 1] 2 5 [:lp 0]] :attrs "Ct0"}]
        }"#),
    )]
    fn test_stereo_resolver_resolve_ring_virtual(
        #[case] molecule: Molecule,
        #[case] expected: Molecule,
        #[values(0, 8, 64)] stereo_bond_min_ring_size: u32,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        failure: StereoFailurePolicy,
    ) {
        let resolver = StereoResolver::with_config(
            &StereoModel {
                stereo_bond_min_ring_size,
                ..Default::default()
            },
            StereoResolveConfig {
                stereo_atom_failure: failure,
                stereo_bond_failure: failure,
                ring_config: RingConfig::default(),
            },
        );
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Solution::Determined(expected.clone())
        );
        let mut in_place = molecule;
        assert_eq!(
            resolver.resolve_into(&mut in_place),
            Solution::Determined(())
        );
        assert_eq!(in_place, expected);
    }

    #[rstest]
    #[case::tetrahedral_keep(
        StereoResolveConfig {
            stereo_atom_failure: StereoFailurePolicy::Keep,
            stereo_bond_failure: StereoFailurePolicy::Error,
            ring_config: Default::default(),
        },
        mol_dsl_concrete!(
            r#"{:atoms ["C #h3" "S #h0 #T+" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ),
        Solution::Determined(Edits::new())
    )]
    #[case::tetrahedral_remove(
        StereoResolveConfig {
            stereo_atom_failure: StereoFailurePolicy::Remove,
            stereo_bond_failure: StereoFailurePolicy::Error,
            ring_config: Default::default(),
        },
        mol_dsl_concrete!(r#"{:atoms ["C #h3" "S #h0 #T+" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        Solution::Determined(Edits::from_iter([Edit::ModifyAtomConstraints {
            key: AtomConstraintKey::TetrahedralStereo,
            changes: vec![(
                AtomHandle::Id(AtomId(1)),
                Some(AtomConstraintForm::TetrahedralStereo(
                    TetrahedralStereoForm::Stereo,
                )),
                None,
            )],
        }])),
    )]
    #[case::tetrahedral_error(
        StereoResolveConfig {
            stereo_atom_failure: StereoFailurePolicy::Error,
            stereo_bond_failure: StereoFailurePolicy::Keep,
            ring_config: Default::default(),
        },
        mol_dsl_concrete!(
            r#"{:atoms ["C #h3" "S #h0 #T+" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ),
        Solution::Contradictory(StereoContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) }
        ))
    )]
    #[case::cis_trans_keep(
        StereoResolveConfig {
            stereo_atom_failure: StereoFailurePolicy::Error,
            stereo_bond_failure: StereoFailurePolicy::Keep,
            ring_config: Default::default(),
        },
        mol_dsl_concrete!(
            r#"{:atoms ["C #h3" "C #h2" "C #h1"] :bonds [[0 1 "1"] [1 2 "2#C+"]]}"#
        ),
        Solution::Determined(Edits::new())
    )]
    #[case::cis_trans_remove(
        StereoResolveConfig {
            stereo_atom_failure: StereoFailurePolicy::Error,
            stereo_bond_failure: StereoFailurePolicy::Remove,
            ring_config: Default::default(),
        },
        mol_dsl_concrete!(r#"{:atoms ["C #h3" "C #h2" "C #h1"] :bonds [[0 1 "1"] [1 2 "2#C+"]]}"#),
        Solution::Determined(Edits::from_iter([Edit::ModifyBondConstraints {
            key: BondConstraintKey::CisTransStereo,
            changes: vec![(
                BondHandle::Id(BondId(1)),
                Some(BondConstraintForm::CisTransStereo(
                    CisTransStereoForm::Stereo,
                )),
                None,
            )],
        }])),
    )]
    #[case::cis_trans_error(
        StereoResolveConfig {
            stereo_atom_failure: StereoFailurePolicy::Keep,
            stereo_bond_failure: StereoFailurePolicy::Error,
            ring_config: Default::default(),
        },
        mol_dsl_concrete!(
            r#"{:atoms ["C #h3" "C #h2" "C #h1"] :bonds [[0 1 "1"] [1 2 "2#C+"]]}"#
        ),
        Solution::Contradictory(StereoContradiction::Inconsistency(
            StereoInconsistency::CisTransStereoFailure { bond: BondId(1) }
        ))
    )]
    fn test_stereo_resolver_plan_resolve_constraint_failure(
        stereo_model: StereoModel,
        #[case] config: StereoResolveConfig,
        #[case] molecule: Molecule,
        #[case] expected: Solution<Edits, StereoContradiction>,
    ) {
        assert_eq!(
            StereoResolver::with_config(&stereo_model, config).plan_resolve(&molecule),
            expected
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
    fn test_stereo_resolver_plan_resolve_stereo_atom_failure(
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
            .plan_resolve(&tetrahedral_entity_failure_molecule),
            expected
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
    fn test_stereo_resolver_plan_resolve_stereo_bond_failure(
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
            .plan_resolve(&cis_trans_entity_failure_molecule),
            expected
        );
    }

    #[rstest]
    fn test_stereo_resolver_plan_resolve_tetrahedral_mismatch(
        stereo_model: StereoModel,
        mut tetrahedral_mismatch_molecule: Molecule,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        failure: StereoFailurePolicy,
    ) {
        let original = tetrahedral_mismatch_molecule.clone();
        let resolver = StereoResolver::with_config(
            &stereo_model,
            StereoResolveConfig {
                stereo_atom_failure: failure,
                stereo_bond_failure: failure,
                ring_config: Default::default(),
            },
        );
        let expected =
            StereoContradiction::Inconsistency(StereoInconsistency::TetrahedralStereoMismatch {
                atom: AtomId(1),
                stereo_atom: StereoAtomId(0),
            });
        assert_eq!(
            resolver.plan_resolve(&tetrahedral_mismatch_molecule),
            Solution::Contradictory(expected.clone())
        );
        assert_eq!(
            resolver.resolve(tetrahedral_mismatch_molecule.clone()),
            Solution::Contradictory(expected.clone())
        );
        assert_eq!(
            resolver.resolve_into(&mut tetrahedral_mismatch_molecule),
            Solution::Contradictory(expected)
        );
        assert_eq!(tetrahedral_mismatch_molecule, original);
    }

    #[rstest]
    fn test_stereo_resolver_plan_resolve_cis_trans_mismatch(
        stereo_model: StereoModel,
        mut cis_trans_mismatch_molecule: Molecule,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        failure: StereoFailurePolicy,
    ) {
        let original = cis_trans_mismatch_molecule.clone();
        let resolver = StereoResolver::with_config(
            &stereo_model,
            StereoResolveConfig {
                stereo_atom_failure: failure,
                stereo_bond_failure: failure,
                ring_config: Default::default(),
            },
        );
        let expected =
            StereoContradiction::Inconsistency(StereoInconsistency::CisTransStereoMismatch {
                bond: BondId(1),
                stereo_bond: StereoBondId(0),
            });
        assert_eq!(
            resolver.plan_resolve(&cis_trans_mismatch_molecule),
            Solution::Contradictory(expected.clone())
        );
        assert_eq!(
            resolver.resolve(cis_trans_mismatch_molecule.clone()),
            Solution::Contradictory(expected.clone())
        );
        assert_eq!(
            resolver.resolve_into(&mut cis_trans_mismatch_molecule),
            Solution::Contradictory(expected)
        );
        assert_eq!(cis_trans_mismatch_molecule, original);
    }

    #[rstest]
    #[case::tetrahedral(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T!" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Th1"}]
        }"#),
        StereoInconsistency::TetrahedralStereoMismatch { atom: AtomId(1), stereo_atom: StereoAtomId(0) },
    )]
    #[case::cis_trans(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C!"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]
        }"#),
        StereoInconsistency::CisTransStereoMismatch { bond: BondId(1), stereo_bond: StereoBondId(0) },
    )]
    fn test_stereo_resolver_plan_resolve_not_stereo_mismatch(
        stereo_model: StereoModel,
        #[case] mut molecule: Molecule,
        #[case] expected: StereoInconsistency,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        failure: StereoFailurePolicy,
    ) {
        let original = molecule.clone();
        let resolver = StereoResolver::with_config(
            &stereo_model,
            StereoResolveConfig {
                stereo_atom_failure: failure,
                stereo_bond_failure: failure,
                ring_config: Default::default(),
            },
        );
        assert_eq!(
            resolver.plan_resolve(&molecule),
            Solution::Contradictory(expected.into())
        );
        assert_eq!(
            resolver.resolve_into(&mut molecule),
            Solution::Contradictory(expected.into())
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::tetrahedral_open(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T+" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Th*"}]
        }"#)
    )]
    #[case::cis_trans_open(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C+"] [2 3 "1"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct*"}]
        }"#)
    )]
    #[case::tetrahedral_hydrogen_and_lone_pair(
        mol_dsl_concrete!(r#"{
            :atoms ["P#h#n#T+" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["P#h#n" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 [:h 0] [:lp 0]] :attrs "Th*"}]
        }"#)
    )]
    #[case::cis_trans_hydrogen_and_lone_pair(
        mol_dsl_concrete!(r#"{
            :atoms ["N#h#n" "N#h#n"] :bonds [[0 1 "2#C+"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["N#h#n" "N#h#n"] :bonds [[0 1 "2"]]
            :stereo-bonds [{:site 0 :ligands [[:h 0] [:lp 0] [:h 1] [:lp 1]] :attrs "Ct*"}]
        }"#)
    )]
    #[case::mixed_stored_and_added(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#T+#v3" "F" "Cl" "Br" "C#h#T+" "F" "Cl" "Br"
                    "C#h" "C#h" "F" "Cl" "C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]
                    [4 5 "1"] [4 6 "1"] [4 7 "1"]
                    [8 9 "2#C+#R!"] [8 10 "1"] [9 11 "1"]
                    [12 13 "2#C+"] [12 14 "1"] [13 15 "1"]]
            :stereo-atoms [{:site 4 :ligands [5 6 7 [:h 4]] :attrs "Th1"}]
            :stereo-bonds [{:site 9 :ligands [14 [:h 12] 15 [:h 13]] :attrs "Ct1"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#v3" "F" "Cl" "Br" "C#h" "F" "Cl" "Br"
                    "C#h" "C#h" "F" "Cl" "C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]
                    [4 5 "1"] [4 6 "1"] [4 7 "1"]
                    [8 9 "2#R!"] [8 10 "1"] [9 11 "1"]
                    [12 13 "2"] [12 14 "1"] [13 15 "1"]]
            :stereo-atoms [{:site 4 :ligands [5 6 7 [:h 4]] :attrs "Th1"}
                           {:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th*"}]
            :stereo-bonds [{:site 9 :ligands [14 [:h 12] 15 [:h 13]] :attrs "Ct1"}
                           {:site 6 :ligands [10 [:h 8] 11 [:h 9]] :attrs "Ct*"}]
        }"#)
    )]
    #[case::virtual_tetrahedral_hydrogen(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:virt 0] 3 2 1] :attrs "Th1#g/#o(0,1)/#p(0,1)!#f(0,1)"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:h 0] 3 2 1] :attrs "Th1#g/#o(0,1)/#p(0,1)!#f(0,1)"}]
        }"#),
    )]
    #[case::virtual_tetrahedral_lone_pair(
        mol_dsl_concrete!(r#"{
            :atoms ["N#n" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [2 [:virt 0] 1 3] :attrs "Th*"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["N#n" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [2 [:lp 0] 1 3] :attrs "Th*"}]
        }"#),
    )]
    #[case::virtual_with_fixed_lone_pair(
        mol_dsl_concrete!(r#"{
            :atoms ["P#h#n#T+" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 [:lp 0] [:virt 0]] :attrs "Th0"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["P#h#n" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 [:lp 0] [:h 0]] :attrs "Th0"}]
        }"#),
    )]
    #[case::virtual_cis_trans_hydrogens(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 [:virt 0] [:virt 1] 3] :attrs "Ct0#g/#o(0,2)/"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 [:h 0] [:h 1] 3] :attrs "Ct0#g/#o(0,2)/"}]
        }"#),
    )]
    #[case::virtual_cis_trans_lone_pairs(
        mol_dsl_concrete!(r#"{
            :atoms ["N#n" "N#n" "C#h3" "C#h3"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [3 [:virt 1] [:virt 0] 2] :attrs "Ct1"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["N#n" "N#n" "C#h3" "C#h3"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [3 [:lp 1] [:lp 0] 2] :attrs "Ct1"}]
        }"#),
    )]
    #[case::virtual_cis_trans_mixed_kinds(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "N#n" "F" "Cl"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:virt 1] 3 2 [:virt 0]] :attrs "Ct*"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "N#n" "F" "Cl"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:lp 1] 3 2 [:h 0]] :attrs "Ct*"}]
        }"#),
    )]
    fn test_stereo_resolver_resolve(
        stereo_model: StereoModel,
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        let resolver = StereoResolver::new(&stereo_model);
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Solution::Determined(expected.clone())
        );
        assert_eq!(
            resolver.resolve_into(&mut molecule),
            Solution::Determined(())
        );
        assert_eq!(molecule, expected);
        assert_eq!(resolver.resolve(molecule), Solution::Determined(expected));
    }

    #[rstest]
    #[case::atom_assertion(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#T+" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        }"#),
        StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(0) },
    )]
    #[case::atom_entity(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        }"#),
        StereoInconsistency::StereoAtomFailure { stereo_atom: StereoAtomId(0) },
    )]
    #[case::bond_assertion(
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C" "F" "Cl" "Br" "I" "H"]
            :bonds [[0 1 "2#C+"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"] [1 6 "1"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C" "F" "Cl" "Br" "I" "H"]
            :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"] [1 6 "1"]]
        }"#),
        StereoInconsistency::CisTransStereoFailure { bond: BondId(0) },
    )]
    #[case::bond_entity(
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C" "F" "Cl" "Br" "I" "H"]
            :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"] [1 6 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct1"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C" "F" "Cl" "Br" "I" "H"]
            :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"] [1 6 "1"]]
        }"#),
        StereoInconsistency::StereoBondFailure { stereo_bond: StereoBondId(0) },
    )]
    fn test_stereo_resolver_resolve_ligands(
        stereo_model: StereoModel,
        #[case] mut molecule: Molecule,
        #[case] removed: Molecule,
        #[case] inconsistency: StereoInconsistency,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        policy: StereoFailurePolicy,
    ) {
        let original = molecule.clone();
        let resolver = StereoResolver::with_config(
            &stereo_model,
            StereoResolveConfig {
                stereo_atom_failure: policy,
                stereo_bond_failure: policy,
                ring_config: Default::default(),
            },
        );
        let expected = if policy == StereoFailurePolicy::Remove {
            removed
        } else {
            original.clone()
        };
        let (consuming, in_place) = if policy == StereoFailurePolicy::Error {
            (
                Solution::Contradictory(inconsistency.into()),
                Solution::Contradictory(inconsistency.into()),
            )
        } else {
            (
                Solution::Determined(expected.clone()),
                Solution::Determined(()),
            )
        };
        assert_eq!(resolver.resolve(original), consuming);
        assert_eq!(resolver.resolve_into(&mut molecule), in_place);
        assert_eq!(molecule, expected);
    }

    #[rstest]
    fn test_stereo_resolver_resolve_underdetermined(
        stereo_model: StereoModel,
        mut underdetermined_molecule: Molecule,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        stereo_atom_failure: StereoFailurePolicy,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        stereo_bond_failure: StereoFailurePolicy,
    ) {
        let original = underdetermined_molecule.clone();
        let resolver = StereoResolver::with_config(
            &stereo_model,
            StereoResolveConfig {
                stereo_atom_failure,
                stereo_bond_failure,
                ring_config: Default::default(),
            },
        );
        assert_eq!(
            resolver.plan_resolve(&underdetermined_molecule),
            Solution::Underdetermined(Edits::new())
        );
        assert_eq!(
            resolver.resolve(underdetermined_molecule.clone()),
            Solution::Underdetermined(())
        );
        assert_eq!(
            resolver.resolve_into(&mut underdetermined_molecule),
            Solution::Underdetermined(())
        );
        assert_eq!(underdetermined_molecule, original);
    }

    #[rstest]
    #[case::atom_failure(
        mol_dsl_concrete!(r#"{
            :atoms ["C" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:virt 0]] :attrs "Th1"}]
        }"#),
        StereoInconsistency::StereoAtomFailure { stereo_atom: StereoAtomId(2) },
        true,
    )]
    #[case::bond_failure(
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C#h1" "F" "Cl"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 [:virt 0] 3 [:virt 1]] :attrs "Ct1"}]
        }"#),
        StereoInconsistency::StereoBondFailure { stereo_bond: StereoBondId(0) },
        true,
    )]
    #[case::tetrahedral_mismatch(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h1#T!" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:virt 0]] :attrs "Th1"}]
        }"#),
        StereoInconsistency::TetrahedralStereoMismatch {
            atom: AtomId(9), stereo_atom: StereoAtomId(2),
        },
        false,
    )]
    #[case::cis_trans_mismatch(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h1" "C#h1" "F" "Cl"]
            :bonds [[0 1 "2#C!"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 [:virt 0] 3 [:virt 1]] :attrs "Ct1"}]
        }"#),
        StereoInconsistency::CisTransStereoMismatch {
            bond: BondId(7), stereo_bond: StereoBondId(0),
        },
        false,
    )]
    fn test_stereo_resolver_resolve_partial(
        stereo_model: StereoModel,
        underdetermined_molecule: Molecule,
        #[case] known: Molecule,
        #[case] inconsistency: StereoInconsistency,
        #[case] recoverable: bool,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        policy: StereoFailurePolicy,
    ) {
        let mut molecule = underdetermined_molecule.combine(&known);
        let original = molecule.clone();
        let resolver = StereoResolver::with_config(
            &stereo_model,
            StereoResolveConfig {
                stereo_atom_failure: policy,
                stereo_bond_failure: policy,
                ring_config: Default::default(),
            },
        );
        if recoverable && policy != StereoFailurePolicy::Error {
            assert_eq!(
                resolver.plan_resolve(&molecule),
                Solution::Underdetermined(Edits::new())
            );
            assert_eq!(
                resolver.resolve(molecule.clone()),
                Solution::Underdetermined(())
            );
            assert_eq!(
                resolver.resolve_into(&mut molecule),
                Solution::Underdetermined(())
            );
        } else {
            let expected = StereoContradiction::Inconsistency(inconsistency);
            assert_eq!(
                resolver.plan_resolve(&molecule),
                Solution::Contradictory(expected.clone())
            );
            assert_eq!(
                resolver.resolve(molecule.clone()),
                Solution::Contradictory(expected.clone())
            );
            assert_eq!(
                resolver.resolve_into(&mut molecule),
                Solution::Contradictory(expected)
            );
        }
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::atom(
        mol_dsl!(r#"{:atoms ["C#T+" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
        None
    )]
    #[case::bond(
        mol_dsl!(r#"{:atoms ["C" "C#h#n0" "F" "Cl"]
            :bonds [[0 1 "2#C+"] [0 2 "1"] [1 3 "1"]]}"#),
        None
    )]
    #[case::mixed_failure(
        mol_dsl!(r#"{:atoms ["C#T+" "F" "Cl" "Br" "C#h0#n0#T+"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
        Some(StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(4) })
    )]
    fn test_stereo_resolver_resolve_partial_requests(
        stereo_model: StereoModel,
        #[case] mut molecule: Molecule,
        #[case] failure: Option<StereoInconsistency>,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        policy: StereoFailurePolicy,
    ) {
        let original = molecule.clone();
        let resolver = StereoResolver::with_config(
            &stereo_model,
            StereoResolveConfig {
                stereo_atom_failure: policy,
                stereo_bond_failure: policy,
                ring_config: Default::default(),
            },
        );
        if let Some(failure) = failure.filter(|_| policy == StereoFailurePolicy::Error) {
            let contradiction = StereoContradiction::Inconsistency(failure);
            assert_eq!(
                resolver.plan_resolve(&molecule),
                Solution::Contradictory(contradiction.clone())
            );
            assert_eq!(
                resolver.resolve(molecule.clone()),
                Solution::Contradictory(contradiction.clone())
            );
            assert_eq!(
                resolver.resolve_into(&mut molecule),
                Solution::Contradictory(contradiction)
            );
        } else {
            assert_eq!(
                resolver.plan_resolve(&molecule),
                Solution::Underdetermined(Edits::new())
            );
            assert_eq!(
                resolver.resolve(molecule.clone()),
                Solution::Underdetermined(())
            );
            assert_eq!(
                resolver.resolve_into(&mut molecule),
                Solution::Underdetermined(())
            );
        }
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::hydrogen(
        NumForm::Lit(1),
        NumForm::Lit(0),
        Solution::Determined(StereoLigandKind::ImplicitHydrogen)
    )]
    #[case::lone_pair(
        NumForm::Lit(0),
        NumForm::Lit(1),
        Solution::Determined(StereoLigandKind::LonePair)
    )]
    #[case::absent(NumForm::Lit(0), NumForm::Lit(0), Solution::Contradictory(()))]
    #[case::both_present(NumForm::Lit(1), NumForm::Lit(1), Solution::Contradictory(()))]
    #[case::two_hydrogens(NumForm::Lit(2), NumForm::Lit(0), Solution::Contradictory(()))]
    #[case::two_lone_pairs(NumForm::Lit(0), NumForm::Lit(2), Solution::Contradictory(()))]
    #[case::open_counts(NumForm::Undetermined, NumForm::Undetermined, Solution::Underdetermined(()))]
    #[case::finite_counts(NumForm::lit_set([0, 1]), NumForm::lit_set([0, 1]), Solution::Underdetermined(()))]
    #[case::excluded_lone_pair(
        NumForm::Undetermined,
        NumForm::Lit(0),
        Solution::Determined(StereoLigandKind::ImplicitHydrogen)
    )]
    #[case::excluded_hydrogen(
        NumForm::Lit(0),
        NumForm::Undetermined,
        Solution::Determined(StereoLigandKind::LonePair)
    )]
    #[case::finite_hydrogen(NumForm::lit_set([0, 1]), NumForm::Lit(0), Solution::Determined(StereoLigandKind::ImplicitHydrogen))]
    #[case::no_supported_count(NumForm::range_from(2), NumForm::Lit(0), Solution::Contradictory(()))]
    fn test_stereo_resolver_resolve_virtual_ligands(
        stereo_model: StereoModel,
        #[case] hydrogens: NumForm,
        #[case] lone_pairs: NumForm,
        #[case] outcome: Solution<StereoLigandKind, (), ()>,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        policy: StereoFailurePolicy,
    ) {
        let mut entries = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C)
                    .with_implicit_hydrogens(hydrogens)
                    .with_lone_pairs(lone_pairs),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
            ],
            bonds: (1..=3)
                .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
                .collect(),
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Virtual),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        };
        let mut molecule = Molecule::from_entries(entries.clone());
        let original = molecule.clone();
        let resolver = StereoResolver::with_config(
            &stereo_model,
            StereoResolveConfig {
                stereo_atom_failure: policy,
                stereo_bond_failure: policy,
                ring_config: Default::default(),
            },
        );
        let expected = match outcome {
            Solution::Determined(kind) => {
                entries.stereo_atoms[0].1[0].kind = kind;
                Solution::Determined(Molecule::from_entries(entries))
            }
            Solution::Underdetermined(()) => Solution::Underdetermined(()),
            Solution::Contradictory(()) => match policy {
                StereoFailurePolicy::Error => Solution::Contradictory(
                    StereoContradiction::Inconsistency(StereoInconsistency::StereoAtomFailure {
                        stereo_atom: StereoAtomId(0),
                    }),
                ),
                StereoFailurePolicy::Keep => Solution::Determined(original.clone()),
                StereoFailurePolicy::Remove => {
                    entries.stereo_atoms.clear();
                    Solution::Determined(Molecule::from_entries(entries))
                }
            },
        };
        assert_eq!(resolver.resolve(molecule.clone()), expected.clone());
        let in_place = match &expected {
            Solution::Determined(_) => Solution::Determined(()),
            Solution::Underdetermined(()) => Solution::Underdetermined(()),
            Solution::Contradictory(cause) => Solution::Contradictory(cause.clone()),
        };
        assert_eq!(resolver.resolve_into(&mut molecule), in_place);
        if let Solution::Determined(expected) = expected {
            assert_eq!(molecule, expected);
            assert_eq!(resolver.resolve(molecule), Solution::Determined(expected));
        } else {
            assert_eq!(molecule, original);
        }
    }

    #[rstest]
    #[case::atom(
        mol_dsl_concrete!(r#"{:atoms ["C #h3" "S #h0 #T+" "C #h3"]
                             :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        StereoContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) }
        )
    )]
    #[case::bond(
        mol_dsl_concrete!(r#"{:atoms ["C #h3" "C #h2" "C #h1"]
                             :bonds [[0 1 "1"] [1 2 "2#C+"]]}"#),
        StereoContradiction::Inconsistency(
            StereoInconsistency::CisTransStereoFailure { bond: BondId(1) }
        )
    )]
    #[case::atom_open(
        mol_dsl_concrete!(r#"{:atoms ["C#h3" "S#h0#T+" "C#h3"]
                             :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        StereoContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) }
        )
    )]
    #[case::bond_open(
        mol_dsl_concrete!(r#"{:atoms ["C#h3" "C#h2" "C#h"]
                             :bonds [[0 1 "1"] [1 2 "2#C+"]]}"#),
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
            StereoResolver::new(&stereo_model).resolve(molecule.clone()),
            Solution::Contradictory(expected.clone())
        );
        assert_eq!(
            StereoResolver::new(&stereo_model).resolve_into(&mut molecule),
            Solution::Contradictory(expected)
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::actual_ligands(
        mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#T+" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#))]
    #[case::implicit_hydrogen(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs "Th0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T+" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#))]
    #[case::explicit_hydrogen(
        mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [4 1 2 3] :attrs "Th0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#T+" "F" "Cl" "Br" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#))]
    #[case::lone_pair(
        mol_dsl_concrete!(r#"{:atoms ["N#n1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:lp 0] 1 2 3] :attrs "Th0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["N#n1#T+" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#))]
    #[case::tetrahedral_open(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs "Th*"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T+" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#))]
    #[case::existing_assertion(
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T+" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs "Th*"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T+" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#))]
    #[case::bond_side_swap(
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [3 2 4 5] :attrs "Ct0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2#C+"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#))]
    #[case::bond_endpoint_swap(
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[1 0 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [4 5 2 3] :attrs "Ct0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2#C+"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#))]
    #[case::bond_hydrogens(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "C#h1" "F" "Cl"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:h 0] 2 3 [:h 1]] :attrs "Ct0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "C#h1" "F" "Cl"] :bonds [[0 1 "2#C+"] [0 2 "1"] [1 3 "1"]]}"#))]
    #[case::bond_lone_pair(
        mol_dsl_concrete!(r#"{:atoms ["N#n1" "C#h1" "C#h3" "F"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:lp 0] 2 3 [:h 1]] :attrs "Ct0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["N#n1" "C#h1" "C#h3" "F"] :bonds [[0 1 "2#C+"] [0 2 "1"] [1 3 "1"]]}"#))]
    #[case::bond_either(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "C#h1" "F" "Cl"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:h 0] 2 3 [:h 1]] :attrs "Ct*"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "C#h1" "F" "Cl"] :bonds [[0 1 "2#C+"] [0 2 "1"] [1 3 "1"]]}"#))]
    #[case::virtual_ligand(
        mol_dsl_concrete!(r#"{:atoms ["N" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:virt 0]] :attrs "Th?configuration"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["N#T+" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#))]
    #[case::unrealizable_hydrogen(
        mol_dsl_concrete!(r#"{:atoms ["C#h0#n0" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h0#n0#T+" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#))]
    #[case::configuration_set(
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct{0,1}"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2#C+"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#))]
    #[case::invalid_coset(
        mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th2"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#T+" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#))]
    #[case::entity_constraints(
        mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#g/#o(0,1)/#p(0,1)!#f(0,1)"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#T+" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#))]
    fn test_stereo_resolver_project(#[case] mut molecule: Molecule, #[case] expected: Molecule) {
        assert_eq!(
            StereoResolver::project(molecule.clone()),
            Ok(expected.clone())
        );
        assert_eq!(StereoResolver::project_into(&mut molecule), Ok(()));
        assert_eq!(molecule, expected);
        assert_eq!(StereoResolver::project_into(&mut molecule), Ok(()));
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::empty(Molecule::new())]
    #[case::assertions(mol_dsl_concrete!(r#"{:atoms ["C#h1#T+" "C#h1#T!"] :bonds [[0 1 "2#C+"]]}"#))]
    fn test_stereo_resolver_project_identity(#[case] mut molecule: Molecule) {
        let original = molecule.clone();
        assert_eq!(
            StereoResolver::project(molecule.clone()),
            Ok(original.clone())
        );
        assert_eq!(StereoResolver::project_into(&mut molecule), Ok(()));
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
    #[case::atom_assertion(
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T!" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}"#),
        StereoProjectError::AtomAssertion { atom: AtomId(0) })]
    #[case::late_bond_assertion(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "F" "Cl" "Br" "C#h1" "C#h1" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [4 5 "2#C!"] [4 6 "1"] [5 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]
            :stereo-bonds [{:site 3 :ligands [6 [:h 4] 7 [:h 5]] :attrs "Ct1"}]}"#),
        StereoProjectError::BondAssertion { bond: BondId(3) })]
    fn test_stereo_resolver_project_error(
        #[case] mut molecule: Molecule,
        #[case] expected: StereoProjectError,
    ) {
        let original = molecule.clone();
        assert_eq!(
            StereoResolver::project(molecule.clone()),
            Err(expected.clone())
        );
        assert_eq!(StereoResolver::project_into(&mut molecule), Err(expected));
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::tetrahedral("[NH2][C@H](F)Cl",
        mol_dsl_concrete!(r#"{:atoms ["N#h2#n1" "C#h1#T+" "F#n3" "Cl#n3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#))]
    #[case::root_hydrogen("[C@H](F)(Cl)Br",
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T+" "F#n3" "Cl#n3" "Br#n3"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#))]
    #[case::explicit_hydrogen("[H][C@](F)(Cl)Br",
        mol_dsl_concrete!(r#"{:atoms ["H" "C#T+" "F#n3" "Cl#n3" "Br#n3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"] [1 4 "1"]]}"#))]
    #[case::lone_pair("[CH3][S@](=[O])[CH2][CH3]",
        mol_dsl_concrete!(r#"{:atoms ["C#h3" "S#n1#T+" "O#n2" "C#h2" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2"] [1 3 "1"] [3 4 "1"]]}"#))]
    #[case::cis_trans("[CH3]/[CH]=[CH]/[CH3]",
        mol_dsl_concrete!(r#"{:atoms ["C#h3" "C#h1" "C#h1" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C+"] [2 3 "1"]]}"#))]
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
            StereoResolver::project(molecule.clone()),
            Ok(expected.clone())
        );
        assert_eq!(StereoResolver::project_into(&mut molecule), Ok(()));
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::tetrahedral(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs "Th0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h1#T+" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
        Edits::from_iter([
            Edit::ModifyAtomConstraints {
                key: AtomConstraintKey::TetrahedralStereo,
                changes: vec![(
                    AtomHandle::Id(AtomId(0)),
                    None,
                    Some(AtomConstraintForm::TetrahedralStereo(
                        TetrahedralStereoForm::Stereo,
                    )),
                )],
            },
            Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(0)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::ImplicitHydrogen),
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                )],
            },
        ]),
    )]
    #[case::cis_trans(
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "C#h1" "F" "Cl"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [[:h 0] 2 3 [:h 1]] :attrs "Ct0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h1" "C#h1" "F" "Cl"]
            :bonds [[0 1 "2#C+"] [0 2 "1"] [1 3 "1"]]}"#),
        Edits::from_iter([
            Edit::ModifyBondConstraints {
                key: BondConstraintKey::CisTransStereo,
                changes: vec![(
                    BondHandle::Id(BondId(0)),
                    None,
                    Some(BondConstraintForm::CisTransStereo(
                        CisTransStereoForm::Stereo,
                    )),
                )],
            },
            Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(0)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::ImplicitHydrogen),
                        (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::ImplicitHydrogen),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                )],
            },
        ]),
    )]
    fn test_stereo_resolver_plan_project(
        #[case] molecule: Molecule,
        #[case] expected: Molecule,
        #[case] expected_edits: Edits,
    ) {
        let edits = StereoResolver::plan_project(&molecule).unwrap();
        assert_eq!(edits, expected_edits);
        let molecule = molecule.apply(edits).unwrap();
        assert_eq!(molecule, expected);
    }

    #[rstest]
    fn test_stereo_resolver_plan_project_error() {
        let molecule = mol_dsl_concrete!(
            r#"{:atoms ["C#h1#T!" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs "Th0"}]}"#
        );
        assert_eq!(
            StereoResolver::plan_project(&molecule),
            Err(StereoProjectError::AtomAssertion { atom: AtomId(0) }),
        );
    }
}
