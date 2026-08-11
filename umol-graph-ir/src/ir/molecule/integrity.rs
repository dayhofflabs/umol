//! Representation-integrity checks for [`Molecule`].

use thiserror::Error;
use umol_perm::Permutation;

use super::super::constraint::{
    StereoAtomConstraintForm, StereoBondConstraintForm, StereoLigandPair,
};
use super::super::electrons::ElectronCountsForm;
use super::super::entity::Entity;
use super::super::stereo::{StereoConfigurationForm, StereoCoset, StereoKind, StereoTerm};
use super::Molecule;

/// Failure of the representation contract required to interpret a [`Molecule`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MoleculeIntegrityError {
    #[error("molecule references unavailable {entity}")]
    InvalidReference { entity: Entity },
    #[error("molecule stores {actual} {kind:?} values for {expected} graph entities")]
    EntityCountMismatch {
        kind: super::super::entity::EntityKind,
        expected: usize,
        actual: usize,
    },
    #[error(
        "{entity}: electron-count vector has length {electron_counts}, expected {participants}"
    )]
    ElectronCountLengthMismatch {
        entity: Entity,
        participants: usize,
        electron_counts: usize,
    },
    #[error("{entity}: stereo kind {kind:?} is not permitted for this entity family")]
    StereoKindNotPermitted { entity: Entity, kind: StereoKind },
    #[error("{entity}: stereo kind is required to interpret positional constraints")]
    StereoKindRequired { entity: Entity },
    #[error("{entity}: stereo frame has {actual} ligands, expected {expected} for {kind:?}")]
    StereoLigandArity {
        entity: Entity,
        kind: StereoKind,
        expected: usize,
        actual: usize,
    },
    #[error("{entity}: coset {coset} is outside 0..{count} for {kind:?}")]
    StereoCosetOutOfRange {
        entity: Entity,
        kind: StereoKind,
        coset: u32,
        count: usize,
    },
    #[error(
        "{entity}: permutation has degree {actual}, expected {expected} for stereo kind {kind:?}"
    )]
    StereoPermutationDegree {
        entity: Entity,
        kind: StereoKind,
        expected: usize,
        actual: usize,
    },
    #[error("{entity}: ligand position {position} is outside 0..{degree}")]
    StereoLigandPositionOutOfRange {
        entity: Entity,
        position: usize,
        degree: usize,
    },
}

impl Molecule {
    /// Check the representation invariants required to interpret this molecule.
    ///
    /// This checks stored references, parallel collection shapes, and kind-dependent stereo
    /// domains. It does not check graph simplicity, entity uniqueness, chemistry, or constraint
    /// satisfaction.
    pub fn check_integrity(&self) -> Result<(), MoleculeIntegrityError> {
        use super::super::entity::EntityKind;

        if self.atoms.len() != self.graph.node_count() {
            return Err(MoleculeIntegrityError::EntityCountMismatch {
                kind: EntityKind::Atom,
                expected: self.graph.node_count(),
                actual: self.atoms.len(),
            });
        }
        if self.bonds.len() != self.graph.edge_count() {
            return Err(MoleculeIntegrityError::EntityCountMismatch {
                kind: EntityKind::Bond,
                expected: self.graph.edge_count(),
                actual: self.bonds.len(),
            });
        }

        let contains = |entity| match entity {
            Entity::Atom(id) => id.index() < self.atoms.len(),
            Entity::Bond(id) => id.index() < self.bonds.len(),
            Entity::DativeBond(id) => id.index() < self.dative_bonds.count(),
            Entity::AromaticSystem(id) => id.index() < self.aromatic_systems.count(),
            Entity::MulticenterBond(id) => id.index() < self.multicenter_bonds.count(),
            Entity::NoncovalentBond(id) => id.index() < self.noncovalent_bonds.count(),
            Entity::StereoAtom(id) => id.index() < self.stereo_atoms.count(),
            Entity::StereoBond(id) => id.index() < self.stereo_bonds.count(),
        };

        for view in self.bonds().iter() {
            for atom in view.atom_ids() {
                require_reference(&contains, Entity::Atom(atom))?;
            }
        }
        for view in self.dative_bonds().iter() {
            require_reference(&contains, Entity::Atom(view.acceptor_id()))?;
            require_references(&contains, view.donor_ids().map(Entity::Atom))?;
        }
        for view in self.aromatic_systems().iter() {
            require_references(&contains, view.atom_ids().map(Entity::Atom))?;
            check_electron_count_length(
                Entity::AromaticSystem(view.id),
                view.atom_ids().count(),
                &view.attributes.electrons,
            )?;
        }
        for view in self.multicenter_bonds().iter() {
            require_references(&contains, view.atom_ids().map(Entity::Atom))?;
            check_electron_count_length(
                Entity::MulticenterBond(view.id),
                view.atom_ids().count(),
                &view.attributes.electrons,
            )?;
        }
        for view in self.noncovalent_bonds().iter() {
            require_references(&contains, view.atom_ids().into_iter().map(Entity::Atom))?;
        }
        for view in self.stereo_atoms().iter() {
            let entity = Entity::StereoAtom(view.id);
            require_reference(&contains, Entity::Atom(view.site_id()))?;
            require_references(
                &contains,
                view.ligand_frame()
                    .into_iter()
                    .map(|ligand| Entity::Atom(ligand.atom_id)),
            )?;
            check_stereo_atom(entity, view.ligands().count(), view.attributes)?;
        }
        for view in self.stereo_bonds().iter() {
            let entity = Entity::StereoBond(view.id);
            require_reference(&contains, Entity::Bond(view.site_id()))?;
            require_references(
                &contains,
                view.ligand_frame()
                    .into_iter()
                    .map(|ligand| Entity::Atom(ligand.atom_id)),
            )?;
            check_stereo_bond(entity, view.ligands().count(), view.attributes)?;
        }
        for constraint in self.constraints.iter() {
            super::validate_constraint_references(constraint, &contains)
                .map_err(|entity| MoleculeIntegrityError::InvalidReference { entity })?;
        }
        Ok(())
    }
}

pub(super) fn require_reference(
    contains: &dyn Fn(Entity) -> bool,
    entity: Entity,
) -> Result<(), MoleculeIntegrityError> {
    if contains(entity) {
        Ok(())
    } else {
        Err(MoleculeIntegrityError::InvalidReference { entity })
    }
}

pub(super) fn require_references(
    contains: &dyn Fn(Entity) -> bool,
    entities: impl IntoIterator<Item = Entity>,
) -> Result<(), MoleculeIntegrityError> {
    for entity in entities {
        require_reference(contains, entity)?;
    }
    Ok(())
}

fn check_electron_count_length(
    entity: Entity,
    participants: usize,
    electrons: &ElectronCountsForm,
) -> Result<(), MoleculeIntegrityError> {
    if let ElectronCountsForm::Lit(counts) = electrons {
        if counts.len() != participants {
            return Err(MoleculeIntegrityError::ElectronCountLengthMismatch {
                entity,
                participants,
                electron_counts: counts.len(),
            });
        }
    }
    Ok(())
}

fn check_stereo_atom(
    entity: Entity,
    ligand_count: usize,
    attributes: &super::super::stereo::StereoAtomForm,
) -> Result<(), MoleculeIntegrityError> {
    let kind = check_configuration(entity, ligand_count, &attributes.configuration, |kind| {
        matches!(
            kind,
            StereoKind::Tetrahedral
                | StereoKind::SquarePlanar
                | StereoKind::TrigonalBipyramidal
                | StereoKind::Octahedral
        )
    })?;
    for constraint in attributes.constraints.iter() {
        match constraint {
            StereoAtomConstraintForm::LigandSymmetry(value) => check_permutation(
                entity,
                kind_required(entity, kind)?,
                value.permutation.permutation.0,
            )?,
            StereoAtomConstraintForm::Fluxionality(value) => {
                check_permutation(entity, kind_required(entity, kind)?, value.permutation.0)?
            }
            StereoAtomConstraintForm::Topicity(value) => {
                check_pair(entity, kind_required(entity, kind)?, value.pair)?
            }
            StereoAtomConstraintForm::Stereogenicity(_) => {}
        }
    }
    Ok(())
}

fn check_stereo_bond(
    entity: Entity,
    ligand_count: usize,
    attributes: &super::super::stereo::StereoBondForm,
) -> Result<(), MoleculeIntegrityError> {
    let kind = check_configuration(entity, ligand_count, &attributes.configuration, |kind| {
        kind == StereoKind::CisTrans
    })?;
    for constraint in attributes.constraints.iter() {
        match constraint {
            StereoBondConstraintForm::LigandSymmetry(value) => check_permutation(
                entity,
                kind_required(entity, kind)?,
                value.permutation.permutation.0,
            )?,
            StereoBondConstraintForm::Fluxionality(value) => {
                check_permutation(entity, kind_required(entity, kind)?, value.permutation.0)?
            }
            StereoBondConstraintForm::Topicity(value) => {
                check_pair(entity, kind_required(entity, kind)?, value.pair)?
            }
            StereoBondConstraintForm::Stereogenicity(_) => {}
        }
    }
    Ok(())
}

fn check_configuration(
    entity: Entity,
    ligand_count: usize,
    configuration: &StereoConfigurationForm,
    permitted: impl FnOnce(StereoKind) -> bool,
) -> Result<Option<StereoKind>, MoleculeIntegrityError> {
    let StereoConfigurationForm::Kinded(kind, coset) = configuration else {
        return Ok(None);
    };
    if !permitted(*kind) {
        return Err(MoleculeIntegrityError::StereoKindNotPermitted {
            entity,
            kind: *kind,
        });
    }
    if ligand_count != kind.degree() {
        return Err(MoleculeIntegrityError::StereoLigandArity {
            entity,
            kind: *kind,
            expected: kind.degree(),
            actual: ligand_count,
        });
    }
    check_coset(entity, *kind, coset)?;
    Ok(Some(*kind))
}

fn check_coset(
    entity: Entity,
    kind: StereoKind,
    coset: &StereoCoset,
) -> Result<(), MoleculeIntegrityError> {
    match coset {
        StereoCoset::Undetermined => Ok(()),
        StereoCoset::Lit(value) => check_coset_index(entity, kind, *value),
        StereoCoset::LitSet(values) => {
            for &value in values {
                check_coset_index(entity, kind, value)?;
            }
            Ok(())
        }
        StereoCoset::Term(term) => check_term(entity, kind, term),
    }
}

fn check_term(
    entity: Entity,
    kind: StereoKind,
    term: &StereoTerm,
) -> Result<(), MoleculeIntegrityError> {
    match term {
        StereoTerm::Var(value) => {
            if let Some(domain) = &value.1 {
                for &coset in domain {
                    check_coset_index(entity, kind, coset)?;
                }
            }
            Ok(())
        }
        StereoTerm::Lit(value) => check_coset_index(entity, kind, *value),
        StereoTerm::LitSet(values) => {
            for &value in values {
                check_coset_index(entity, kind, value)?;
            }
            Ok(())
        }
        StereoTerm::Swap(inner) | StereoTerm::Mirror(inner) => check_term(entity, kind, inner),
        StereoTerm::Apply(inner, permutation) => {
            check_permutation(entity, kind, *permutation)?;
            check_term(entity, kind, inner)
        }
    }
}

fn check_coset_index(
    entity: Entity,
    kind: StereoKind,
    coset: u32,
) -> Result<(), MoleculeIntegrityError> {
    if coset as usize >= kind.count() {
        Err(MoleculeIntegrityError::StereoCosetOutOfRange {
            entity,
            kind,
            coset,
            count: kind.count(),
        })
    } else {
        Ok(())
    }
}

fn kind_required(
    entity: Entity,
    kind: Option<StereoKind>,
) -> Result<StereoKind, MoleculeIntegrityError> {
    kind.ok_or(MoleculeIntegrityError::StereoKindRequired { entity })
}

fn check_permutation(
    entity: Entity,
    kind: StereoKind,
    permutation: Permutation,
) -> Result<(), MoleculeIntegrityError> {
    if permutation.degree() != kind.degree() {
        Err(MoleculeIntegrityError::StereoPermutationDegree {
            entity,
            kind,
            expected: kind.degree(),
            actual: permutation.degree(),
        })
    } else {
        Ok(())
    }
}

fn check_pair(
    entity: Entity,
    kind: StereoKind,
    pair: StereoLigandPair,
) -> Result<(), MoleculeIntegrityError> {
    for position in [pair.first(), pair.second()] {
        if position.index() >= kind.degree() {
            return Err(MoleculeIntegrityError::StereoLigandPositionOutOfRange {
                entity,
                position: position.index(),
                degree: kind.degree(),
            });
        }
    }
    Ok(())
}
