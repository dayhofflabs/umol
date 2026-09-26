//! Local stereo integrity checks shared by Molecule and Reaction.

use umol_perm::{Permutation, MAX_DEGREE};

use super::StereoKind;
use crate::ir::constraint::{StereoAtomConstraintForm, StereoBondConstraintForm, StereoLigandPair};
use crate::ir::entity::Entity;
use crate::ir::id::AtomId;
use crate::ir::ligand::{StereoLigand, StereoLigandKind};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum StereoIntegrityError {
    DuplicateAtom {
        entity: Entity,
        atom: AtomId,
    },
    DuplicateStereoLigand {
        entity: Entity,
        ligand: StereoLigand,
    },
    StereoFrameDegreeTooLarge {
        entity: Entity,
        degree: usize,
        maximum: usize,
    },
    StereoKindSiteMismatch {
        entity: Entity,
        kind: StereoKind,
    },
    StereoLigandArity {
        entity: Entity,
        kind: StereoKind,
        expected: usize,
        actual: usize,
    },
    StereoPermutationDegree {
        entity: Entity,
        expected: usize,
        actual: usize,
    },
    StereoLigandPositionOutOfRange {
        entity: Entity,
        position: usize,
        degree: usize,
    },
}

fn check_stereo_frame(
    entity: Entity,
    ligand_frame: &[StereoLigand],
) -> Result<(), StereoIntegrityError> {
    if ligand_frame.len() > MAX_DEGREE {
        return Err(StereoIntegrityError::StereoFrameDegreeTooLarge {
            entity,
            degree: ligand_frame.len(),
            maximum: MAX_DEGREE,
        });
    }

    for (position, &ligand) in ligand_frame.iter().enumerate() {
        if ligand_frame[..position].contains(&ligand) {
            return Err(StereoIntegrityError::DuplicateStereoLigand { entity, ligand });
        }
    }
    Ok(())
}

pub(crate) fn check_stereo_atom_entry(
    entity: Entity,
    site: AtomId,
    ligand_frame: &[StereoLigand],
) -> Result<(), StereoIntegrityError> {
    check_stereo_frame(entity, ligand_frame)?;
    if ligand_frame
        .iter()
        .any(|ligand| ligand.kind == StereoLigandKind::Atom && ligand.atom_id == site)
    {
        return Err(StereoIntegrityError::DuplicateAtom { entity, atom: site });
    }
    Ok(())
}

pub(crate) fn check_stereo_bond_entry(
    entity: Entity,
    ligand_frame: &[StereoLigand],
) -> Result<(), StereoIntegrityError> {
    check_stereo_frame(entity, ligand_frame)
}

pub(crate) fn check_stereo_atom_kind(
    entity: Entity,
    kind: StereoKind,
) -> Result<(), StereoIntegrityError> {
    check_stereo_site_kind(entity, kind, StereoSite::Atom)
}

pub(crate) fn check_stereo_bond_kind(
    entity: Entity,
    kind: StereoKind,
) -> Result<(), StereoIntegrityError> {
    check_stereo_site_kind(entity, kind, StereoSite::Bond)
}

pub(crate) fn check_stereo_atom_constraint_on_frame(
    entity: Entity,
    ligand_count: usize,
    kind: StereoKind,
    constraint: &StereoAtomConstraintForm,
) -> Result<(), StereoIntegrityError> {
    check_stereo_site_kind(entity, kind, StereoSite::Atom)?;
    check_stereo_frame_arity(entity, ligand_count, kind)?;
    check_stereo_atom_constraint(entity, ligand_count, constraint)
}

pub(crate) fn check_stereo_bond_constraint_on_frame(
    entity: Entity,
    ligand_count: usize,
    kind: StereoKind,
    constraint: &StereoBondConstraintForm,
) -> Result<(), StereoIntegrityError> {
    check_stereo_site_kind(entity, kind, StereoSite::Bond)?;
    check_stereo_frame_arity(entity, ligand_count, kind)?;
    check_stereo_bond_constraint(entity, ligand_count, constraint)
}

/// Which site a stereo entry sits on. Local to the admissibility check; the entity id already
/// carries the distinction everywhere else.
#[derive(Clone, Copy)]
enum StereoSite {
    Atom,
    Bond,
}

fn check_stereo_site_kind(
    entity: Entity,
    kind: StereoKind,
    site: StereoSite,
) -> Result<(), StereoIntegrityError> {
    let admissible = match (site, kind) {
        (
            StereoSite::Atom,
            StereoKind::Tetrahedral
            | StereoKind::SquarePlanar
            | StereoKind::TrigonalBipyramidal
            | StereoKind::Octahedral
            | StereoKind::Axial,
        ) => true,
        (StereoSite::Atom, StereoKind::CisTrans) => false,
        (StereoSite::Bond, StereoKind::CisTrans | StereoKind::Axial) => true,
        (
            StereoSite::Bond,
            StereoKind::Tetrahedral
            | StereoKind::SquarePlanar
            | StereoKind::TrigonalBipyramidal
            | StereoKind::Octahedral,
        ) => false,
    };
    if admissible {
        Ok(())
    } else {
        Err(StereoIntegrityError::StereoKindSiteMismatch { entity, kind })
    }
}

fn check_stereo_frame_arity(
    entity: Entity,
    ligand_count: usize,
    kind: StereoKind,
) -> Result<(), StereoIntegrityError> {
    if ligand_count != kind.degree() {
        return Err(StereoIntegrityError::StereoLigandArity {
            entity,
            kind,
            expected: kind.degree(),
            actual: ligand_count,
        });
    }
    Ok(())
}

fn check_stereo_atom_constraint(
    entity: Entity,
    ligand_count: usize,
    constraint: &StereoAtomConstraintForm,
) -> Result<(), StereoIntegrityError> {
    match constraint {
        StereoAtomConstraintForm::LigandSymmetry(value) => {
            check_permutation(entity, ligand_count, value.permutation.permutation.0)
        }
        StereoAtomConstraintForm::Fluxionality(value) => {
            check_permutation(entity, ligand_count, value.permutation.0)
        }
        StereoAtomConstraintForm::Topicity(value) => check_pair(entity, ligand_count, value.pair),
        StereoAtomConstraintForm::Stereogenicity(_) => Ok(()),
    }
}

fn check_stereo_bond_constraint(
    entity: Entity,
    ligand_count: usize,
    constraint: &StereoBondConstraintForm,
) -> Result<(), StereoIntegrityError> {
    match constraint {
        StereoBondConstraintForm::LigandSymmetry(value) => {
            check_permutation(entity, ligand_count, value.permutation.permutation.0)
        }
        StereoBondConstraintForm::Fluxionality(value) => {
            check_permutation(entity, ligand_count, value.permutation.0)
        }
        StereoBondConstraintForm::Topicity(value) => check_pair(entity, ligand_count, value.pair),
        StereoBondConstraintForm::Stereogenicity(_) => Ok(()),
    }
}

fn check_permutation(
    entity: Entity,
    expected: usize,
    permutation: Permutation,
) -> Result<(), StereoIntegrityError> {
    if permutation.degree() != expected {
        Err(StereoIntegrityError::StereoPermutationDegree {
            entity,
            expected,
            actual: permutation.degree(),
        })
    } else {
        Ok(())
    }
}

fn check_pair(
    entity: Entity,
    degree: usize,
    pair: StereoLigandPair,
) -> Result<(), StereoIntegrityError> {
    for position in [pair.first(), pair.second()] {
        if position.index() >= degree {
            return Err(StereoIntegrityError::StereoLigandPositionOutOfRange {
                entity,
                position: position.index(),
                degree,
            });
        }
    }
    Ok(())
}
