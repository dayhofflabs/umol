//! Local stereo integrity checks shared by Molecule and Reaction.

use umol_perm::{Permutation, MAX_DEGREE};

use super::{
    StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset, StereoKind, StereoTerm,
};
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
    StereoCosetOutOfRange {
        entity: Entity,
        kind: StereoKind,
        coset: u32,
        count: usize,
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
    attributes: &StereoAtomForm,
) -> Result<(), StereoIntegrityError> {
    check_stereo_frame(entity, ligand_frame)?;
    if ligand_frame
        .iter()
        .any(|ligand| ligand.kind == StereoLigandKind::Atom && ligand.atom_id == site)
    {
        return Err(StereoIntegrityError::DuplicateAtom { entity, atom: site });
    }
    check_stereo_atom(entity, ligand_frame.len(), attributes)
}

pub(crate) fn check_stereo_bond_entry(
    entity: Entity,
    ligand_frame: &[StereoLigand],
    attributes: &StereoBondForm,
) -> Result<(), StereoIntegrityError> {
    check_stereo_frame(entity, ligand_frame)?;
    check_stereo_bond(entity, ligand_frame.len(), attributes)
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

fn check_stereo_atom(
    entity: Entity,
    ligand_count: usize,
    attributes: &StereoAtomForm,
) -> Result<(), StereoIntegrityError> {
    check_stereo_atom_configuration_on_frame(entity, ligand_count, &attributes.configuration)?;
    for constraint in attributes.constraints.iter() {
        check_stereo_atom_constraint(entity, ligand_count, constraint)?;
    }
    Ok(())
}

fn check_stereo_bond(
    entity: Entity,
    ligand_count: usize,
    attributes: &StereoBondForm,
) -> Result<(), StereoIntegrityError> {
    check_stereo_bond_configuration_on_frame(entity, ligand_count, &attributes.configuration)?;
    for constraint in attributes.constraints.iter() {
        check_stereo_bond_constraint(entity, ligand_count, constraint)?;
    }
    Ok(())
}

pub(crate) fn check_stereo_atom_configuration_on_frame(
    entity: Entity,
    ligand_count: usize,
    configuration: &StereoConfigurationForm,
) -> Result<(), StereoIntegrityError> {
    check_configuration_site_kind(entity, configuration, StereoSite::Atom)?;
    check_configuration(entity, ligand_count, configuration)
}

pub(crate) fn check_stereo_bond_configuration_on_frame(
    entity: Entity,
    ligand_count: usize,
    configuration: &StereoConfigurationForm,
) -> Result<(), StereoIntegrityError> {
    check_configuration_site_kind(entity, configuration, StereoSite::Bond)?;
    check_configuration(entity, ligand_count, configuration)
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

/// A stereo kind describes a coordination geometry, and a geometry belongs to an atom or to a bond.
/// Arity cannot separate them: `Tetrahedral`, `CisTrans`, `Axial`, and `SquarePlanar` all have
/// degree 4. `Axial` is admissible on both, since axial chirality arises at an allene's central
/// atom and about an atropisomeric biaryl bond.
///
/// Matched exhaustively so that a new stereo kind must decide its site here.
fn check_configuration_site_kind(
    entity: Entity,
    configuration: &StereoConfigurationForm,
    site: StereoSite,
) -> Result<(), StereoIntegrityError> {
    let Some(kind) = configuration.kind() else {
        return Ok(());
    };
    check_stereo_site_kind(entity, kind, site)
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

pub(crate) fn check_stereo_frame_arity(
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

pub(crate) fn check_stereo_atom_constraint(
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

pub(crate) fn check_stereo_bond_constraint(
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

fn check_configuration(
    entity: Entity,
    ligand_count: usize,
    configuration: &StereoConfigurationForm,
) -> Result<(), StereoIntegrityError> {
    let StereoConfigurationForm::Kinded(kind, coset) = configuration else {
        return Ok(());
    };
    check_stereo_frame_arity(entity, ligand_count, *kind)?;
    check_coset(entity, *kind, coset)?;
    Ok(())
}

fn check_coset(
    entity: Entity,
    kind: StereoKind,
    coset: &StereoCoset,
) -> Result<(), StereoIntegrityError> {
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
) -> Result<(), StereoIntegrityError> {
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
            check_permutation(entity, kind.degree(), *permutation)?;
            check_term(entity, kind, inner)
        }
    }
}

fn check_coset_index(
    entity: Entity,
    kind: StereoKind,
    coset: u32,
) -> Result<(), StereoIntegrityError> {
    if coset as usize >= kind.count() {
        Err(StereoIntegrityError::StereoCosetOutOfRange {
            entity,
            kind,
            coset,
            count: kind.count(),
        })
    } else {
        Ok(())
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
