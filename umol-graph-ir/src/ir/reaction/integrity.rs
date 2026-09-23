//! Reaction representation-integrity checks.

use std::collections::HashMap;
use std::iter;

use thiserror::Error;

use super::super::constraint::{Constraint, MoleculeConstraint, RelationalConstraint};
use super::super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, ConstraintDelta, DativeBondDelta, Delta, Deltas,
    MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta,
};
use super::super::edit::{
    AromaticSystemFieldChange, MulticenterBondFieldChange, StereoAtomFieldChange,
    StereoBondFieldChange,
};
use super::super::electrons::ElectronCountsForm;
use super::super::entity::Entity;
use super::super::id::AtomId;
use super::super::ligand::StereoLigand;
use super::super::molecule::Molecule;
use super::super::stereo::integrity::{
    check_stereo_atom_entry, check_stereo_atom_kind, check_stereo_bond_entry,
    check_stereo_bond_kind, StereoIntegrityError,
};
use super::super::stereo::StereoKind;
use super::Reaction;

/// Internal implementation of reaction integrity checking.
#[derive(Clone, Copy, Debug, Default)]
struct ReactionIntegrityCheck;

/// Failure of the representation contract required to interpret a [`Reaction`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionIntegrityError {
    /// A delta refers to an entity unavailable from either the lhs or the reaction's additions.
    #[error("reaction references unavailable entity {entity:?}")]
    InvalidReference { entity: Entity },
    /// An addition reuses an entity ID from the lhs or an earlier addition.
    #[error("reaction adds duplicate entity reference {entity:?}")]
    DuplicateReference { entity: Entity },
    /// A literal electron-count vector has a different length from its participant frame.
    #[error(
        "{entity}: electron-count vector has length {electron_counts}, expected {participants}"
    )]
    ElectronCountLengthMismatch {
        entity: Entity,
        participants: usize,
        electron_counts: usize,
    },
    /// An atom occurs twice in one stereo entity's atom references.
    #[error("{entity}: participant atom {atom:?} is duplicated")]
    DuplicateAtom { entity: Entity, atom: AtomId },
    /// A stereo frame repeats a complete ligand value.
    #[error("{entity}: stereo ligand {ligand:?} is duplicated in the frame")]
    DuplicateStereoLigand {
        entity: Entity,
        ligand: StereoLigand,
    },
    /// A stereo frame exceeds the supported permutation degree.
    #[error(
        "{entity}: stereo frame has degree {degree}, exceeding the supported maximum {maximum}"
    )]
    StereoFrameDegreeTooLarge {
        entity: Entity,
        degree: usize,
        maximum: usize,
    },
    /// A stereo kind cannot be borne by its site type.
    #[error("{entity}: stereo kind {kind:?} is not admissible for this site type")]
    StereoKindSiteMismatch { entity: Entity, kind: StereoKind },
    /// A stereo frame's length differs from its declared kind's degree.
    #[error("{entity}: stereo frame has {actual} ligands, expected {expected} for {kind:?}")]
    StereoLigandArity {
        entity: Entity,
        kind: StereoKind,
        expected: usize,
        actual: usize,
    },
    /// A coset index is outside the declared kind's range.
    #[error("{entity}: coset {coset} is outside 0..{count} for {kind:?}")]
    StereoCosetOutOfRange {
        entity: Entity,
        kind: StereoKind,
        coset: u32,
        count: usize,
    },
    /// A permutation has a different degree from the frame or kind it acts on.
    #[error(
        "{entity}: permutation has degree {actual}, expected {expected} for the stored ligand frame"
    )]
    StereoPermutationDegree {
        entity: Entity,
        expected: usize,
        actual: usize,
    },
    /// A topicity pair names a position outside its stereo frame.
    #[error("{entity}: ligand position {position} is outside 0..{degree}")]
    StereoLigandPositionOutOfRange {
        entity: Entity,
        position: usize,
        degree: usize,
    },
    /// A removal records incidence incompatible with its source entity's participant structure.
    #[error("reaction incidence does not match source entity {entity:?}")]
    IncidenceMismatch { entity: Entity },
    /// A configuration change replaces one stereo kind with another within a single entity.
    #[error("{entity:?}: configuration change replaces stereo kind {old:?} with {new:?}")]
    StereoKindModified {
        entity: Entity,
        old: StereoKind,
        new: StereoKind,
    },
}

impl From<StereoIntegrityError> for ReactionIntegrityError {
    fn from(error: StereoIntegrityError) -> Self {
        match error {
            StereoIntegrityError::DuplicateAtom { entity, atom } => {
                Self::DuplicateAtom { entity, atom }
            }
            StereoIntegrityError::DuplicateStereoLigand { entity, ligand } => {
                Self::DuplicateStereoLigand { entity, ligand }
            }
            StereoIntegrityError::StereoFrameDegreeTooLarge {
                entity,
                degree,
                maximum,
            } => Self::StereoFrameDegreeTooLarge {
                entity,
                degree,
                maximum,
            },
            StereoIntegrityError::StereoKindSiteMismatch { entity, kind } => {
                Self::StereoKindSiteMismatch { entity, kind }
            }
            StereoIntegrityError::StereoLigandArity {
                entity,
                kind,
                expected,
                actual,
            } => Self::StereoLigandArity {
                entity,
                kind,
                expected,
                actual,
            },
            StereoIntegrityError::StereoCosetOutOfRange {
                entity,
                kind,
                coset,
                count,
            } => Self::StereoCosetOutOfRange {
                entity,
                kind,
                coset,
                count,
            },
            StereoIntegrityError::StereoPermutationDegree {
                entity,
                expected,
                actual,
            } => Self::StereoPermutationDegree {
                entity,
                expected,
                actual,
            },
            StereoIntegrityError::StereoLigandPositionOutOfRange {
                entity,
                position,
                degree,
            } => Self::StereoLigandPositionOutOfRange {
                entity,
                position,
                degree,
            },
        }
    }
}

/// A stereo entity keeps its kind across a configuration change: the kind names the coordination
/// geometry, so replacing it replaces the stereogenic unit rather than its configuration. That is
/// expressed as removal plus addition, where the two entities carry different ids.
///
/// An undetermined side asserts no geometry and so restricts nothing.
fn check_delta_stereo_kind(
    entity: Entity,
    old: Option<StereoKind>,
    new: Option<StereoKind>,
) -> Result<(), ReactionIntegrityError> {
    match (old, new) {
        (Some(old), Some(new)) if old != new => {
            Err(ReactionIntegrityError::StereoKindModified { entity, old, new })
        }
        _ => Ok(()),
    }
}

impl ReactionIntegrityCheck {
    fn check(&self, lhs: &Molecule, deltas: &Deltas) -> Result<(), ReactionIntegrityError> {
        let mut added = HashMap::new();
        for delta in deltas.iter() {
            if let Some(entity) = added_entity(delta) {
                if contains_entity(lhs, entity) || added.insert(entity, delta).is_some() {
                    return Err(ReactionIntegrityError::DuplicateReference { entity });
                }
            }
        }

        for delta in deltas.iter() {
            self.validate_references(lhs, &added, delta)?;
        }
        for delta in deltas.iter() {
            self.validate_electron_count_delta(lhs, &added, delta)?;
        }
        for delta in deltas.iter() {
            self.validate_stereo_delta(delta)?;
        }
        for delta in deltas.iter() {
            self.validate_removal_incidence(lhs, &added, delta)?;
        }
        Ok(())
    }

    fn validate_references(
        &self,
        lhs: &Molecule,
        added: &HashMap<Entity, &Delta>,
        delta: &Delta,
    ) -> Result<(), ReactionIntegrityError> {
        match delta {
            Delta::Atom(AtomDelta::Add { .. }) => Ok(()),
            Delta::Atom(
                AtomDelta::Remove { id, .. }
                | AtomDelta::ModifyField { id, .. }
                | AtomDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, added, Entity::Atom(*id)),
            Delta::Bond(BondDelta::Add { atoms, .. }) => self.require_atoms(lhs, added, *atoms),
            Delta::Bond(BondDelta::Remove { id, atoms, .. }) => {
                self.require_available(lhs, added, Entity::Bond(*id))?;
                self.require_atoms(lhs, added, *atoms)
            }
            Delta::Bond(
                BondDelta::ModifyField { id, .. } | BondDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, added, Entity::Bond(*id)),
            Delta::DativeBond(DativeBondDelta::Add {
                donors, acceptor, ..
            }) => self.require_atoms(lhs, added, donors.iter().copied().chain([*acceptor])),
            Delta::DativeBond(DativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                ..
            }) => {
                self.require_available(lhs, added, Entity::DativeBond(*id))?;
                self.require_atoms(lhs, added, donors.iter().copied().chain([*acceptor]))
            }
            Delta::DativeBond(
                DativeBondDelta::ModifyField { id, .. }
                | DativeBondDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, added, Entity::DativeBond(*id)),
            Delta::AromaticSystem(AromaticSystemDelta::Add { atoms, .. }) => {
                self.require_atoms(lhs, added, atoms.iter().copied())
            }
            Delta::AromaticSystem(AromaticSystemDelta::Remove { id, atoms, .. }) => {
                self.require_available(lhs, added, Entity::AromaticSystem(*id))?;
                self.require_atoms(lhs, added, atoms.iter().copied())
            }
            Delta::AromaticSystem(
                AromaticSystemDelta::ModifyField { id, .. }
                | AromaticSystemDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, added, Entity::AromaticSystem(*id)),
            Delta::MulticenterBond(MulticenterBondDelta::Add { atoms, .. }) => {
                self.require_atoms(lhs, added, atoms.iter().copied())
            }
            Delta::MulticenterBond(MulticenterBondDelta::Remove { id, atoms, .. }) => {
                self.require_available(lhs, added, Entity::MulticenterBond(*id))?;
                self.require_atoms(lhs, added, atoms.iter().copied())
            }
            Delta::MulticenterBond(
                MulticenterBondDelta::ModifyField { id, .. }
                | MulticenterBondDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, added, Entity::MulticenterBond(*id)),
            Delta::NoncovalentBond(NoncovalentBondDelta::Add { atoms, .. }) => {
                self.require_atoms(lhs, added, *atoms)
            }
            Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, atoms, .. }) => {
                self.require_available(lhs, added, Entity::NoncovalentBond(*id))?;
                self.require_atoms(lhs, added, *atoms)
            }
            Delta::NoncovalentBond(
                NoncovalentBondDelta::ModifyField { id, .. }
                | NoncovalentBondDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, added, Entity::NoncovalentBond(*id)),
            Delta::StereoAtom(StereoAtomDelta::Add { site, ligands, .. }) => self.require_atoms(
                lhs,
                added,
                iter::once(*site).chain(ligands.iter().map(|ligand| ligand.atom_id)),
            ),
            Delta::StereoAtom(StereoAtomDelta::Remove {
                id, site, ligands, ..
            }) => {
                self.require_available(lhs, added, Entity::StereoAtom(*id))?;
                self.require_atoms(
                    lhs,
                    added,
                    iter::once(*site).chain(ligands.iter().map(|ligand| ligand.atom_id)),
                )
            }
            Delta::StereoAtom(
                StereoAtomDelta::ModifyField { id, .. }
                | StereoAtomDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, added, Entity::StereoAtom(*id)),
            Delta::StereoBond(StereoBondDelta::Add { site, ligands, .. }) => {
                self.require_available(lhs, added, Entity::Bond(*site))?;
                self.require_atoms(lhs, added, ligands.iter().map(|ligand| ligand.atom_id))
            }
            Delta::StereoBond(StereoBondDelta::Remove {
                id, site, ligands, ..
            }) => {
                self.require_available(lhs, added, Entity::StereoBond(*id))?;
                self.require_available(lhs, added, Entity::Bond(*site))?;
                self.require_atoms(lhs, added, ligands.iter().map(|ligand| ligand.atom_id))
            }
            Delta::StereoBond(
                StereoBondDelta::ModifyField { id, .. }
                | StereoBondDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, added, Entity::StereoBond(*id)),
            Delta::Constraint(ConstraintDelta::Add(constraint))
            | Delta::Constraint(ConstraintDelta::Remove(constraint)) => {
                self.validate_constraint(lhs, added, constraint)
            }
        }
    }

    fn require_available(
        &self,
        lhs: &Molecule,
        added: &HashMap<Entity, &Delta>,
        entity: Entity,
    ) -> Result<(), ReactionIntegrityError> {
        if contains_entity(lhs, entity) || added.contains_key(&entity) {
            Ok(())
        } else {
            Err(ReactionIntegrityError::InvalidReference { entity })
        }
    }

    fn require_atoms(
        &self,
        lhs: &Molecule,
        added: &HashMap<Entity, &Delta>,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Result<(), ReactionIntegrityError> {
        for atom in atoms {
            self.require_available(lhs, added, Entity::Atom(atom))?;
        }
        Ok(())
    }

    fn validate_electron_count_delta(
        &self,
        lhs: &Molecule,
        added: &HashMap<Entity, &Delta>,
        delta: &Delta,
    ) -> Result<(), ReactionIntegrityError> {
        match delta {
            Delta::AromaticSystem(
                AromaticSystemDelta::Add {
                    id,
                    atoms,
                    attributes,
                }
                | AromaticSystemDelta::Remove {
                    id,
                    atoms,
                    attributes,
                },
            ) => check_electron_count_length(
                Entity::AromaticSystem(*id),
                atoms.len(),
                &attributes.electrons,
            ),
            Delta::MulticenterBond(
                MulticenterBondDelta::Add {
                    id,
                    atoms,
                    attributes,
                }
                | MulticenterBondDelta::Remove {
                    id,
                    atoms,
                    attributes,
                },
            ) => check_electron_count_length(
                Entity::MulticenterBond(*id),
                atoms.len(),
                &attributes.electrons,
            ),
            Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
                id,
                change: AromaticSystemFieldChange::Electrons { old, new },
            }) => {
                let entity = Entity::AromaticSystem(*id);
                let participants = if let Some(view) = lhs.aromatic_systems().get(*id) {
                    view.atom_ids().len()
                } else {
                    let Delta::AromaticSystem(AromaticSystemDelta::Add { atoms, .. }) =
                        added[&entity]
                    else {
                        unreachable!("reference check found the added aromatic system")
                    };
                    atoms.len()
                };
                check_electron_count_length(entity, participants, old)?;
                check_electron_count_length(entity, participants, new)
            }
            Delta::MulticenterBond(MulticenterBondDelta::ModifyField {
                id,
                change: MulticenterBondFieldChange::Electrons { old, new },
            }) => {
                let entity = Entity::MulticenterBond(*id);
                let participants = if let Some(view) = lhs.multicenter_bonds().get(*id) {
                    view.atom_ids().len()
                } else {
                    let Delta::MulticenterBond(MulticenterBondDelta::Add { atoms, .. }) =
                        added[&entity]
                    else {
                        unreachable!("reference check found the added multicenter bond")
                    };
                    atoms.len()
                };
                check_electron_count_length(entity, participants, old)?;
                check_electron_count_length(entity, participants, new)
            }
            _ => Ok(()),
        }
    }

    fn validate_stereo_delta(&self, delta: &Delta) -> Result<(), ReactionIntegrityError> {
        let result = match delta {
            Delta::StereoAtom(StereoAtomDelta::Add {
                id,
                site,
                ligands,
                attributes,
            }) => check_stereo_atom_entry(Entity::StereoAtom(*id), *site, ligands, attributes),
            Delta::StereoAtom(StereoAtomDelta::ModifyConstraint {
                id,
                kind: Some(kind),
                ..
            }) => check_stereo_atom_kind(Entity::StereoAtom(*id), *kind),
            Delta::StereoBond(StereoBondDelta::Add {
                id,
                ligands,
                attributes,
                ..
            }) => check_stereo_bond_entry(Entity::StereoBond(*id), ligands, attributes),
            Delta::StereoBond(StereoBondDelta::ModifyConstraint {
                id,
                kind: Some(kind),
                ..
            }) => check_stereo_bond_kind(Entity::StereoBond(*id), *kind),
            _ => Ok(()),
        };
        result.map_err(ReactionIntegrityError::from)?;

        match delta {
            Delta::StereoAtom(StereoAtomDelta::ModifyField { id, change }) => {
                let StereoAtomFieldChange::Configuration { old, new } = change;
                check_delta_stereo_kind(Entity::StereoAtom(*id), old.kind(), new.kind())
            }
            Delta::StereoBond(StereoBondDelta::ModifyField { id, change }) => {
                let StereoBondFieldChange::Configuration { old, new } = change;
                check_delta_stereo_kind(Entity::StereoBond(*id), old.kind(), new.kind())
            }
            _ => Ok(()),
        }
    }

    fn validate_removal_incidence(
        &self,
        lhs: &Molecule,
        added: &HashMap<Entity, &Delta>,
        delta: &Delta,
    ) -> Result<(), ReactionIntegrityError> {
        let (entity, matches) = match delta {
            Delta::Bond(BondDelta::Remove { id, atoms, .. }) => {
                let entity = Entity::Bond(*id);
                let source = if let Some(view) = lhs.bonds().get(*id) {
                    view.atom_ids()
                } else {
                    let Delta::Bond(BondDelta::Add { atoms, .. }) = added[&entity] else {
                        unreachable!("reference check found the added bond")
                    };
                    *atoms
                };
                (entity, unordered_pair(source) == unordered_pair(*atoms))
            }
            Delta::DativeBond(DativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                ..
            }) => {
                let entity = Entity::DativeBond(*id);
                let (source_acceptor, source_donors) =
                    if let Some(view) = lhs.dative_bonds().get(*id) {
                        (view.acceptor_id(), unordered_ids(view.donor_ids()))
                    } else {
                        let Delta::DativeBond(DativeBondDelta::Add {
                            donors, acceptor, ..
                        }) = added[&entity]
                        else {
                            unreachable!("reference check found the added dative bond")
                        };
                        (*acceptor, unordered_ids(donors.iter().copied()))
                    };
                (
                    entity,
                    source_acceptor == *acceptor
                        && source_donors == unordered_ids(donors.iter().copied()),
                )
            }
            Delta::AromaticSystem(AromaticSystemDelta::Remove { id, atoms, .. }) => {
                let entity = Entity::AromaticSystem(*id);
                let source = if let Some(view) = lhs.aromatic_systems().get(*id) {
                    unordered_ids(view.atom_ids())
                } else {
                    let Delta::AromaticSystem(AromaticSystemDelta::Add { atoms, .. }) =
                        added[&entity]
                    else {
                        unreachable!("reference check found the added aromatic system")
                    };
                    unordered_ids(atoms.iter().copied())
                };
                (entity, source == unordered_ids(atoms.iter().copied()))
            }
            Delta::MulticenterBond(MulticenterBondDelta::Remove { id, atoms, .. }) => {
                let entity = Entity::MulticenterBond(*id);
                let source = if let Some(view) = lhs.multicenter_bonds().get(*id) {
                    unordered_ids(view.atom_ids())
                } else {
                    let Delta::MulticenterBond(MulticenterBondDelta::Add { atoms, .. }) =
                        added[&entity]
                    else {
                        unreachable!("reference check found the added multicenter bond")
                    };
                    unordered_ids(atoms.iter().copied())
                };
                (entity, source == unordered_ids(atoms.iter().copied()))
            }
            Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, atoms, .. }) => {
                let entity = Entity::NoncovalentBond(*id);
                let source = if let Some(view) = lhs.noncovalent_bonds().get(*id) {
                    view.atom_ids()
                } else {
                    let Delta::NoncovalentBond(NoncovalentBondDelta::Add { atoms, .. }) =
                        added[&entity]
                    else {
                        unreachable!("reference check found the added noncovalent bond")
                    };
                    *atoms
                };
                (entity, unordered_pair(source) == unordered_pair(*atoms))
            }
            Delta::StereoAtom(StereoAtomDelta::Remove {
                id, site, ligands, ..
            }) => {
                let entity = Entity::StereoAtom(*id);
                let matches = if let Some(view) = lhs.stereo_atoms().get(*id) {
                    view.site_id() == *site
                        && unordered_ligands(view.ligand_frame())
                            == unordered_ligands(ligands.iter().copied())
                } else {
                    let Delta::StereoAtom(StereoAtomDelta::Add {
                        site: added_site,
                        ligands: added_ligands,
                        ..
                    }) = added[&entity]
                    else {
                        unreachable!("reference check found the added stereo atom")
                    };
                    *added_site == *site
                        && unordered_ligands(added_ligands.iter().copied())
                            == unordered_ligands(ligands.iter().copied())
                };
                (entity, matches)
            }
            Delta::StereoBond(StereoBondDelta::Remove {
                id, site, ligands, ..
            }) => {
                let entity = Entity::StereoBond(*id);
                let matches = if let Some(view) = lhs.stereo_bonds().get(*id) {
                    view.site_id() == *site
                        && stereo_bond_frames_match(&view.ligand_frame(), ligands)
                } else {
                    let Delta::StereoBond(StereoBondDelta::Add {
                        site: added_site,
                        ligands: added_ligands,
                        ..
                    }) = added[&entity]
                    else {
                        unreachable!("reference check found the added stereo bond")
                    };
                    *added_site == *site && stereo_bond_frames_match(added_ligands, ligands)
                };
                (entity, matches)
            }
            _ => return Ok(()),
        };
        if matches {
            Ok(())
        } else {
            Err(ReactionIntegrityError::IncidenceMismatch { entity })
        }
    }

    fn validate_constraint(
        &self,
        lhs: &Molecule,
        added: &HashMap<Entity, &Delta>,
        constraint: &Constraint,
    ) -> Result<(), ReactionIntegrityError> {
        match constraint {
            Constraint::Atom(id, _) => self.require_available(lhs, added, Entity::Atom(*id)),
            Constraint::Bond(id, _) => self.require_available(lhs, added, Entity::Bond(*id)),
            Constraint::DativeBond(id, _) => {
                self.require_available(lhs, added, Entity::DativeBond(*id))
            }
            Constraint::AromaticSystem(id, _) => {
                self.require_available(lhs, added, Entity::AromaticSystem(*id))
            }
            Constraint::MulticenterBond(id, _) => {
                self.require_available(lhs, added, Entity::MulticenterBond(*id))
            }
            Constraint::NoncovalentBond(id, _) => {
                self.require_available(lhs, added, Entity::NoncovalentBond(*id))
            }
            Constraint::StereoAtom(id, kind, _) => {
                self.require_available(lhs, added, Entity::StereoAtom(*id))?;
                check_stereo_atom_kind(Entity::StereoAtom(*id), *kind)
                    .map_err(ReactionIntegrityError::from)
            }
            Constraint::StereoBond(id, kind, _) => {
                self.require_available(lhs, added, Entity::StereoBond(*id))?;
                check_stereo_bond_kind(Entity::StereoBond(*id), *kind)
                    .map_err(ReactionIntegrityError::from)
            }
            Constraint::Relational(constraint) => {
                self.validate_relational_constraint(lhs, added, constraint)
            }
            Constraint::Molecule(constraint) => {
                self.validate_molecule_constraint(lhs, added, constraint)
            }
            Constraint::And(constraints) | Constraint::Or(constraints) => {
                for constraint in constraints {
                    self.validate_constraint(lhs, added, constraint)?;
                }
                Ok(())
            }
            Constraint::Not(constraint) => self.validate_constraint(lhs, added, constraint),
        }
    }

    fn validate_relational_constraint(
        &self,
        lhs: &Molecule,
        added: &HashMap<Entity, &Delta>,
        constraint: &RelationalConstraint,
    ) -> Result<(), ReactionIntegrityError> {
        match constraint {
            RelationalConstraint::DativeBondDonors { bond, atoms }
            | RelationalConstraint::DativeBondContainsAllDonors { bond, atoms } => {
                self.require_available(lhs, added, Entity::DativeBond(*bond))?;
                self.require_atoms(lhs, added, atoms.iter().copied())
            }
            RelationalConstraint::DativeBondDonor { bond, atom }
            | RelationalConstraint::DativeBondAcceptor { bond, atom } => {
                self.require_available(lhs, added, Entity::DativeBond(*bond))?;
                self.require_available(lhs, added, Entity::Atom(*atom))
            }
            RelationalConstraint::DativeBondAllDonors { bond, .. }
            | RelationalConstraint::DativeBondAnyDonor { bond, .. }
            | RelationalConstraint::DativeBondAcceptorSatisfies { bond, .. } => {
                self.require_available(lhs, added, Entity::DativeBond(*bond))
            }
            RelationalConstraint::DativeBondParallels { dative, parallel } => {
                self.require_available(lhs, added, Entity::DativeBond(*dative))?;
                self.require_available(lhs, added, Entity::Bond(*parallel))
            }
            RelationalConstraint::AromaticSystemAtoms { system, atoms }
            | RelationalConstraint::AromaticSystemContainsAll { system, atoms } => {
                self.require_available(lhs, added, Entity::AromaticSystem(*system))?;
                self.require_atoms(lhs, added, atoms.iter().copied())
            }
            RelationalConstraint::AromaticSystemContains { system, atom } => {
                self.require_available(lhs, added, Entity::AromaticSystem(*system))?;
                self.require_available(lhs, added, Entity::Atom(*atom))
            }
            RelationalConstraint::AromaticSystemAllAtoms { system, .. }
            | RelationalConstraint::AromaticSystemAnyAtom { system, .. } => {
                self.require_available(lhs, added, Entity::AromaticSystem(*system))
            }
            RelationalConstraint::MulticenterBondAtoms { bond, atoms }
            | RelationalConstraint::MulticenterBondContainsAll { bond, atoms } => {
                self.require_available(lhs, added, Entity::MulticenterBond(*bond))?;
                self.require_atoms(lhs, added, atoms.iter().copied())
            }
            RelationalConstraint::MulticenterBondContains { bond, atom } => {
                self.require_available(lhs, added, Entity::MulticenterBond(*bond))?;
                self.require_available(lhs, added, Entity::Atom(*atom))
            }
            RelationalConstraint::MulticenterBondAllAtoms { bond, .. }
            | RelationalConstraint::MulticenterBondAnyAtom { bond, .. } => {
                self.require_available(lhs, added, Entity::MulticenterBond(*bond))
            }
            RelationalConstraint::NoncovalentBondEnds { bond, atoms } => {
                self.require_available(lhs, added, Entity::NoncovalentBond(*bond))?;
                self.require_atoms(lhs, added, *atoms)
            }
            RelationalConstraint::NoncovalentBondContains { bond, atom } => {
                self.require_available(lhs, added, Entity::NoncovalentBond(*bond))?;
                self.require_available(lhs, added, Entity::Atom(*atom))
            }
            RelationalConstraint::NoncovalentBondEndsSatisfy { bond, .. } => {
                self.require_available(lhs, added, Entity::NoncovalentBond(*bond))
            }
            RelationalConstraint::StereoAtomSite { stereo_atom, atom }
            | RelationalConstraint::StereoAtomContains { stereo_atom, atom } => {
                self.require_available(lhs, added, Entity::StereoAtom(*stereo_atom))?;
                self.require_available(lhs, added, Entity::Atom(*atom))
            }
            RelationalConstraint::StereoAtomLigands { stereo_atom, atoms } => {
                self.require_available(lhs, added, Entity::StereoAtom(*stereo_atom))?;
                self.require_atoms(lhs, added, atoms.iter().copied())
            }
            RelationalConstraint::StereoAtomAllLigands { stereo_atom, .. }
            | RelationalConstraint::StereoAtomAnyLigand { stereo_atom, .. } => {
                self.require_available(lhs, added, Entity::StereoAtom(*stereo_atom))
            }
            RelationalConstraint::StereoBondSite { stereo_bond, bond } => {
                self.require_available(lhs, added, Entity::StereoBond(*stereo_bond))?;
                self.require_available(lhs, added, Entity::Bond(*bond))
            }
            RelationalConstraint::StereoBondContains { stereo_bond, atom } => {
                self.require_available(lhs, added, Entity::StereoBond(*stereo_bond))?;
                self.require_available(lhs, added, Entity::Atom(*atom))
            }
            RelationalConstraint::StereoBondLigands { stereo_bond, atoms } => {
                self.require_available(lhs, added, Entity::StereoBond(*stereo_bond))?;
                self.require_atoms(lhs, added, atoms.iter().copied())
            }
            RelationalConstraint::StereoBondAllLigands { stereo_bond, .. }
            | RelationalConstraint::StereoBondAnyLigand { stereo_bond, .. } => {
                self.require_available(lhs, added, Entity::StereoBond(*stereo_bond))
            }
        }
    }

    fn validate_molecule_constraint(
        &self,
        lhs: &Molecule,
        added: &HashMap<Entity, &Delta>,
        constraint: &MoleculeConstraint,
    ) -> Result<(), ReactionIntegrityError> {
        match constraint {
            MoleculeConstraint::ChargeSum { atoms, .. }
            | MoleculeConstraint::UnpairedElectronCoupling { atoms, .. }
            | MoleculeConstraint::Connected { atoms } => {
                self.require_atoms(lhs, added, atoms.iter().flatten().copied())
            }
            MoleculeConstraint::BondOrderSum { bonds, .. } => {
                for &bond in bonds.iter().flatten() {
                    self.require_available(lhs, added, Entity::Bond(bond))?;
                }
                Ok(())
            }
        }
    }
}

impl Reaction {
    /// Check the representation invariants required to interpret this reaction.
    ///
    /// The check covers delta references, added-id uniqueness, positional electron-count lengths,
    /// local stereo data carried by additions and constraint wrappers, and the source incidence and
    /// participant structure recorded by removals. The closed lhs already satisfies molecule
    /// integrity. This check does not impose DPO or chemistry semantics.
    pub(crate) fn check_integrity(&self) -> Result<(), ReactionIntegrityError> {
        ReactionIntegrityCheck.check(&self.lhs, &self.deltas)
    }
}

fn added_entity(delta: &Delta) -> Option<Entity> {
    match delta {
        Delta::Atom(AtomDelta::Add { id, .. }) => Some(Entity::Atom(*id)),
        Delta::Bond(BondDelta::Add { id, .. }) => Some(Entity::Bond(*id)),
        Delta::DativeBond(DativeBondDelta::Add { id, .. }) => Some(Entity::DativeBond(*id)),
        Delta::AromaticSystem(AromaticSystemDelta::Add { id, .. }) => {
            Some(Entity::AromaticSystem(*id))
        }
        Delta::MulticenterBond(MulticenterBondDelta::Add { id, .. }) => {
            Some(Entity::MulticenterBond(*id))
        }
        Delta::NoncovalentBond(NoncovalentBondDelta::Add { id, .. }) => {
            Some(Entity::NoncovalentBond(*id))
        }
        Delta::StereoAtom(StereoAtomDelta::Add { id, .. }) => Some(Entity::StereoAtom(*id)),
        Delta::StereoBond(StereoBondDelta::Add { id, .. }) => Some(Entity::StereoBond(*id)),
        _ => None,
    }
}

fn check_electron_count_length(
    entity: Entity,
    participants: usize,
    electrons: &ElectronCountsForm,
) -> Result<(), ReactionIntegrityError> {
    if let ElectronCountsForm::Lit(counts) = electrons {
        if counts.len() != participants {
            return Err(ReactionIntegrityError::ElectronCountLengthMismatch {
                entity,
                participants,
                electron_counts: counts.len(),
            });
        }
    }
    Ok(())
}

fn contains_entity(molecule: &Molecule, entity: Entity) -> bool {
    match entity {
        Entity::Atom(id) => molecule.atoms().contains(id),
        Entity::Bond(id) => molecule.bonds().contains(id),
        Entity::DativeBond(id) => molecule.dative_bonds().contains(id),
        Entity::AromaticSystem(id) => molecule.aromatic_systems().contains(id),
        Entity::MulticenterBond(id) => molecule.multicenter_bonds().contains(id),
        Entity::NoncovalentBond(id) => molecule.noncovalent_bonds().contains(id),
        Entity::StereoAtom(id) => molecule.stereo_atoms().contains(id),
        Entity::StereoBond(id) => molecule.stereo_bonds().contains(id),
    }
}

fn unordered_ids(ids: impl IntoIterator<Item = AtomId>) -> Vec<AtomId> {
    let mut ids: Vec<AtomId> = ids.into_iter().collect();
    ids.sort_unstable();
    ids
}

fn unordered_ligands(ligands: impl IntoIterator<Item = StereoLigand>) -> Vec<StereoLigand> {
    let mut ligands: Vec<StereoLigand> = ligands.into_iter().collect();
    ligands.sort_unstable();
    ligands
}

/// Stereo-bond ligand frames have two endpoint blocks. A compatible local frame may reorder each
/// block and may swap the two complete blocks, but may not move one ligand across the partition.
fn stereo_bond_frames_match(source: &[StereoLigand], local: &[StereoLigand]) -> bool {
    if source.len() != 4 || local.len() != 4 {
        return false;
    }
    let source_first = unordered_ligands(source[..2].iter().copied());
    let source_second = unordered_ligands(source[2..].iter().copied());
    let local_first = unordered_ligands(local[..2].iter().copied());
    let local_second = unordered_ligands(local[2..].iter().copied());
    (source_first == local_first && source_second == local_second)
        || (source_first == local_second && source_second == local_first)
}

fn unordered_pair(mut atoms: [AtomId; 2]) -> [AtomId; 2] {
    atoms.sort_unstable();
    atoms
}
