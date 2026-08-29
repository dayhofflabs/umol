//! Reaction representation-integrity checks.

use std::collections::{HashMap, HashSet};
use std::iter;

use thiserror::Error;

use super::super::constraint::{Constraint, MoleculeConstraint, RelationalConstraint};
use super::super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, ConstraintDelta, DativeBondDelta, Delta, Deltas,
    MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta,
};
use super::super::edit::{StereoAtomFieldChange, StereoBondFieldChange};
use super::super::entity::Entity;
use super::super::id::{AtomId, BondId};
use super::super::ligand::StereoLigand;
use super::super::molecule::integrity::{
    check_stereo_atom_entry, check_stereo_atom_kind, check_stereo_bond_entry,
    check_stereo_bond_kind,
};
use super::super::molecule::{Molecule, MoleculeIntegrityError};
use super::super::stereo::StereoKind;
use super::Reaction;

/// Internal implementation of reaction integrity checking.
#[derive(Clone, Copy, Debug, Default)]
struct ReactionIntegrityCheck;

#[derive(Clone, Debug)]
enum OverlayFrame {
    Dative {
        donors: Vec<AtomId>,
        acceptor: AtomId,
    },
    Aromatic(Vec<AtomId>),
    Multicenter(Vec<AtomId>),
    Noncovalent([AtomId; 2]),
    StereoAtom {
        site: AtomId,
        ligands: Vec<StereoLigand>,
    },
    StereoBond {
        site: BondId,
        ligands: Vec<StereoLigand>,
    },
}

/// Failure of the representation contract required to interpret a [`Reaction`].
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionIntegrityError {
    /// A delta refers to an entity unavailable from either the lhs or the reaction's additions.
    #[error("reaction references unavailable entity {entity:?}")]
    InvalidReference { entity: Entity },
    /// Stereo data carried by the reaction violates a local representation invariant.
    #[error("reaction stereo representation is invalid: {0}")]
    StereoIntegrityError(MoleculeIntegrityError),
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
        let mut source_frames = source_frames(lhs);
        let mut created = HashSet::new();
        for delta in deltas.iter() {
            if let Some((entity, frame)) = added_entity_and_frame(delta) {
                if contains_entity(lhs, entity) || !created.insert(entity) {
                    return Err(ReactionIntegrityError::InvalidReference { entity });
                }
                if let Some(frame) = frame {
                    source_frames.insert(entity, frame);
                }
            }
        }

        for delta in deltas.iter() {
            self.validate_references(lhs, &created, delta)?;
        }
        for delta in deltas.iter() {
            self.validate_stereo_delta(delta)?;
        }
        for delta in deltas.iter() {
            self.validate_removal_incidence(lhs, &source_frames, delta)?;
        }
        Ok(())
    }

    fn validate_references(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &Delta,
    ) -> Result<(), ReactionIntegrityError> {
        match delta {
            Delta::Atom(AtomDelta::Add { .. }) => Ok(()),
            Delta::Atom(
                AtomDelta::Remove { id, .. }
                | AtomDelta::ModifyField { id, .. }
                | AtomDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, created, Entity::Atom(*id)),
            Delta::Bond(BondDelta::Add { atoms, .. }) => self.require_atoms(lhs, created, *atoms),
            Delta::Bond(BondDelta::Remove { id, atoms, .. }) => {
                self.require_available(lhs, created, Entity::Bond(*id))?;
                self.require_atoms(lhs, created, *atoms)
            }
            Delta::Bond(
                BondDelta::ModifyField { id, .. } | BondDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, created, Entity::Bond(*id)),
            Delta::DativeBond(DativeBondDelta::Add {
                donors, acceptor, ..
            }) => self.require_atoms(lhs, created, donors.iter().copied().chain([*acceptor])),
            Delta::DativeBond(DativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                ..
            }) => {
                self.require_available(lhs, created, Entity::DativeBond(*id))?;
                self.require_atoms(lhs, created, donors.iter().copied().chain([*acceptor]))
            }
            Delta::DativeBond(
                DativeBondDelta::ModifyField { id, .. }
                | DativeBondDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, created, Entity::DativeBond(*id)),
            Delta::AromaticSystem(AromaticSystemDelta::Add { atoms, .. }) => {
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            Delta::AromaticSystem(AromaticSystemDelta::Remove { id, atoms, .. }) => {
                self.require_available(lhs, created, Entity::AromaticSystem(*id))?;
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            Delta::AromaticSystem(
                AromaticSystemDelta::ModifyField { id, .. }
                | AromaticSystemDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, created, Entity::AromaticSystem(*id)),
            Delta::MulticenterBond(MulticenterBondDelta::Add { atoms, .. }) => {
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            Delta::MulticenterBond(MulticenterBondDelta::Remove { id, atoms, .. }) => {
                self.require_available(lhs, created, Entity::MulticenterBond(*id))?;
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            Delta::MulticenterBond(
                MulticenterBondDelta::ModifyField { id, .. }
                | MulticenterBondDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, created, Entity::MulticenterBond(*id)),
            Delta::NoncovalentBond(NoncovalentBondDelta::Add { atoms, .. }) => {
                self.require_atoms(lhs, created, *atoms)
            }
            Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, atoms, .. }) => {
                self.require_available(lhs, created, Entity::NoncovalentBond(*id))?;
                self.require_atoms(lhs, created, *atoms)
            }
            Delta::NoncovalentBond(
                NoncovalentBondDelta::ModifyField { id, .. }
                | NoncovalentBondDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, created, Entity::NoncovalentBond(*id)),
            Delta::StereoAtom(StereoAtomDelta::Add { site, ligands, .. }) => self.require_atoms(
                lhs,
                created,
                iter::once(*site).chain(ligands.iter().map(|ligand| ligand.atom_id)),
            ),
            Delta::StereoAtom(StereoAtomDelta::Remove {
                id, site, ligands, ..
            }) => {
                self.require_available(lhs, created, Entity::StereoAtom(*id))?;
                self.require_atoms(
                    lhs,
                    created,
                    iter::once(*site).chain(ligands.iter().map(|ligand| ligand.atom_id)),
                )
            }
            Delta::StereoAtom(
                StereoAtomDelta::ModifyField { id, .. }
                | StereoAtomDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, created, Entity::StereoAtom(*id)),
            Delta::StereoBond(StereoBondDelta::Add { site, ligands, .. }) => {
                self.require_available(lhs, created, Entity::Bond(*site))?;
                self.require_atoms(lhs, created, ligands.iter().map(|ligand| ligand.atom_id))
            }
            Delta::StereoBond(StereoBondDelta::Remove {
                id, site, ligands, ..
            }) => {
                self.require_available(lhs, created, Entity::StereoBond(*id))?;
                self.require_available(lhs, created, Entity::Bond(*site))?;
                self.require_atoms(lhs, created, ligands.iter().map(|ligand| ligand.atom_id))
            }
            Delta::StereoBond(
                StereoBondDelta::ModifyField { id, .. }
                | StereoBondDelta::ModifyConstraint { id, .. },
            ) => self.require_available(lhs, created, Entity::StereoBond(*id)),
            Delta::Constraint(ConstraintDelta::Add(constraint))
            | Delta::Constraint(ConstraintDelta::Remove(constraint)) => {
                self.validate_constraint(lhs, created, constraint)
            }
        }
    }

    fn require_available(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        entity: Entity,
    ) -> Result<(), ReactionIntegrityError> {
        if contains_entity(lhs, entity) || created.contains(&entity) {
            Ok(())
        } else {
            Err(ReactionIntegrityError::InvalidReference { entity })
        }
    }

    fn require_atoms(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Result<(), ReactionIntegrityError> {
        for atom in atoms {
            self.require_available(lhs, created, Entity::Atom(atom))?;
        }
        Ok(())
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
        result.map_err(ReactionIntegrityError::StereoIntegrityError)?;

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
        source_frames: &HashMap<Entity, OverlayFrame>,
        delta: &Delta,
    ) -> Result<(), ReactionIntegrityError> {
        let (entity, matches) = match delta {
            Delta::Bond(BondDelta::Remove { id, atoms, .. }) => {
                let entity = Entity::Bond(*id);
                return if !lhs.bonds().contains(*id)
                    || unordered_pair(lhs.bonds().get(*id).expect("checked lhs bond").atom_ids())
                        == unordered_pair(*atoms)
                {
                    Ok(())
                } else {
                    Err(ReactionIntegrityError::IncidenceMismatch { entity })
                };
            }
            Delta::DativeBond(DativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                ..
            }) => {
                let entity = Entity::DativeBond(*id);
                let OverlayFrame::Dative {
                    donors: source,
                    acceptor: source_acceptor,
                } = &source_frames[&entity]
                else {
                    unreachable!("entity kind fixes its source-frame variant")
                };
                (
                    entity,
                    *source_acceptor == *acceptor
                        && unordered_ids(source.iter().copied())
                            == unordered_ids(donors.iter().copied()),
                )
            }
            Delta::AromaticSystem(AromaticSystemDelta::Remove { id, atoms, .. }) => {
                let entity = Entity::AromaticSystem(*id);
                let OverlayFrame::Aromatic(source) = &source_frames[&entity] else {
                    unreachable!("entity kind fixes its source-frame variant")
                };
                (
                    entity,
                    unordered_ids(source.iter().copied()) == unordered_ids(atoms.iter().copied()),
                )
            }
            Delta::MulticenterBond(MulticenterBondDelta::Remove { id, atoms, .. }) => {
                let entity = Entity::MulticenterBond(*id);
                let OverlayFrame::Multicenter(source) = &source_frames[&entity] else {
                    unreachable!("entity kind fixes its source-frame variant")
                };
                (
                    entity,
                    unordered_ids(source.iter().copied()) == unordered_ids(atoms.iter().copied()),
                )
            }
            Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, atoms, .. }) => {
                let entity = Entity::NoncovalentBond(*id);
                let OverlayFrame::Noncovalent(source) = source_frames[&entity] else {
                    unreachable!("entity kind fixes its source-frame variant")
                };
                (entity, unordered_pair(source) == unordered_pair(*atoms))
            }
            Delta::StereoAtom(StereoAtomDelta::Remove {
                id, site, ligands, ..
            }) => {
                let entity = Entity::StereoAtom(*id);
                let OverlayFrame::StereoAtom {
                    site: source_site,
                    ligands: source,
                } = &source_frames[&entity]
                else {
                    unreachable!("entity kind fixes its source-frame variant")
                };
                (
                    entity,
                    *source_site == *site
                        && unordered_ligands(source.iter().copied())
                            == unordered_ligands(ligands.iter().copied()),
                )
            }
            Delta::StereoBond(StereoBondDelta::Remove {
                id, site, ligands, ..
            }) => {
                let entity = Entity::StereoBond(*id);
                let OverlayFrame::StereoBond {
                    site: source_site,
                    ligands: source,
                } = &source_frames[&entity]
                else {
                    unreachable!("entity kind fixes its source-frame variant")
                };
                (
                    entity,
                    *source_site == *site && stereo_bond_frames_match(source, ligands),
                )
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
        created: &HashSet<Entity>,
        constraint: &Constraint,
    ) -> Result<(), ReactionIntegrityError> {
        match constraint {
            Constraint::Atom(id, _) => self.require_available(lhs, created, Entity::Atom(*id)),
            Constraint::Bond(id, _) => self.require_available(lhs, created, Entity::Bond(*id)),
            Constraint::DativeBond(id, _) => {
                self.require_available(lhs, created, Entity::DativeBond(*id))
            }
            Constraint::AromaticSystem(id, _) => {
                self.require_available(lhs, created, Entity::AromaticSystem(*id))
            }
            Constraint::MulticenterBond(id, _) => {
                self.require_available(lhs, created, Entity::MulticenterBond(*id))
            }
            Constraint::NoncovalentBond(id, _) => {
                self.require_available(lhs, created, Entity::NoncovalentBond(*id))
            }
            Constraint::StereoAtom(id, kind, _) => {
                self.require_available(lhs, created, Entity::StereoAtom(*id))?;
                check_stereo_atom_kind(Entity::StereoAtom(*id), *kind)
                    .map_err(ReactionIntegrityError::StereoIntegrityError)
            }
            Constraint::StereoBond(id, kind, _) => {
                self.require_available(lhs, created, Entity::StereoBond(*id))?;
                check_stereo_bond_kind(Entity::StereoBond(*id), *kind)
                    .map_err(ReactionIntegrityError::StereoIntegrityError)
            }
            Constraint::Relational(constraint) => {
                self.validate_relational_constraint(lhs, created, constraint)
            }
            Constraint::Molecule(constraint) => {
                self.validate_molecule_constraint(lhs, created, constraint)
            }
            Constraint::And(constraints) | Constraint::Or(constraints) => {
                for constraint in constraints {
                    self.validate_constraint(lhs, created, constraint)?;
                }
                Ok(())
            }
            Constraint::Not(constraint) => self.validate_constraint(lhs, created, constraint),
        }
    }

    fn validate_relational_constraint(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        constraint: &RelationalConstraint,
    ) -> Result<(), ReactionIntegrityError> {
        match constraint {
            RelationalConstraint::DativeBondDonors { bond, atoms }
            | RelationalConstraint::DativeBondContainsAllDonors { bond, atoms } => {
                self.require_available(lhs, created, Entity::DativeBond(*bond))?;
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            RelationalConstraint::DativeBondDonor { bond, atom }
            | RelationalConstraint::DativeBondAcceptor { bond, atom } => {
                self.require_available(lhs, created, Entity::DativeBond(*bond))?;
                self.require_available(lhs, created, Entity::Atom(*atom))
            }
            RelationalConstraint::DativeBondAllDonors { bond, .. }
            | RelationalConstraint::DativeBondAnyDonor { bond, .. }
            | RelationalConstraint::DativeBondAcceptorSatisfies { bond, .. } => {
                self.require_available(lhs, created, Entity::DativeBond(*bond))
            }
            RelationalConstraint::DativeBondParallels { dative, parallel } => {
                self.require_available(lhs, created, Entity::DativeBond(*dative))?;
                self.require_available(lhs, created, Entity::Bond(*parallel))
            }
            RelationalConstraint::AromaticSystemAtoms { system, atoms }
            | RelationalConstraint::AromaticSystemContainsAll { system, atoms } => {
                self.require_available(lhs, created, Entity::AromaticSystem(*system))?;
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            RelationalConstraint::AromaticSystemContains { system, atom } => {
                self.require_available(lhs, created, Entity::AromaticSystem(*system))?;
                self.require_available(lhs, created, Entity::Atom(*atom))
            }
            RelationalConstraint::AromaticSystemAllAtoms { system, .. }
            | RelationalConstraint::AromaticSystemAnyAtom { system, .. } => {
                self.require_available(lhs, created, Entity::AromaticSystem(*system))
            }
            RelationalConstraint::MulticenterBondAtoms { bond, atoms }
            | RelationalConstraint::MulticenterBondContainsAll { bond, atoms } => {
                self.require_available(lhs, created, Entity::MulticenterBond(*bond))?;
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            RelationalConstraint::MulticenterBondContains { bond, atom } => {
                self.require_available(lhs, created, Entity::MulticenterBond(*bond))?;
                self.require_available(lhs, created, Entity::Atom(*atom))
            }
            RelationalConstraint::MulticenterBondAllAtoms { bond, .. }
            | RelationalConstraint::MulticenterBondAnyAtom { bond, .. } => {
                self.require_available(lhs, created, Entity::MulticenterBond(*bond))
            }
            RelationalConstraint::NoncovalentBondEnds { bond, atoms } => {
                self.require_available(lhs, created, Entity::NoncovalentBond(*bond))?;
                self.require_atoms(lhs, created, *atoms)
            }
            RelationalConstraint::NoncovalentBondContains { bond, atom } => {
                self.require_available(lhs, created, Entity::NoncovalentBond(*bond))?;
                self.require_available(lhs, created, Entity::Atom(*atom))
            }
            RelationalConstraint::NoncovalentBondEndsSatisfy { bond, .. } => {
                self.require_available(lhs, created, Entity::NoncovalentBond(*bond))
            }
            RelationalConstraint::StereoAtomSite { stereo_atom, atom }
            | RelationalConstraint::StereoAtomContains { stereo_atom, atom } => {
                self.require_available(lhs, created, Entity::StereoAtom(*stereo_atom))?;
                self.require_available(lhs, created, Entity::Atom(*atom))
            }
            RelationalConstraint::StereoAtomLigands { stereo_atom, atoms } => {
                self.require_available(lhs, created, Entity::StereoAtom(*stereo_atom))?;
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            RelationalConstraint::StereoAtomAllLigands { stereo_atom, .. }
            | RelationalConstraint::StereoAtomAnyLigand { stereo_atom, .. } => {
                self.require_available(lhs, created, Entity::StereoAtom(*stereo_atom))
            }
            RelationalConstraint::StereoBondSite { stereo_bond, bond } => {
                self.require_available(lhs, created, Entity::StereoBond(*stereo_bond))?;
                self.require_available(lhs, created, Entity::Bond(*bond))
            }
            RelationalConstraint::StereoBondContains { stereo_bond, atom } => {
                self.require_available(lhs, created, Entity::StereoBond(*stereo_bond))?;
                self.require_available(lhs, created, Entity::Atom(*atom))
            }
            RelationalConstraint::StereoBondLigands { stereo_bond, atoms } => {
                self.require_available(lhs, created, Entity::StereoBond(*stereo_bond))?;
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            RelationalConstraint::StereoBondAllLigands { stereo_bond, .. }
            | RelationalConstraint::StereoBondAnyLigand { stereo_bond, .. } => {
                self.require_available(lhs, created, Entity::StereoBond(*stereo_bond))
            }
        }
    }

    fn validate_molecule_constraint(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        constraint: &MoleculeConstraint,
    ) -> Result<(), ReactionIntegrityError> {
        match constraint {
            MoleculeConstraint::ChargeSum { atoms, .. }
            | MoleculeConstraint::UnpairedElectronCoupling { atoms, .. }
            | MoleculeConstraint::Connected { atoms } => {
                self.require_atoms(lhs, created, atoms.iter().flatten().copied())
            }
            MoleculeConstraint::BondOrderSum { bonds, .. } => {
                for &bond in bonds.iter().flatten() {
                    self.require_available(lhs, created, Entity::Bond(bond))?;
                }
                Ok(())
            }
        }
    }
}

impl Reaction {
    /// Check the representation invariants required to interpret this reaction.
    ///
    /// The check covers delta references, created-id uniqueness, local stereo data carried by
    /// additions and constraint wrappers, and the source incidence and participant structure
    /// recorded by removals. The closed lhs already satisfies molecule integrity. This check does
    /// not impose DPO or chemistry semantics.
    pub(crate) fn check_integrity(&self) -> Result<(), ReactionIntegrityError> {
        ReactionIntegrityCheck.check(&self.lhs, &self.deltas)
    }
}

fn source_frames(lhs: &Molecule) -> HashMap<Entity, OverlayFrame> {
    let mut frames = HashMap::new();
    for view in lhs.dative_bonds().iter() {
        frames.insert(
            Entity::DativeBond(view.id),
            OverlayFrame::Dative {
                donors: view.donor_ids().collect(),
                acceptor: view.acceptor_id(),
            },
        );
    }
    for view in lhs.aromatic_systems().iter() {
        frames.insert(
            Entity::AromaticSystem(view.id),
            OverlayFrame::Aromatic(view.atom_ids().collect()),
        );
    }
    for view in lhs.multicenter_bonds().iter() {
        frames.insert(
            Entity::MulticenterBond(view.id),
            OverlayFrame::Multicenter(view.atom_ids().collect()),
        );
    }
    for view in lhs.noncovalent_bonds().iter() {
        frames.insert(
            Entity::NoncovalentBond(view.id),
            OverlayFrame::Noncovalent(view.atom_ids()),
        );
    }
    for view in lhs.stereo_atoms().iter() {
        frames.insert(
            Entity::StereoAtom(view.id),
            OverlayFrame::StereoAtom {
                site: view.site_id(),
                ligands: view.ligand_frame(),
            },
        );
    }
    for view in lhs.stereo_bonds().iter() {
        frames.insert(
            Entity::StereoBond(view.id),
            OverlayFrame::StereoBond {
                site: view.site_id(),
                ligands: view.ligand_frame(),
            },
        );
    }
    frames
}

fn added_entity_and_frame(delta: &Delta) -> Option<(Entity, Option<OverlayFrame>)> {
    match delta {
        Delta::Atom(AtomDelta::Add { id, .. }) => Some((Entity::Atom(*id), None)),
        Delta::Bond(BondDelta::Add { id, .. }) => Some((Entity::Bond(*id), None)),
        Delta::DativeBond(DativeBondDelta::Add {
            id,
            donors,
            acceptor,
            ..
        }) => Some((
            Entity::DativeBond(*id),
            Some(OverlayFrame::Dative {
                donors: donors.clone(),
                acceptor: *acceptor,
            }),
        )),
        Delta::AromaticSystem(AromaticSystemDelta::Add { id, atoms, .. }) => Some((
            Entity::AromaticSystem(*id),
            Some(OverlayFrame::Aromatic(atoms.clone())),
        )),
        Delta::MulticenterBond(MulticenterBondDelta::Add { id, atoms, .. }) => Some((
            Entity::MulticenterBond(*id),
            Some(OverlayFrame::Multicenter(atoms.clone())),
        )),
        Delta::NoncovalentBond(NoncovalentBondDelta::Add { id, atoms, .. }) => Some((
            Entity::NoncovalentBond(*id),
            Some(OverlayFrame::Noncovalent(*atoms)),
        )),
        Delta::StereoAtom(StereoAtomDelta::Add {
            id, site, ligands, ..
        }) => Some((
            Entity::StereoAtom(*id),
            Some(OverlayFrame::StereoAtom {
                site: *site,
                ligands: ligands.clone(),
            }),
        )),
        Delta::StereoBond(StereoBondDelta::Add {
            id, site, ligands, ..
        }) => Some((
            Entity::StereoBond(*id),
            Some(OverlayFrame::StereoBond {
                site: *site,
                ligands: ligands.clone(),
            }),
        )),
        _ => None,
    }
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
