//! Reaction integrity validation for entity references and structural incidence.

use std::collections::HashSet;
use std::iter;

use thiserror::Error;

use super::super::constraint::{Constraint, MoleculeConstraint, RelationalConstraint};
use super::super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, ConstraintDelta, DativeBondDelta, Delta, Deltas,
    MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta,
};
use super::super::entity::Entity;
use super::super::id::AtomId;
use super::super::ligand::StereoLigand;
use super::super::molecule::{Molecule, MoleculeIntegrityError};
use super::super::reaction::Reaction;

/// Internal implementation of reaction integrity checking.
#[derive(Clone, Copy, Debug, Default)]
struct ReactionIntegrityCheck;

#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionIntegrityError {
    #[error("reaction lhs is not a valid molecule representation: {0}")]
    Lhs(#[from] MoleculeIntegrityError),
    #[error("reaction references unavailable entity {entity:?}")]
    InvalidReference { entity: Entity },
    #[error("reaction incidence does not match lhs entity {entity:?}")]
    IncidenceMismatch { entity: Entity },
}

impl ReactionIntegrityCheck {
    fn check(&self, lhs: &Molecule, deltas: &Deltas) -> Result<(), ReactionIntegrityError> {
        let mut created = HashSet::new();
        for delta in deltas.iter() {
            let entity = match delta {
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
            };
            if let Some(entity) = entity {
                if contains_entity(lhs, entity) || !created.insert(entity) {
                    return Err(ReactionIntegrityError::InvalidReference { entity });
                }
            }
        }
        for delta in deltas.iter() {
            match delta {
                Delta::Atom(delta) => self.validate_atom(lhs, &created, delta)?,
                Delta::Bond(delta) => self.validate_bond(lhs, &created, delta)?,
                Delta::DativeBond(delta) => self.validate_dative(lhs, &created, delta)?,
                Delta::AromaticSystem(delta) => self.validate_aromatic(lhs, &created, delta)?,
                Delta::MulticenterBond(delta) => self.validate_multicenter(lhs, &created, delta)?,
                Delta::NoncovalentBond(delta) => self.validate_noncovalent(lhs, &created, delta)?,
                Delta::StereoAtom(delta) => self.validate_stereo_atom(lhs, &created, delta)?,
                Delta::StereoBond(delta) => self.validate_stereo_bond(lhs, &created, delta)?,
                Delta::Constraint(ConstraintDelta::Add(constraint))
                | Delta::Constraint(ConstraintDelta::Remove(constraint)) => {
                    self.validate_constraint(lhs, &created, constraint)?;
                }
            }
        }
        Ok(())
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

    fn incidence_mismatch(entity: Entity) -> ReactionIntegrityError {
        ReactionIntegrityError::IncidenceMismatch { entity }
    }

    fn validate_atom(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &AtomDelta,
    ) -> Result<(), ReactionIntegrityError> {
        match delta {
            AtomDelta::Add { .. } => Ok(()),
            AtomDelta::Remove { id, .. }
            | AtomDelta::ModifyField { id, .. }
            | AtomDelta::ModifyConstraint { id, .. } => {
                self.require_available(lhs, created, Entity::Atom(*id))
            }
        }
    }

    fn validate_bond(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &BondDelta,
    ) -> Result<(), ReactionIntegrityError> {
        match delta {
            BondDelta::Add { atoms, .. } => self.require_atoms(lhs, created, *atoms),
            BondDelta::Remove { id, atoms, .. } => {
                let entity = Entity::Bond(*id);
                self.require_available(lhs, created, entity)?;
                self.require_atoms(lhs, created, *atoms)?;
                if !lhs.bonds().contains(*id)
                    || unordered_pair(lhs.bonds().get(*id).expect("checked lhs bond").atom_ids())
                        == unordered_pair(*atoms)
                {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            BondDelta::ModifyField { id, .. } | BondDelta::ModifyConstraint { id, .. } => {
                self.require_available(lhs, created, Entity::Bond(*id))
            }
        }
    }

    fn validate_dative(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &DativeBondDelta,
    ) -> Result<(), ReactionIntegrityError> {
        match delta {
            DativeBondDelta::Add {
                donors, acceptor, ..
            } => self.require_atoms(lhs, created, donors.iter().copied().chain([*acceptor])),
            DativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                ..
            } => {
                let entity = Entity::DativeBond(*id);
                self.require_available(lhs, created, entity)?;
                self.require_atoms(lhs, created, donors.iter().copied().chain([*acceptor]))?;
                if !lhs.dative_bonds().contains(*id) || {
                    let view = lhs
                        .dative_bonds()
                        .get(*id)
                        .expect("checked lhs dative bond");
                    view.acceptor_id() == *acceptor
                        && unordered_ids(view.donor_ids()) == unordered_ids(donors.iter().copied())
                } {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            DativeBondDelta::ModifyField { id, .. }
            | DativeBondDelta::ModifyConstraint { id, .. } => {
                self.require_available(lhs, created, Entity::DativeBond(*id))
            }
        }
    }

    fn validate_aromatic(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &AromaticSystemDelta,
    ) -> Result<(), ReactionIntegrityError> {
        match delta {
            AromaticSystemDelta::Add { atoms, .. } => {
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            AromaticSystemDelta::Remove { id, atoms, .. } => {
                let entity = Entity::AromaticSystem(*id);
                self.require_available(lhs, created, entity)?;
                self.require_atoms(lhs, created, atoms.iter().copied())?;
                if !lhs.aromatic_systems().contains(*id)
                    || unordered_ids(
                        lhs.aromatic_systems()
                            .get(*id)
                            .expect("checked lhs aromatic system")
                            .atom_ids(),
                    ) == unordered_ids(atoms.iter().copied())
                {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            AromaticSystemDelta::ModifyField { id, .. }
            | AromaticSystemDelta::ModifyConstraint { id, .. } => {
                self.require_available(lhs, created, Entity::AromaticSystem(*id))
            }
        }
    }

    fn validate_multicenter(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &MulticenterBondDelta,
    ) -> Result<(), ReactionIntegrityError> {
        match delta {
            MulticenterBondDelta::Add { atoms, .. } => {
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            MulticenterBondDelta::Remove { id, atoms, .. } => {
                let entity = Entity::MulticenterBond(*id);
                self.require_available(lhs, created, entity)?;
                self.require_atoms(lhs, created, atoms.iter().copied())?;
                if !lhs.multicenter_bonds().contains(*id)
                    || unordered_ids(
                        lhs.multicenter_bonds()
                            .get(*id)
                            .expect("checked lhs multicenter bond")
                            .atom_ids(),
                    ) == unordered_ids(atoms.iter().copied())
                {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            MulticenterBondDelta::ModifyField { id, .. }
            | MulticenterBondDelta::ModifyConstraint { id, .. } => {
                self.require_available(lhs, created, Entity::MulticenterBond(*id))
            }
        }
    }

    fn validate_noncovalent(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &NoncovalentBondDelta,
    ) -> Result<(), ReactionIntegrityError> {
        match delta {
            NoncovalentBondDelta::Add { atoms, .. } => self.require_atoms(lhs, created, *atoms),
            NoncovalentBondDelta::Remove { id, atoms, .. } => {
                let entity = Entity::NoncovalentBond(*id);
                self.require_available(lhs, created, entity)?;
                self.require_atoms(lhs, created, *atoms)?;
                if !lhs.noncovalent_bonds().contains(*id)
                    || unordered_pair(
                        lhs.noncovalent_bonds()
                            .get(*id)
                            .expect("checked lhs noncovalent bond")
                            .atom_ids(),
                    ) == unordered_pair(*atoms)
                {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            NoncovalentBondDelta::ModifyField { id, .. }
            | NoncovalentBondDelta::ModifyConstraint { id, .. } => {
                self.require_available(lhs, created, Entity::NoncovalentBond(*id))
            }
        }
    }

    fn validate_stereo_atom(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &StereoAtomDelta,
    ) -> Result<(), ReactionIntegrityError> {
        match delta {
            StereoAtomDelta::Add { site, ligands, .. } => self.require_atoms(
                lhs,
                created,
                iter::once(*site).chain(ligands.iter().map(|l| l.atom_id)),
            ),
            StereoAtomDelta::Remove {
                id, site, ligands, ..
            } => {
                let entity = Entity::StereoAtom(*id);
                self.require_available(lhs, created, entity)?;
                self.require_atoms(
                    lhs,
                    created,
                    iter::once(*site).chain(ligands.iter().map(|l| l.atom_id)),
                )?;
                if !lhs.stereo_atoms().contains(*id) || {
                    let view = lhs
                        .stereo_atoms()
                        .get(*id)
                        .expect("checked lhs stereo atom");
                    let actual: Vec<StereoLigand> = view
                        .ligands()
                        .map(|ligand| StereoLigand::new(ligand.atom_id(), ligand.kind()))
                        .collect();
                    view.site_id() == *site && actual == *ligands
                } {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            StereoAtomDelta::ModifyField { id, .. }
            | StereoAtomDelta::ModifyConstraint { id, .. }
            | StereoAtomDelta::Apply { id, .. }
            | StereoAtomDelta::Swap { id, .. }
            | StereoAtomDelta::Mirror { id, .. } => {
                self.require_available(lhs, created, Entity::StereoAtom(*id))
            }
        }
    }

    fn validate_stereo_bond(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &StereoBondDelta,
    ) -> Result<(), ReactionIntegrityError> {
        match delta {
            StereoBondDelta::Add { site, ligands, .. } => {
                self.require_available(lhs, created, Entity::Bond(*site))?;
                self.require_atoms(lhs, created, ligands.iter().map(|l| l.atom_id))
            }
            StereoBondDelta::Remove {
                id, site, ligands, ..
            } => {
                let entity = Entity::StereoBond(*id);
                self.require_available(lhs, created, entity)?;
                self.require_available(lhs, created, Entity::Bond(*site))?;
                self.require_atoms(lhs, created, ligands.iter().map(|l| l.atom_id))?;
                if !lhs.stereo_bonds().contains(*id) || {
                    let view = lhs
                        .stereo_bonds()
                        .get(*id)
                        .expect("checked lhs stereo bond");
                    let actual: Vec<StereoLigand> = view
                        .ligands()
                        .map(|ligand| StereoLigand::new(ligand.atom_id(), ligand.kind()))
                        .collect();
                    view.site_id() == *site && actual == *ligands
                } {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            StereoBondDelta::ModifyField { id, .. }
            | StereoBondDelta::ModifyConstraint { id, .. }
            | StereoBondDelta::Apply { id, .. }
            | StereoBondDelta::Swap { id, .. }
            | StereoBondDelta::Mirror { id, .. } => {
                self.require_available(lhs, created, Entity::StereoBond(*id))
            }
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
            Constraint::StereoAtom(id, _, _) => {
                self.require_available(lhs, created, Entity::StereoAtom(*id))
            }
            Constraint::StereoBond(id, _, _) => {
                self.require_available(lhs, created, Entity::StereoBond(*id))
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
    /// The check covers the lhs molecule, delta references, created-id uniqueness, and the
    /// recorded incidence of removed entities. It does not impose DPO or chemistry semantics.
    pub fn check_integrity(&self) -> Result<(), ReactionIntegrityError> {
        self.lhs.check_integrity()?;
        ReactionIntegrityCheck.check(&self.lhs, &self.deltas)
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

fn unordered_pair(mut atoms: [AtomId; 2]) -> [AtomId; 2] {
    atoms.sort_unstable();
    atoms
}
