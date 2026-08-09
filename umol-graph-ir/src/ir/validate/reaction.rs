//! Reaction integrity validation for entity references and structural incidence.

use std::collections::HashSet;
use std::iter;

use thiserror::Error;
use umol_utils::solution::Solution;

use super::super::constraint::{Constraint, MoleculeConstraint, RelationalConstraint};
use super::super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, ConstraintDelta, DativeBondDelta, Delta, Deltas,
    MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta,
};
use super::super::entity::Entity;
use super::super::id::AtomId;
use super::super::ligand::StereoLigand;
use super::super::molecule::Molecule;

/// Tier-1 integrity validator for reaction delta references and structural incidence.
#[derive(Clone, Copy, Debug, Default)]
pub struct ReactionIntegrityValidator;

#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionIntegrityContradiction {
    #[error("reaction references unavailable entity {entity:?}")]
    InvalidReference { entity: Entity },
    #[error("reaction incidence does not match lhs entity {entity:?}")]
    IncidenceMismatch { entity: Entity },
}

#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionIntegrityError {}

impl ReactionIntegrityValidator {
    pub fn validate(
        &self,
        lhs: &Molecule,
        deltas: &Deltas,
    ) -> Result<Solution<(), ReactionIntegrityContradiction>, ReactionIntegrityError> {
        match self.validate_inner(lhs, deltas) {
            Ok(()) => Ok(Solution::Determined(())),
            Err(contradiction) => Ok(Solution::Contradictory(contradiction)),
        }
    }

    fn validate_inner(
        &self,
        lhs: &Molecule,
        deltas: &Deltas,
    ) -> Result<(), ReactionIntegrityContradiction> {
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
                    return Err(ReactionIntegrityContradiction::InvalidReference { entity });
                }
            }
        }
        for delta in deltas.iter() {
            match delta {
                Delta::Atom(delta) => self.validate_atom(lhs, delta)?,
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

    fn require_lhs(
        &self,
        lhs: &Molecule,
        entity: Entity,
    ) -> Result<(), ReactionIntegrityContradiction> {
        if contains_entity(lhs, entity) {
            Ok(())
        } else {
            Err(ReactionIntegrityContradiction::InvalidReference { entity })
        }
    }

    fn require_available(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        entity: Entity,
    ) -> Result<(), ReactionIntegrityContradiction> {
        if contains_entity(lhs, entity) || created.contains(&entity) {
            Ok(())
        } else {
            Err(ReactionIntegrityContradiction::InvalidReference { entity })
        }
    }

    fn require_atoms(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Result<(), ReactionIntegrityContradiction> {
        for atom in atoms {
            self.require_available(lhs, created, Entity::Atom(atom))?;
        }
        Ok(())
    }

    fn incidence_mismatch(entity: Entity) -> ReactionIntegrityContradiction {
        ReactionIntegrityContradiction::IncidenceMismatch { entity }
    }

    fn validate_atom(
        &self,
        lhs: &Molecule,
        delta: &AtomDelta,
    ) -> Result<(), ReactionIntegrityContradiction> {
        match delta {
            AtomDelta::Add { .. } => Ok(()),
            AtomDelta::Remove { id, .. }
            | AtomDelta::ModifyField { id, .. }
            | AtomDelta::ModifyConstraint { id, .. } => self.require_lhs(lhs, Entity::Atom(*id)),
        }
    }

    fn validate_bond(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &BondDelta,
    ) -> Result<(), ReactionIntegrityContradiction> {
        match delta {
            BondDelta::Add { atoms, .. } => self.require_atoms(lhs, created, *atoms),
            BondDelta::Remove { id, atoms, .. } => {
                let entity = Entity::Bond(*id);
                self.require_lhs(lhs, entity)?;
                self.require_atoms(lhs, created, *atoms)?;
                let actual = lhs
                    .bonds()
                    .get(*id)
                    .ok_or(ReactionIntegrityContradiction::InvalidReference { entity })?
                    .atom_ids();
                if unordered_pair(actual) == unordered_pair(*atoms) {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            BondDelta::ModifyField { id, .. } | BondDelta::ModifyConstraint { id, .. } => {
                self.require_lhs(lhs, Entity::Bond(*id))
            }
        }
    }

    fn validate_dative(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &DativeBondDelta,
    ) -> Result<(), ReactionIntegrityContradiction> {
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
                self.require_lhs(lhs, entity)?;
                self.require_atoms(lhs, created, donors.iter().copied().chain([*acceptor]))?;
                let view = lhs
                    .dative_bonds()
                    .get(*id)
                    .ok_or(ReactionIntegrityContradiction::InvalidReference { entity })?;
                if view.acceptor_id() == *acceptor
                    && unordered_ids(view.donor_ids()) == unordered_ids(donors.iter().copied())
                {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            DativeBondDelta::ModifyField { id, .. }
            | DativeBondDelta::ModifyConstraint { id, .. } => {
                self.require_lhs(lhs, Entity::DativeBond(*id))
            }
        }
    }

    fn validate_aromatic(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &AromaticSystemDelta,
    ) -> Result<(), ReactionIntegrityContradiction> {
        match delta {
            AromaticSystemDelta::Add { atoms, .. } => {
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            AromaticSystemDelta::Remove { id, atoms, .. } => {
                let entity = Entity::AromaticSystem(*id);
                self.require_lhs(lhs, entity)?;
                self.require_atoms(lhs, created, atoms.iter().copied())?;
                let view = lhs
                    .aromatic_systems()
                    .get(*id)
                    .ok_or(ReactionIntegrityContradiction::InvalidReference { entity })?;
                if unordered_ids(view.atom_ids()) == unordered_ids(atoms.iter().copied()) {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            AromaticSystemDelta::ModifyField { id, .. }
            | AromaticSystemDelta::ModifyConstraint { id, .. } => {
                self.require_lhs(lhs, Entity::AromaticSystem(*id))
            }
        }
    }

    fn validate_multicenter(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &MulticenterBondDelta,
    ) -> Result<(), ReactionIntegrityContradiction> {
        match delta {
            MulticenterBondDelta::Add { atoms, .. } => {
                self.require_atoms(lhs, created, atoms.iter().copied())
            }
            MulticenterBondDelta::Remove { id, atoms, .. } => {
                let entity = Entity::MulticenterBond(*id);
                self.require_lhs(lhs, entity)?;
                self.require_atoms(lhs, created, atoms.iter().copied())?;
                let view = lhs
                    .multicenter_bonds()
                    .get(*id)
                    .ok_or(ReactionIntegrityContradiction::InvalidReference { entity })?;
                if unordered_ids(view.atom_ids()) == unordered_ids(atoms.iter().copied()) {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            MulticenterBondDelta::ModifyField { id, .. }
            | MulticenterBondDelta::ModifyConstraint { id, .. } => {
                self.require_lhs(lhs, Entity::MulticenterBond(*id))
            }
        }
    }

    fn validate_noncovalent(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &NoncovalentBondDelta,
    ) -> Result<(), ReactionIntegrityContradiction> {
        match delta {
            NoncovalentBondDelta::Add { atoms, .. } => self.require_atoms(lhs, created, *atoms),
            NoncovalentBondDelta::Remove { id, atoms, .. } => {
                let entity = Entity::NoncovalentBond(*id);
                self.require_lhs(lhs, entity)?;
                self.require_atoms(lhs, created, *atoms)?;
                let actual = lhs
                    .noncovalent_bonds()
                    .get(*id)
                    .ok_or(ReactionIntegrityContradiction::InvalidReference { entity })?
                    .atom_ids();
                if unordered_pair(actual) == unordered_pair(*atoms) {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            NoncovalentBondDelta::ModifyField { id, .. }
            | NoncovalentBondDelta::ModifyConstraint { id, .. } => {
                self.require_lhs(lhs, Entity::NoncovalentBond(*id))
            }
        }
    }

    fn validate_stereo_atom(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &StereoAtomDelta,
    ) -> Result<(), ReactionIntegrityContradiction> {
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
                self.require_lhs(lhs, entity)?;
                self.require_atoms(
                    lhs,
                    created,
                    iter::once(*site).chain(ligands.iter().map(|l| l.atom_id)),
                )?;
                let view = lhs
                    .stereo_atoms()
                    .get(*id)
                    .ok_or(ReactionIntegrityContradiction::InvalidReference { entity })?;
                let actual: Vec<StereoLigand> = view
                    .ligands()
                    .map(|ligand| StereoLigand::new(ligand.atom_id(), ligand.kind()))
                    .collect();
                if view.site_id() == *site && actual == *ligands {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            StereoAtomDelta::ModifyField { id, .. }
            | StereoAtomDelta::ModifyConstraint { id, .. }
            | StereoAtomDelta::Apply { id, .. }
            | StereoAtomDelta::Swap { id, .. }
            | StereoAtomDelta::Mirror { id, .. } => self.require_lhs(lhs, Entity::StereoAtom(*id)),
        }
    }

    fn validate_stereo_bond(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        delta: &StereoBondDelta,
    ) -> Result<(), ReactionIntegrityContradiction> {
        match delta {
            StereoBondDelta::Add { site, ligands, .. } => {
                self.require_available(lhs, created, Entity::Bond(*site))?;
                self.require_atoms(lhs, created, ligands.iter().map(|l| l.atom_id))
            }
            StereoBondDelta::Remove {
                id, site, ligands, ..
            } => {
                let entity = Entity::StereoBond(*id);
                self.require_lhs(lhs, entity)?;
                self.require_lhs(lhs, Entity::Bond(*site))?;
                self.require_atoms(lhs, created, ligands.iter().map(|l| l.atom_id))?;
                let view = lhs
                    .stereo_bonds()
                    .get(*id)
                    .ok_or(ReactionIntegrityContradiction::InvalidReference { entity })?;
                let actual: Vec<StereoLigand> = view
                    .ligands()
                    .map(|ligand| StereoLigand::new(ligand.atom_id(), ligand.kind()))
                    .collect();
                if view.site_id() == *site && actual == *ligands {
                    Ok(())
                } else {
                    Err(Self::incidence_mismatch(entity))
                }
            }
            StereoBondDelta::ModifyField { id, .. }
            | StereoBondDelta::ModifyConstraint { id, .. }
            | StereoBondDelta::Apply { id, .. }
            | StereoBondDelta::Swap { id, .. }
            | StereoBondDelta::Mirror { id, .. } => self.require_lhs(lhs, Entity::StereoBond(*id)),
        }
    }

    fn validate_constraint(
        &self,
        lhs: &Molecule,
        created: &HashSet<Entity>,
        constraint: &Constraint,
    ) -> Result<(), ReactionIntegrityContradiction> {
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
    ) -> Result<(), ReactionIntegrityContradiction> {
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
    ) -> Result<(), ReactionIntegrityContradiction> {
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
            MoleculeConstraint::SubPattern { anchor, pattern } => {
                let pattern_created = HashSet::new();
                for &(target, pattern_id) in anchor.atoms() {
                    self.require_available(lhs, created, Entity::Atom(target))?;
                    self.require_lhs(pattern, Entity::Atom(pattern_id))?;
                }
                for &(target, pattern_id) in anchor.bonds() {
                    self.require_available(lhs, created, Entity::Bond(target))?;
                    self.require_lhs(pattern, Entity::Bond(pattern_id))?;
                }
                for &(target, pattern_id) in anchor.dative_bonds() {
                    self.require_available(lhs, created, Entity::DativeBond(target))?;
                    self.require_lhs(pattern, Entity::DativeBond(pattern_id))?;
                }
                for &(target, pattern_id) in anchor.aromatic_systems() {
                    self.require_available(lhs, created, Entity::AromaticSystem(target))?;
                    self.require_lhs(pattern, Entity::AromaticSystem(pattern_id))?;
                }
                for &(target, pattern_id) in anchor.multicenter_bonds() {
                    self.require_available(lhs, created, Entity::MulticenterBond(target))?;
                    self.require_lhs(pattern, Entity::MulticenterBond(pattern_id))?;
                }
                for &(target, pattern_id) in anchor.noncovalent_bonds() {
                    self.require_available(lhs, created, Entity::NoncovalentBond(target))?;
                    self.require_lhs(pattern, Entity::NoncovalentBond(pattern_id))?;
                }
                for &(target, pattern_id) in anchor.stereo_atoms() {
                    self.require_available(lhs, created, Entity::StereoAtom(target))?;
                    self.require_lhs(pattern, Entity::StereoAtom(pattern_id))?;
                }
                for &(target, pattern_id) in anchor.stereo_bonds() {
                    self.require_available(lhs, created, Entity::StereoBond(target))?;
                    self.require_lhs(pattern, Entity::StereoBond(pattern_id))?;
                }
                for constraint in pattern.constraints().iter() {
                    self.validate_constraint(pattern, &pattern_created, constraint)?;
                }
                Ok(())
            }
        }
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
