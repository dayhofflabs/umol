//! Checked mutation of molecule-level constraints.

use std::slice::Iter;

use super::integrity::{check_constraint_references, check_molecule_constraint};
use super::{Constraint, Constraints, Entity, Molecule, MoleculeIntegrityError};

/// Mutable access to a molecule's constraints, preserving representation integrity.
///
/// Obtained through [`Molecule::constraints_mut`]. Incoming constraints are checked against
/// the molecule's entity ids and stereo frames before storage changes. Constraint satisfaction
/// is not checked. The view does not clone the molecule or check unrelated storage.
///
/// # Semantic properties
///
/// - Successful writes preserve the supplied order and duplicate entries exactly.
/// - A rejected write leaves the entire molecule unchanged.
/// - Removing constraints preserves representation integrity.
///
/// These properties are exercised by exact cases and generated constraint sequences compared
/// with publication through the checked editor boundary.
pub struct ConstraintsViewMut<'a> {
    molecule: &'a mut Molecule,
}

impl Molecule {
    /// Borrow molecule-level constraints for checked mutation.
    ///
    /// Insertions and replacement check incoming references and stereo frames. No molecule
    /// copy or finalization is required; each successful write preserves integrity.
    pub fn constraints_mut(&mut self) -> ConstraintsViewMut<'_> {
        ConstraintsViewMut { molecule: self }
    }
}

impl ConstraintsViewMut<'_> {
    /// Number of stored constraints, including duplicate entries.
    pub fn len(&self) -> usize {
        self.molecule.constraints.len()
    }

    /// Whether the constraint collection is empty.
    pub fn is_empty(&self) -> bool {
        self.molecule.constraints.is_empty()
    }

    /// Borrow the constraints in stored order.
    pub fn as_slice(&self) -> &[Constraint] {
        self.molecule.constraints.as_slice()
    }

    /// Iterate over constraints in stored order.
    pub fn iter(&self) -> Iter<'_, Constraint> {
        self.molecule.constraints.iter()
    }

    /// Append one constraint after checking its references and stereo frames.
    ///
    /// # Errors
    ///
    /// Returns the incoming constraint's integrity error without changing the molecule.
    pub fn push(&mut self, constraint: Constraint) -> Result<(), MoleculeIntegrityError> {
        self.check_constraint(&constraint)?;
        self.molecule.constraints.push(constraint);
        Ok(())
    }

    /// Append constraints in their supplied order, preserving duplicates.
    ///
    /// All incoming constraints are checked before any are appended.
    ///
    /// # Errors
    ///
    /// Returns the first incoming constraint integrity error without changing the molecule.
    pub fn extend(&mut self, constraints: Constraints) -> Result<(), MoleculeIntegrityError> {
        for constraint in constraints.iter() {
            self.check_constraint(constraint)?;
        }
        for constraint in constraints {
            self.molecule.constraints.push(constraint);
        }
        Ok(())
    }

    /// Replace the entire collection, preserving the supplied order and duplicates.
    ///
    /// # Errors
    ///
    /// Returns the first incoming constraint integrity error without changing the molecule.
    pub fn replace(&mut self, constraints: Constraints) -> Result<(), MoleculeIntegrityError> {
        for constraint in constraints.iter() {
            self.check_constraint(constraint)?;
        }
        self.molecule.constraints = constraints;
        Ok(())
    }

    /// Remove and return the constraint at a position, preserving the remaining order.
    ///
    /// # Panics
    ///
    /// Panics if `position >= self.len()`.
    pub fn remove_at(&mut self, position: usize) -> Constraint {
        self.molecule.constraints.remove_at(position)
    }

    /// Remove every constraint.
    pub fn clear(&mut self) {
        self.molecule.constraints.clear();
    }

    fn check_constraint(&self, constraint: &Constraint) -> Result<(), MoleculeIntegrityError> {
        let molecule = &self.molecule;
        check_constraint_references(constraint, &|entity| match entity {
            Entity::Atom(id) => id.index() < molecule.atoms.len(),
            Entity::Bond(id) => id.index() < molecule.bonds.len(),
            Entity::DativeBond(id) => id.index() < molecule.dative_bonds.count(),
            Entity::AromaticSystem(id) => id.index() < molecule.aromatic_systems.count(),
            Entity::MulticenterBond(id) => id.index() < molecule.multicenter_bonds.count(),
            Entity::NoncovalentBond(id) => id.index() < molecule.noncovalent_bonds.count(),
            Entity::StereoAtom(id) => id.index() < molecule.stereo_atoms.count(),
            Entity::StereoBond(id) => id.index() < molecule.stereo_bonds.count(),
        })
        .map_err(|entity| MoleculeIntegrityError::InvalidReference { entity })?;
        check_molecule_constraint(molecule, constraint)
    }
}

#[cfg(test)]
mod tests;
