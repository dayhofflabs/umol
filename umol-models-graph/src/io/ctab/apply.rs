//! Apply property entries to molecule

use super::properties::{
    AtomAliasEntry, AtomValueEntry, ChargeEntry, IsotopeEntry, RadicalEntry, SGroupAtomListEntry,
    SGroupBondListEntry, SGroupLabelEntry, SGroupTypeEntry,
};
use crate::molecule::Molecule;
use crate::sgroup::{SGroup, SGroupType};
use umol::error::{DataError, ValidationError};
use umol::Result;

/// Trait for applying property entries to molecule
pub trait Apply {
    fn apply(self, molecule: &mut Molecule) -> Result<()>;
}

impl Apply for Vec<ChargeEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            if entry.atom_index >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(entry.atom_index).into());
            }

            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                // Check for conflicts - if atom already has a charge set
                if atom.charge != 0 && atom.charge != entry.charge {
                    return Err(ValidationError::InvalidComponent(format!(
                        "Charge conflict for atom {}: existing {} vs new {}",
                        entry.atom_index, atom.charge, entry.charge
                    ))
                    .into());
                }
                atom.charge = entry.charge;
            }
        }
        Ok(())
    }
}

impl Apply for Vec<RadicalEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            if entry.atom_index >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(entry.atom_index).into());
            }

            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                // Validate radical type
                let radical_value = match entry.radical_type {
                    0 => None,    // No radical
                    1 => Some(1), // Singlet (:)
                    2 => Some(2), // Doublet (. or ^)
                    3 => Some(3), // Triplet (^^)
                    _ => {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Invalid radical type for atom {}: {}",
                            entry.atom_index, entry.radical_type
                        ))
                        .into())
                    }
                };

                // Check for conflicts - only conflict if both are Some and different
                if let (Some(existing), Some(new)) = (atom.radical, radical_value) {
                    if existing != new {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Radical conflict for atom {}: existing Some({}) vs new Some({})",
                            entry.atom_index, existing, new
                        ))
                        .into());
                    }
                }
                atom.radical = radical_value;
            }
        }
        Ok(())
    }
}

impl Apply for Vec<IsotopeEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            if entry.atom_index >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(entry.atom_index).into());
            }

            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                // Check for conflicts
                if let Some(existing) = atom.isotope_mass {
                    if existing != entry.mass {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Isotope conflict for atom {}: existing {} vs new {}",
                            entry.atom_index, existing, entry.mass
                        ))
                        .into());
                    }
                }
                atom.isotope_mass = Some(entry.mass);
            }
        }
        Ok(())
    }
}

impl Apply for Vec<SGroupTypeEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            // Parse SGroup type
            let sgroup_type = SGroup::get_type(&entry.sgroup_type)?;

            // Ensure SGroup exists
            ensure_sgroup(molecule, entry.sgroup_index)?;

            // Apply the type
            if let Some(sgroup) = molecule.sgroups.get_mut(entry.sgroup_index) {
                // Check for conflicts
                if sgroup.group_type != SGroupType::Generic && sgroup.group_type != sgroup_type {
                    return Err(ValidationError::InvalidComponent(format!(
                        "SGroup type conflict for SGroup {}: existing {:?} vs new {:?}",
                        entry.sgroup_index, sgroup.group_type, sgroup_type
                    ))
                    .into());
                }
                sgroup.group_type = sgroup_type;
            }
        }
        Ok(())
    }
}

impl Apply for Vec<SGroupLabelEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            // Ensure SGroup exists
            ensure_sgroup(molecule, entry.sgroup_index)?;

            // Apply the label
            if let Some(sgroup) = molecule.sgroups.get_mut(entry.sgroup_index) {
                // Check for conflicts
                if let Some(ref existing_label) = sgroup.label {
                    if existing_label != &entry.label {
                        return Err(ValidationError::InvalidComponent(format!(
                            "SGroup label conflict for SGroup {}: existing '{}' vs new '{}'",
                            entry.sgroup_index, existing_label, entry.label
                        ))
                        .into());
                    }
                }
                sgroup.label = Some(entry.label);
            }
        }
        Ok(())
    }
}

impl Apply for SGroupAtomListEntry {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        // Validate all atom indices exist
        for &atom_index in &self.atom_indices {
            if atom_index >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(atom_index).into());
            }
        }

        // Ensure SGroup exists
        ensure_sgroup(molecule, self.sgroup_index)?;

        // Apply the atom list
        if let Some(sgroup) = molecule.sgroups.get_mut(self.sgroup_index) {
            // Check if atoms are already assigned
            if !sgroup.atom_indices.is_empty() && sgroup.atom_indices != self.atom_indices {
                return Err(ValidationError::InvalidComponent(format!(
                    "SGroup atom list conflict for SGroup {}",
                    self.sgroup_index
                ))
                .into());
            }
            sgroup.atom_indices = self.atom_indices;
        }
        Ok(())
    }
}

impl Apply for SGroupBondListEntry {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        // Validate all bond indices exist
        for &bond_index in &self.bond_indices {
            if bond_index >= molecule.bond_count() {
                return Err(DataError::MissingBondIndex(bond_index).into());
            }
        }

        // Ensure SGroup exists
        ensure_sgroup(molecule, self.sgroup_index)?;

        // Apply the bond list
        if let Some(sgroup) = molecule.sgroups.get_mut(self.sgroup_index) {
            // Check for conflicts
            if !sgroup.bond_indices.is_empty() && sgroup.bond_indices != self.bond_indices {
                return Err(ValidationError::InvalidComponent(format!(
                    "SGroup bond list conflict for SGroup {}",
                    self.sgroup_index
                ))
                .into());
            }
            sgroup.bond_indices = self.bond_indices;
        }
        Ok(())
    }
}

impl Apply for AtomAliasEntry {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        if self.atom_index >= molecule.atom_count() {
            return Err(DataError::MissingAtomIndex(self.atom_index).into());
        }

        if let Some(atom) = molecule.atom_mut(self.atom_index) {
            // Check for conflicts
            if let Some(existing_alias) = atom.properties.get("molFileAlias") {
                if existing_alias != &self.alias {
                    return Err(ValidationError::InvalidComponent(format!(
                        "Atom alias conflict for atom {}: existing '{}' vs new '{}'",
                        self.atom_index, existing_alias, self.alias
                    ))
                    .into());
                }
            }
            atom.properties
                .insert("molFileAlias".to_string(), self.alias);
        }
        Ok(())
    }
}

impl Apply for AtomValueEntry {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        if self.atom_index >= molecule.atom_count() {
            return Err(DataError::MissingAtomIndex(self.atom_index).into());
        }

        if let Some(atom) = molecule.atom_mut(self.atom_index) {
            // Check for conflicts
            if let Some(existing_value) = atom.properties.get("molFileValue") {
                if existing_value != &self.value {
                    return Err(ValidationError::InvalidComponent(format!(
                        "Atom value conflict for atom {}: existing '{}' vs new '{}'",
                        self.atom_index, existing_value, self.value
                    ))
                    .into());
                }
            }
            atom.properties
                .insert("molFileValue".to_string(), self.value);
        }
        Ok(())
    }
}

/// Ensure an SGroup with given index exists
fn ensure_sgroup(molecule: &mut Molecule, sgroup_index: usize) -> Result<()> {
    // Extend the vector to include the required index
    while molecule.sgroups.len() <= sgroup_index {
        let new_id = molecule.sgroups.len();
        molecule
            .sgroups
            .push(SGroup::new(new_id, SGroupType::Generic));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
