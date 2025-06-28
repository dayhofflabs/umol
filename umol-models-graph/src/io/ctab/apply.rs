//! Apply property entries to molecule

use super::properties::{
    AtomAliasEntry, AtomAttachmentOrderEntry, AtomListEntry, AtomValueEntry, AttachmentPointEntry,
    ChargeEntry, IsotopeEntry, LinkAtomEntry, RadicalEntry, RingBondCountEntry, SGroupAtomListEntry,
    SGroupBondListEntry, SGroupLabelEntry, SGroupTypeEntry, SubstitutionCountEntry, UnsaturatedAtomEntry,
};
use crate::atom::{AtomList, AtomSymbol, LinkAtomSpec};
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

impl Apply for AtomListEntry {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        if self.atom_index >= molecule.atom_count() {
            return Err(DataError::MissingAtomIndex(self.atom_index).into());
        }

        if let Some(atom) = molecule.atom_mut(self.atom_index) {
            // Create AtomList from elements
            let atom_list = AtomList {
                elements: self.elements,
            };

            // Check for conflicts - if atom already has a symbol set
            match &atom.symbol {
                AtomSymbol::Element(_) => {
                    // Replace element with atom list
                    atom.symbol = AtomSymbol::AtomList(atom_list);
                }
                AtomSymbol::AtomList(existing) => {
                    // Check for conflict
                    if existing.elements != atom_list.elements {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Atom list conflict for atom {}: existing vs new atom list",
                            self.atom_index
                        ))
                        .into());
                    }
                }
                _ => {
                    // Replace other symbols with atom list
                    atom.symbol = AtomSymbol::AtomList(atom_list);
                }
            }
        }
        Ok(())
    }
}

impl Apply for Vec<AttachmentPointEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            if entry.atom_index >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(entry.atom_index).into());
            }

            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                // Check for conflicts
                if let Some(existing) = atom.attachment_point {
                    if existing != entry.attachment_type {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Attachment point conflict for atom {}: existing {} vs new {}",
                            entry.atom_index, existing, entry.attachment_type
                        ))
                        .into());
                    }
                }
                atom.attachment_point = Some(entry.attachment_type);
            }
        }
        Ok(())
    }
}

impl Apply for AtomAttachmentOrderEntry {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        if self.atom_index >= molecule.atom_count() {
            return Err(DataError::MissingAtomIndex(self.atom_index).into());
        }

        if let Some(atom) = molecule.atom_mut(self.atom_index) {
            // Check for conflicts
            if let Some(ref existing) = atom.attachment_order {
                if existing != &self.attachments {
                    return Err(ValidationError::InvalidComponent(format!(
                        "Attachment order conflict for atom {}: existing vs new order",
                        self.atom_index
                    ))
                    .into());
                }
            }
            atom.attachment_order = Some(self.attachments);
        }
        Ok(())
    }
}

impl Apply for Vec<RingBondCountEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            if entry.atom_index >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(entry.atom_index).into());
            }

            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                // Check for conflicts
                if let Some(existing) = atom.ring_bond_count {
                    if existing != entry.ring_bond_count {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Ring bond count conflict for atom {}: existing {} vs new {}",
                            entry.atom_index, existing, entry.ring_bond_count
                        ))
                        .into());
                    }
                }
                atom.ring_bond_count = Some(entry.ring_bond_count);
            }
        }
        Ok(())
    }
}

impl Apply for Vec<SubstitutionCountEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            if entry.atom_index >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(entry.atom_index).into());
            }

            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                // Check for conflicts
                if let Some(existing) = atom.substitution_count {
                    if existing != entry.substitution_count {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Substitution count conflict for atom {}: existing {} vs new {}",
                            entry.atom_index, existing, entry.substitution_count
                        ))
                        .into());
                    }
                }
                atom.substitution_count = Some(entry.substitution_count);
            }
        }
        Ok(())
    }
}

impl Apply for Vec<UnsaturatedAtomEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            if entry.atom_index >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(entry.atom_index).into());
            }

            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                // Convert integer to boolean
                let unsaturated_value = match entry.unsaturated {
                    0 => None,    // Off
                    1 => Some(true), // On
                    _ => {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Invalid unsaturated value for atom {}: {}",
                            entry.atom_index, entry.unsaturated
                        ))
                        .into())
                    }
                };

                // Check for conflicts
                if let Some(existing) = atom.unsaturated {
                    if let Some(new) = unsaturated_value {
                        if existing != new {
                            return Err(ValidationError::InvalidComponent(format!(
                                "Unsaturated conflict for atom {}: existing {} vs new {}",
                                entry.atom_index, existing, new
                            ))
                            .into());
                        }
                    }
                }
                atom.unsaturated = unsaturated_value;
            }
        }
        Ok(())
    }
}

impl Apply for Vec<LinkAtomEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            if entry.atom_index >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(entry.atom_index).into());
            }

            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                let link_spec = LinkAtomSpec {
                    repeat_count: entry.repeat_count,
                    bond1: entry.bond1,
                    bond2: entry.bond2,
                };

                // Check for conflicts
                if let Some(ref existing) = atom.link_atom {
                    if *existing != link_spec {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Link atom conflict for atom {}: existing vs new link specification",
                            entry.atom_index
                        ))
                        .into());
                    }
                }
                atom.link_atom = Some(link_spec);
            }
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
