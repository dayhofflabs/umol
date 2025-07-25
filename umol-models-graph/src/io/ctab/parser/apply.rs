//! Apply property entries to molecule

use super::convert::{
    convert_atom_isotope_mass_number, convert_attachment_point_code, convert_radical_type_code,
    convert_ring_bond_count_code, convert_substitution_count_code, convert_unsaturated_atom_code,
};
use super::properties::{
    AtomAliasEntry, AtomAttachmentOrderEntry, AtomListEntry, AtomValueEntry, AttachmentPointEntry,
    ChargeEntry, IsotopeEntry, LinkAtomEntry, PropertyEntries, RGroupLabelEntry, RGroupLogicEntry,
    RadicalEntry, RingBondCountEntry, SGroupAtomListEntry, SGroupBondListEntry, SGroupLabelEntry,
    SGroupSubtypeEntry, SGroupTypeEntry, SubstitutionCountEntry, UnsaturatedAtomEntry,
};
use crate::io::ctab::atom::{AtomList, AtomSymbol, LinkAtom};
use crate::io::ctab::molecule::Molecule;
use crate::io::ctab::rgroup::RGroup;
use crate::io::ctab::sgroup::{SGroup, SGroupType};
use umol::error::{DataError, ValidationError};
use umol::{Error, Result};

/// Trait for applying property entries to molecule
pub trait Apply {
    fn apply(self, molecule: &mut Molecule) -> Result<()>;
}

/// Implementation for PropertyEntries enum - dispatches to specific implementations
impl Apply for PropertyEntries {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        match self {
            PropertyEntries::AtomAliasEntry(entry) => entry.apply(molecule),
            PropertyEntries::AtomValueEntry(entry) => entry.apply(molecule),
            PropertyEntries::ChargeEntries(entries) => entries.apply(molecule),
            PropertyEntries::RadicalEntries(entries) => entries.apply(molecule),
            PropertyEntries::IsotopeEntries(entries) => entries.apply(molecule),
            PropertyEntries::RingBondCountEntries(entries) => entries.apply(molecule),
            PropertyEntries::SubstitutionCountEntries(entries) => entries.apply(molecule),
            PropertyEntries::UnsaturatedAtomEntries(entries) => entries.apply(molecule),
            PropertyEntries::LinkAtomEntries(entries) => entries.apply(molecule),
            PropertyEntries::AtomListEntry(entry) => entry.apply(molecule),
            PropertyEntries::AttachmentPointEntries(entries) => entries.apply(molecule),
            PropertyEntries::AtomAttachmentOrderEntry(entry) => entry.apply(molecule),
            PropertyEntries::RGroupLabelEntries(entries) => entries.apply(molecule),
            PropertyEntries::RGroupLogicEntry(entry) => entry.apply(molecule),
            PropertyEntries::SGroupTypeEntries(entries) => entries.apply(molecule),
            PropertyEntries::SGroupSubtypeEntries(entries) => entries.apply(molecule),
            PropertyEntries::SGroupLabelEntries(entries) => entries.apply(molecule),
            // [SCN]
            // [SDS EXP]
            PropertyEntries::SGroupAtomListEntry(entry) => entry.apply(molecule),
            PropertyEntries::SGroupBondListEntry(entry) => entry.apply(molecule),
            PropertyEntries::End => Ok(()),
        }
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

impl Apply for Vec<ChargeEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            if entry.atom_index >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(entry.atom_index).into());
            }

            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                // Overwrite existing charge and radical values
                atom.charge = entry.charge;
                atom.radical = None;
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
                let radical_value = convert_radical_type_code(entry.radical_type)?;
                // Overwrite existing radical and charge values
                atom.radical = radical_value;
                atom.charge = 0;
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
                let element = match atom.symbol {
                    AtomSymbol::Element(element) => Ok::<_, Error>(element),
                    AtomSymbol::NamedIsotope(isotope) => Ok::<_, Error>(isotope.element()),
                    _ => Err(ValidationError::InvalidComponent(format!(
                        "Cannot set isotope for atom {}: {:?}",
                        entry.atom_index, atom.symbol
                    ))
                    .into()),
                }?;
                let mass_number = convert_atom_isotope_mass_number(element, entry.mass)?;
                if let Some(existing) = atom.isotope_mass {
                    if Some(existing) != mass_number {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Isotope conflict for atom {}: existing {:?} vs new {:?}",
                            entry.atom_index, existing, mass_number
                        ))
                        .into());
                    }
                }
                atom.isotope_mass = mass_number;
            }
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

            // Conflicting data is an error
            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                let ring_bond_count = convert_ring_bond_count_code(entry.ring_bond_count)?;
                if let Some(existing) = atom.ring_bond_count {
                    if Some(existing) != ring_bond_count {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Ring bond count conflict for atom {}: existing {:?} vs new {:?}",
                            entry.atom_index, existing, ring_bond_count
                        ))
                        .into());
                    }
                }
                atom.ring_bond_count = ring_bond_count;
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

            // Conflicting data is an error
            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                let substitution_count = convert_substitution_count_code(entry.substitution_count)?;
                if let Some(existing) = atom.substitution_count {
                    if Some(existing) != substitution_count {
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

            // Conflicting data is an error
            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                let unsaturated_value = convert_unsaturated_atom_code(entry.unsaturated)?;
                if let Some(existing) = atom.unsaturated {
                    if Some(existing) != unsaturated_value {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Unsaturated conflict for atom {}: existing {:?} vs new {:?}",
                            entry.atom_index, existing, unsaturated_value
                        ))
                        .into());
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
            if entry.subs_index1 >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(entry.subs_index1).into());
            }
            if entry.subs_index1 == entry.atom_index {
                return Err(ValidationError::InvalidComponent(format!(
                    "Invalid link atom: subs index 1 cannot be the same as the atom index",
                ))
                .into());
            }
            if let Some(subs_index2) = entry.subs_index2 {
                if subs_index2 >= molecule.atom_count() {
                    return Err(DataError::MissingAtomIndex(subs_index2).into());
                }
                if subs_index2 == entry.atom_index {
                    return Err(ValidationError::InvalidComponent(format!(
                        "Invalid link atom: subs index 2 cannot be the same as the atom index",
                    ))
                    .into());
                }
            }

            if entry.repeat_count < 2 {
                return Err(ValidationError::InvalidComponent(format!(
                    "Invalid link atom: repeat count must be >= 2, got {}",
                    entry.repeat_count
                ))
                .into());
            }

            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                let link_atom = LinkAtom {
                    repeat_count: entry.repeat_count,
                    subs_index1: entry.subs_index1,
                    subs_index2: entry.subs_index2,
                };

                // Check for conflicts
                if let Some(existing) = atom.link_atom {
                    if existing != link_atom {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Link atom conflict for atom {}: existing vs new link specification",
                            entry.atom_index
                        ))
                        .into());
                    }
                }
                atom.link_atom = Some(link_atom);
            }
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
            let atom_list = AtomList {
                elements: self.elements,
                exclusion: self.exclusion,
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

            if let Some(attachment_type) = convert_attachment_point_code(entry.attachment_type)? {
                if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                    // Check for conflicts
                    if let Some(existing) = atom.attachment_point {
                        if existing != attachment_type {
                            return Err(ValidationError::InvalidComponent(format!(
                                "Attachment point conflict for atom {}: existing {:?} vs new {:?}",
                                entry.atom_index, existing, attachment_type
                            ))
                            .into());
                        }
                    }
                    atom.attachment_point = Some(attachment_type);
                }
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

        if self.attachments.len() > 2 {
            return Err(ValidationError::InvalidComponent(format!(
                "Attachment order invalid for atom {}: more than 2 attachments",
                self.atom_index
            ))
            .into());
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

impl Apply for Vec<RGroupLabelEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            if entry.atom_index >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(entry.atom_index).into());
            }

            // Check for RGroup label conflicts
            for atom in molecule.atoms() {
                if let AtomSymbol::RGroup(ref rgroup) = atom.symbol {
                    if let Some(label) = rgroup.label {
                        if label == entry.label {
                            return Err(ValidationError::InvalidComponent(format!(
                                "RGroup label conflict: label '{}' is not unique",
                                entry.label
                            ))
                            .into());
                        }
                    }
                }
            }

            // Apply RGroup label
            if let Some(atom) = molecule.atom_mut(entry.atom_index) {
                match atom.symbol {
                    AtomSymbol::RGroup(ref mut rgroup) => {
                        if rgroup.label.is_some() && rgroup.label.unwrap() != entry.label {
                            return Err(ValidationError::InvalidComponent(format!(
                                "RGroup label conflict: existing '{}' vs new '{:?}'",
                                rgroup.label.unwrap(),
                                entry.label
                            ))
                            .into());
                        }
                        rgroup.label = Some(entry.label);
                        rgroup.explicit = true;
                    }
                    AtomSymbol::Element(_) | AtomSymbol::NamedIsotope(_) => {
                        atom.symbol = AtomSymbol::RGroup(RGroup::new(Some(entry.label)));
                    }
                    _ => {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Cannot set RGroup label for atom {} with symbol {:?}",
                            entry.atom_index, atom.symbol
                        ))
                        .into());
                    }
                }
            }
        }
        Ok(())
    }
}

impl Apply for RGroupLogicEntry {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        let mut index = 0;
        for atom_index in molecule.atom_indices() {
            if let Some(atom) = molecule.atom(atom_index) {
                if let AtomSymbol::RGroup(ref rgroup) = atom.symbol {
                    if let Some(label) = rgroup.label {
                        if label == self.label {
                            index = atom_index;
                            break;
                        }
                    }
                }
            }
        }
        if let Some(atom) = molecule.atom_mut(index) {
            if let AtomSymbol::RGroup(ref mut rgroup) = atom.symbol {
                rgroup.dependent_label = self.dependent_label;
                rgroup.rgroup_or_h = self.rgroup_or_h;
                rgroup.occurrence = self.occurrence;
            }
        }
        Ok(())
    }
}

impl Apply for Vec<SGroupTypeEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            let sgroup_index = entry.sgroup_index;
            make_sgroup(molecule, entry.sgroup_index)?;
            molecule.sgroups.get_mut(&sgroup_index).unwrap().group_type = entry.sgroup_type;
        }
        Ok(())
    }
}

impl Apply for Vec<SGroupSubtypeEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            let sgroup_index = entry.sgroup_index;
            ensure_sgroup(molecule, entry.sgroup_index)?;
            let sgroup = molecule.sgroups.get_mut(&sgroup_index).unwrap();

            if sgroup.group_subtype.is_some() && sgroup.group_subtype != Some(entry.sgroup_subtype)
            {
                return Err(ValidationError::InvalidComponent(format!(
                    "SGroup subtype conflict for SGroup {}: existing {:?} vs new {:?}",
                    sgroup_index, sgroup.group_subtype, entry.sgroup_subtype
                ))
                .into());
            }
            sgroup.group_subtype = Some(entry.sgroup_subtype);
        }
        Ok(())
    }
}

impl Apply for Vec<SGroupLabelEntry> {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        for entry in self {
            // Ensure SGroup exists
            ensure_sgroup(molecule, entry.sgroup_index)?;

            // Check if label is unique
            if molecule
                .sgroups
                .values()
                .any(|s| s.label == Some(entry.label))
            {
                return Err(ValidationError::InvalidComponent(format!(
                    "SGroup label conflict: duplicate label '{}'",
                    entry.label
                ))
                .into());
            }

            // Apply the label
            if let Some(sgroup) = molecule.sgroups.get_mut(&entry.sgroup_index) {
                // Having multiple SLB conflicting entries for the same SGroup is invalid
                if sgroup.label.is_some() && sgroup.label != Some(entry.label) {
                    return Err(ValidationError::InvalidComponent(format!(
                        "SGroup label conflict {}: existing '{}' vs new '{}'",
                        entry.sgroup_index,
                        sgroup.label.unwrap(),
                        entry.label
                    ))
                    .into());
                }
                sgroup.label = Some(entry.label);
            }
        }
        Ok(())
    }
}

// [SCN]
// [SDS EXP]

impl Apply for SGroupAtomListEntry {
    fn apply(self, molecule: &mut Molecule) -> Result<()> {
        // Ensure SGroup exists
        ensure_sgroup(molecule, self.sgroup_index)?;

        // Validate all atom indices exist
        for &atom_index in &self.atom_indices {
            if atom_index >= molecule.atom_count() {
                return Err(DataError::MissingAtomIndex(atom_index).into());
            }
        }

        // Apply the atom list
        if let Some(sgroup) = molecule.sgroups.get_mut(&self.sgroup_index) {
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
        // Ensure SGroup exists
        ensure_sgroup(molecule, self.sgroup_index)?;

        // Validate all bond indices exist
        for &bond_index in &self.bond_indices {
            if bond_index >= molecule.bond_count() {
                return Err(DataError::MissingBondIndex(bond_index).into());
            }
        }

        // Apply the bond list
        if let Some(sgroup) = molecule.sgroups.get_mut(&self.sgroup_index) {
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

// Utilities

/// Create a new SGroup with the given index if it doesn't exist
fn make_sgroup(molecule: &mut Molecule, sgroup_index: usize) -> Result<()> {
    if molecule.sgroups.contains_key(&sgroup_index) {
        return Err(ValidationError::InvalidComponent(format!(
            "SGroup index conflict: {} already exists",
            sgroup_index
        ))
        .into());
    }
    molecule
        .sgroups
        .insert(sgroup_index, SGroup::new(SGroupType::Generic));
    Ok(())
}

/// Ensure SGroup with given index exists
fn ensure_sgroup(molecule: &mut Molecule, sgroup_index: usize) -> Result<()> {
    if !molecule.sgroups.contains_key(&sgroup_index) {
        return Err(ValidationError::InvalidComponent(format!(
            "Invalid SGroup index: {}",
            sgroup_index
        ))
        .into());
    }
    Ok(())
}

#[cfg(test)]
mod tests;
