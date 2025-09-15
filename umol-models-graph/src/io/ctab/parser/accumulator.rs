//! Accumulator for molecular properties from a CTAB file

use std::collections::BTreeMap;

use umol::error::{DataError, ValidationError};
use umol::Result;
use umol_data::Element;

use super::context::Context;
use super::convert::{
    convert_atom_isotope_mass_number, convert_attachment_point_code, convert_bondlike_type_code,
    convert_radical_type_code, convert_ring_bond_count_code, convert_substitution_count_code,
    convert_unsaturated_atom_code,
};
use crate::io::config::ParseFlags;
use crate::io::ctab::atom::{
    AtomLike, AtomList, AtomRadical, AtomSymbol, AttachmentPointType, LinkAtom, RingBondCount,
    SubstitutionCount, UnsaturatedAtom,
};
use crate::io::ctab::molecule::{Molecule, MoleculeLike};
use crate::io::ctab::parser::properties::{PropertyEntries, SGroupDataEntry};
use crate::io::ctab::rgroup::{RGroup, RGroupOccurrence};
use crate::io::ctab::sgroup::{
    SGroup, SGroupBracketCoords, SGroupConnectingBond, SGroupConnectivity, SGroupData,
    SGroupDataDisplay, SGroupMultiplier, SGroupSubtype, SGroupType,
};

// Accumulator for properties of a single atom
#[derive(Debug, Default)]
pub struct AtomProperties {
    pub alias: Option<String>,
    pub value: Option<String>,
    pub charge: Option<i8>,
    pub radical: Option<AtomRadical>,
    pub isotope_mass: Option<u32>,
    pub hydrogen_count: Option<u8>,
    pub ring_bond_count: Option<RingBondCount>,
    pub substitution_count: Option<SubstitutionCount>,
    pub unsaturated: Option<UnsaturatedAtom>,
    pub link_atom: Option<LinkAtom>,
    pub attachment_point: Option<AttachmentPointType>,
    pub attachment_order: Option<Vec<(usize, u8)>>,
    pub rgroup_label: Option<u32>,
    pub atom_list_elements: Option<Vec<Element>>,
    pub atom_list_exclusion: Option<bool>,
}

// Accumulator for properties of a single bond
#[derive(Debug, Default)]
pub struct BondProperties {
    pub order_override: Option<u8>,
}

// Accumulator for properties of a single R-Group
#[derive(Debug, Default)]
pub struct RGroupProperties {
    pub dependent_label: Option<u32>,
    pub rgroup_or_h: Option<bool>,
    pub occurrence: Option<Vec<RGroupOccurrence>>,
}

// Accumulator for properties of a single S-Group
#[derive(Debug, Default)]
pub struct SGroupProperties {
    pub group_type: Option<SGroupType>,
    pub group_subtype: Option<SGroupSubtype>,
    pub label: Option<u32>,
    pub connectivity: Option<SGroupConnectivity>,
    pub expansion: Option<bool>,
    pub atom_indices: Option<Vec<usize>>,
    pub bond_indices: Option<Vec<usize>>,
    pub parent_atom_indices: Option<Vec<usize>>,
    pub multiplier: Option<SGroupMultiplier>,
    pub subscript: Option<String>,
    pub correspondence: Option<Vec<usize>>,
    pub bracket_coords: Option<SGroupBracketCoords>,
    pub connecting_bond: Option<SGroupConnectingBond>,
    pub hierarchy_parent: Option<usize>,
    pub component_number: Option<u32>,
    pub data: BTreeMap<String, SGroupData>,
    pub display: Option<SGroupDataDisplay>,
}

impl SGroupProperties {
    pub fn new(sgroup_type: SGroupType) -> Self {
        Self {
            group_type: Some(sgroup_type),
            ..Default::default()
        }
    }
}

/// Validate compatibility of SGroup type with subscript data type
fn sgroup_accepts_subscript(sgroup_type: SGroupType) -> bool {
    matches!(sgroup_type, SGroupType::Superatom)
}

/// Validate compatibility of SGroup type with multiplier data type
fn sgroup_accepts_multiplier(sgroup_type: SGroupType) -> bool {
    matches!(
        sgroup_type,
        SGroupType::MultipleGroup | SGroupType::RepeatingUnit
    )
}

/// Accumulator for molecular properties
#[derive(Debug)]
pub struct MoleculeProperties {
    context: Context,
    pub atom_properties: BTreeMap<usize, AtomProperties>,
    pub bond_properties: BTreeMap<usize, BondProperties>,
    pub rgroup_properties: BTreeMap<usize, RGroupProperties>,
    pub sgroup_properties: BTreeMap<usize, SGroupProperties>,
}

impl MoleculeProperties {
    pub fn new() -> Self {
        Self {
            context: Context::new(),
            atom_properties: BTreeMap::new(),
            bond_properties: BTreeMap::new(),
            rgroup_properties: BTreeMap::new(),
            sgroup_properties: BTreeMap::new(),
        }
    }

    pub fn add_entry(&mut self, entry: PropertyEntries, flags: ParseFlags) -> Result<()> {
        let extended_range = flags.contains(ParseFlags::EXTENDED_RANGE);
        match entry {
            PropertyEntries::AtomAliasEntry(e) => {
                let props = self.atom_properties.entry(e.atom_index).or_default();
                props.alias = Some(e.alias);
            }
            PropertyEntries::AtomValueEntry(e) => {
                let props = self.atom_properties.entry(e.atom_index).or_default();
                props.value = Some(e.value);
            }
            PropertyEntries::ChargeEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    props.charge = Some(entry.charge);
                }
            }
            PropertyEntries::RadicalEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    props.radical = convert_radical_type_code(entry.radical_type)?;
                }
            }
            PropertyEntries::IsotopeEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    props.isotope_mass = Some(entry.mass);
                }
            }
            PropertyEntries::RingBondCountEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    props.ring_bond_count = convert_ring_bond_count_code(entry.ring_bond_count)?;
                }
            }
            PropertyEntries::SubstitutionCountEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    props.substitution_count =
                        convert_substitution_count_code(entry.substitution_count, extended_range)?;
                }
            }
            PropertyEntries::UnsaturatedAtomEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    props.unsaturated = convert_unsaturated_atom_code(entry.unsaturated)?;
                }
            }
            PropertyEntries::LinkAtomEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    if props.link_atom.is_some() {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Duplicate link atom property for atom {}",
                            entry.atom_index
                        ))
                        .into());
                    }
                    props.link_atom = Some(LinkAtom {
                        repeat_count: entry.repeat_count,
                        subs_index1: entry.subs_index1,
                        subs_index2: entry.subs_index2,
                    });
                }
            }
            PropertyEntries::AtomListEntry(entry) => {
                let props = self.atom_properties.entry(entry.atom_index).or_default();
                props.atom_list_elements = Some(entry.elements);
                props.atom_list_exclusion = Some(entry.exclusion);
            }
            PropertyEntries::AttachmentPointEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    if props.attachment_point.is_some() {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Duplicate attachment point property for atom {}",
                            entry.atom_index
                        ))
                        .into());
                    }
                    props.attachment_point = convert_attachment_point_code(entry.attachment_type)?;
                }
            }
            PropertyEntries::AtomAttachmentOrderEntry(entry) => {
                let props = self.atom_properties.entry(entry.atom_index).or_default();
                props.attachment_order = Some(entry.attachments);
            }
            PropertyEntries::RGroupLabelEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    props.rgroup_label = Some(entry.label);
                }
            }
            PropertyEntries::RGroupLogicEntry(entry) => {
                let props = self
                    .rgroup_properties
                    .entry(entry.label as usize)
                    .or_default();
                props.dependent_label = entry.dependent_label;
                props.rgroup_or_h = Some(entry.rgroup_or_h);
                props.occurrence = Some(entry.occurrence);
            }
            PropertyEntries::SGroupTypeEntries(entries) => {
                for entry in entries {
                    if self.sgroup_properties.contains_key(&entry.sgroup_index) {
                        return Err(ValidationError::InvalidComponent(format!(
                            "Duplicate S-group type for index {}",
                            entry.sgroup_index
                        ))
                        .into());
                    }
                    let props = SGroupProperties::new(entry.sgroup_type);
                    self.sgroup_properties.insert(entry.sgroup_index, props);
                    self.context
                        .sgroup_types
                        .insert(entry.sgroup_index, entry.sgroup_type);
                }
            }
            PropertyEntries::SGroupSubtypeEntries(entries) => {
                for entry in entries {
                    let props = self
                        .sgroup_properties
                        .get_mut(&entry.sgroup_index)
                        .ok_or_else(|| {
                            DataError::InvalidFragment(format!(
                                "S-group subtype for undefined S-group {}",
                                entry.sgroup_index
                            ))
                        })?;
                    props.group_subtype = Some(entry.sgroup_subtype);
                }
            }
            PropertyEntries::SGroupLabelEntries(entries) => {
                for entry in entries {
                    let props = self
                        .sgroup_properties
                        .get_mut(&entry.sgroup_index)
                        .ok_or_else(|| {
                            DataError::InvalidFragment(format!(
                                "S-group label for undefined S-group {}",
                                entry.sgroup_index
                            ))
                        })?;
                    props.label = Some(entry.label);
                }
            }
            PropertyEntries::SGroupConnectivityEntries(entries) => {
                for entry in entries {
                    let props = self
                        .sgroup_properties
                        .get_mut(&entry.sgroup_index)
                        .ok_or_else(|| {
                            DataError::InvalidFragment(format!(
                                "S-group connectivity for undefined S-group {}",
                                entry.sgroup_index
                            ))
                        })?;
                    props.connectivity = Some(entry.connectivity);
                }
            }
            PropertyEntries::SGroupExpansionEntries(entries) => {
                for entry in entries {
                    let props = self
                        .sgroup_properties
                        .get_mut(&entry.sgroup_index)
                        .ok_or_else(|| {
                            DataError::InvalidFragment(format!(
                                "S-group expansion for undefined S-group {}",
                                entry.sgroup_index
                            ))
                        })?;
                    props.expansion = Some(true);
                }
            }
            PropertyEntries::SGroupAtomListEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        DataError::InvalidFragment(format!(
                            "S-group atom list for undefined S-group {}",
                            entry.sgroup_index
                        ))
                    })?;
                props.atom_indices = Some(entry.atom_indices);
            }
            PropertyEntries::SGroupBondListEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        DataError::InvalidFragment(format!(
                            "S-group bond list for undefined S-group {}",
                            entry.sgroup_index
                        ))
                    })?;
                props.bond_indices = Some(entry.bond_indices);
            }
            PropertyEntries::SGroupParentAtomEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        DataError::InvalidFragment(format!(
                            "S-group parent atom list for undefined S-group {}",
                            entry.sgroup_index
                        ))
                    })?;
                props.parent_atom_indices = Some(entry.atom_indices);
            }
            PropertyEntries::SGroupSubscriptEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        DataError::InvalidFragment(format!(
                            "S-group subscript for undefined S-group {}",
                            entry.sgroup_index
                        ))
                    })?;
                let sgroup_type = self
                    .context
                    .sgroup_types
                    .get(&entry.sgroup_index)
                    .ok_or_else(|| {
                        DataError::InvalidFragment(format!(
                            "S-group subscript for undefined S-group {}",
                            entry.sgroup_index
                        ))
                    })?;

                if !sgroup_accepts_subscript(*sgroup_type)
                    && !sgroup_accepts_multiplier(*sgroup_type)
                {
                    return Err(DataError::InvalidFragment(format!(
                        "S-group subscript and multiplier not allowed for S-group type {:?}",
                        sgroup_type,
                    ))
                    .into());
                } else if sgroup_accepts_subscript(*sgroup_type) {
                    if entry.subscript.is_none() {
                        return Err(DataError::InvalidFragment(format!(
                            "No subscript found for S-group type {:?}",
                            sgroup_type,
                        ))
                        .into());
                    }
                    props.subscript = entry.subscript;
                } else if sgroup_accepts_multiplier(*sgroup_type) {
                    if entry.multiplier.is_none() {
                        return Err(DataError::InvalidFragment(format!(
                            "No multiplier found for S-group type {:?}",
                            sgroup_type,
                        ))
                        .into());
                    }
                    props.multiplier = entry.multiplier;
                } else {
                    return Err(DataError::InvalidFragment(format!(
                        "S-group type {:?} cannot have a subscript and a multiplier",
                        sgroup_type,
                    ))
                    .into());
                }
            }
            PropertyEntries::SGroupCorrespondenceEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        DataError::InvalidFragment(format!(
                            "S-group correspondence for undefined S-group {}",
                            entry.sgroup_index
                        ))
                    })?;
                props.correspondence = Some(entry.bond_indices);
            }
            PropertyEntries::SGroupDisplayInfoEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        DataError::InvalidFragment(format!(
                            "S-group display info for undefined S-group {}",
                            entry.sgroup_index
                        ))
                    })?;
                props.bracket_coords = Some(SGroupBracketCoords {
                    bracket1: (entry.bracket_coords[0], entry.bracket_coords[1]),
                    bracket2: (entry.bracket_coords[2], entry.bracket_coords[3]),
                });
            }
            PropertyEntries::SGroupConnectingBondEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        DataError::InvalidFragment(format!(
                            "S-group connecting bond for undefined S-group {}",
                            entry.sgroup_index
                        ))
                    })?;
                props.connecting_bond = Some(SGroupConnectingBond {
                    bond_index: entry.bond_index,
                    bond_vector: entry.bond_vector,
                });
            }
            PropertyEntries::SGroupDataDescriptionEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        DataError::InvalidFragment(format!(
                            "S-group data description for undefined S-group {}",
                            entry.sgroup_index
                        ))
                    })?;
                props.data.insert(
                    entry.field_name.clone(),
                    SGroupData {
                        field_type: entry.field_type,
                        field_units: entry.field_units,
                        query_identifier: entry.query_identifier,
                        data_query_operator: entry.data_query_operator,
                        data_content: None,
                    },
                );

                self.context.current_sgroup_index = Some(entry.sgroup_index);
                self.context.current_data_field = Some(entry.field_name);
            }
            PropertyEntries::SGroupDataDisplayEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        DataError::InvalidFragment(format!(
                            "S-group data description for undefined S-group {}",
                            entry.sgroup_index
                        ))
                    })?;
                props.display = Some(SGroupDataDisplay {
                    coords: entry.coords,
                    display_type: entry.display_type,
                    display_placement: entry.display_placement,
                    display_units: entry.display_units,
                    display_chars: entry.display_chars,
                });
            }
            PropertyEntries::SGroupDataEntry(entry) => match entry {
                SGroupDataEntry::Continuation {
                    sgroup_index,
                    data_content,
                } => {
                    if self.context.current_sgroup_index.is_none() {
                        return Err(DataError::InvalidFragment(
                            "SCD entry found without SDT context".to_string(),
                        )
                        .into());
                    }

                    let context_sgroup = self.context.current_sgroup_index.unwrap();
                    if sgroup_index != context_sgroup {
                        return Err(DataError::InvalidFragment(format!(
                            "SCD sgroup_index {} doesn't match context sgroup {}",
                            sgroup_index, context_sgroup
                        ))
                        .into());
                    }

                    if self.context.current_data_content.is_none() {
                        self.context.current_data_content = Some(Vec::new());
                    }

                    let buffer = self.context.current_data_content.as_mut().unwrap();
                    if buffer.is_empty() {
                        buffer.push(data_content.clone());
                    } else {
                        let last_line = buffer.last_mut().unwrap();
                        last_line.push_str(&data_content);
                    }
                }

                SGroupDataEntry::EndWithData {
                    sgroup_index,
                    data_content,
                } => {
                    if let Some(context_sgroup) = self.context.current_sgroup_index {
                        if sgroup_index != context_sgroup {
                            return Err(DataError::InvalidFragment(format!(
                                "SED sgroup_index {} doesn't match context sgroup {}",
                                sgroup_index, context_sgroup
                            ))
                            .into());
                        }

                        if self.context.current_data_content.is_none() {
                            self.context.current_data_content = Some(Vec::new());
                        }

                        let buffer = self.context.current_data_content.as_mut().unwrap();
                        if buffer.is_empty() {
                            buffer.push(data_content.clone());
                        } else {
                            let last_line = buffer.last_mut().unwrap();
                            last_line.push_str(&data_content);
                        }

                        self.finalize_sgroup_data(sgroup_index, None)?;
                    } else {
                        self.finalize_sgroup_data(sgroup_index, Some(data_content.clone()))?;
                    }
                }

                SGroupDataEntry::EndBlank { sgroup_index } => {
                    if self.context.current_sgroup_index.is_none() {
                        return Err(DataError::InvalidFragment(
                            "Blank SED entry found without SDT context".to_string(),
                        )
                        .into());
                    }

                    let context_sgroup = self.context.current_sgroup_index.unwrap();
                    if sgroup_index != context_sgroup {
                        return Err(DataError::InvalidFragment(format!(
                            "Blank SED sgroup_index {} doesn't match context sgroup {}",
                            sgroup_index, context_sgroup
                        ))
                        .into());
                    }

                    self.finalize_sgroup_data(sgroup_index, None)?;
                }
            },
            PropertyEntries::SGroupHierarchyEntries(entries) => {
                for entry in entries {
                    let props = self
                        .sgroup_properties
                        .get_mut(&entry.sgroup_index)
                        .ok_or_else(|| {
                            DataError::InvalidFragment(format!(
                                "S-group hierarchy for undefined S-group {}",
                                entry.sgroup_index
                            ))
                        })?;
                    props.hierarchy_parent = Some(entry.parent_sgroup_index);
                }
            }
            PropertyEntries::SGroupComponentEntries(entries) => {
                for entry in entries {
                    let props = self
                        .sgroup_properties
                        .get_mut(&entry.sgroup_index)
                        .ok_or_else(|| {
                            DataError::InvalidFragment(format!(
                                "S-group component for undefined S-group {}",
                                entry.sgroup_index
                            ))
                        })?;
                    props.component_number = Some(entry.component_number);
                }
            }
            PropertyEntries::ZeroBondOrderEntries(entries) => {
                for entry in entries {
                    self.bond_properties.insert(
                        entry.bond_index,
                        BondProperties {
                            order_override: Some(entry.bond_order),
                        },
                    );
                }
            }
            PropertyEntries::ZeroAtomChargeEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    props.charge = Some(entry.charge);
                }
            }
            PropertyEntries::AtomHydrogenCountEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    props.hydrogen_count = Some(entry.hydrogen_count);
                }
            }
            PropertyEntries::End => {}
        }
        Ok(())
    }

    /// Apply all properties to Molecule
    pub fn update_molecule(&mut self, molecule: &mut Molecule, flags: ParseFlags) -> Result<()> {
        let extended_isotopes = flags.contains(ParseFlags::EXTENDED_ISOTOPES);
        for (atom_index, props) in &self.atom_properties {
            if props.alias.is_some() {
                self.apply_atom_alias(*atom_index, props, molecule)?;
            }
            if props.value.is_some() {
                self.apply_atom_value(*atom_index, props, molecule)?;
            }
            if props.charge.is_some() {
                self.apply_atom_charge(*atom_index, props, molecule)?;
            }
            if props.radical.is_some() {
                self.apply_atom_radical(*atom_index, props, molecule)?;
            }
            if props.isotope_mass.is_some() {
                self.apply_atom_isotope(*atom_index, props, molecule, extended_isotopes)?;
            }
            if props.hydrogen_count.is_some() {
                self.apply_atom_hydrogen_count(*atom_index, props, molecule)?;
            }
        }
        for (bond_index, props) in &self.bond_properties {
            if props.order_override.is_some() {
                self.apply_atom_zero_order_bond(*bond_index, props, molecule)?;
            }
        }
        Ok(())
    }

    /// Apply all properties to MoleculeLike
    pub fn update_moleculelike(
        &mut self,
        molecule: &mut MoleculeLike,
        flags: ParseFlags,
    ) -> Result<()> {
        let extended_isotopes = flags.contains(ParseFlags::EXTENDED_ISOTOPES);
        for (&atom_index, props) in &self.atom_properties {
            let atom = molecule
                .atom_mut(atom_index)
                .ok_or(DataError::MissingAtomIndex(atom_index))?;

            if props.alias.is_some() {
                self.apply_atomlike_alias(props, atom)?;
            }
            if props.value.is_some() {
                self.apply_atomlike_value(props, atom)?;
            }
            if props.charge.is_some() {
                self.apply_atomlike_charge(props, atom)?;
            }
            if props.radical.is_some() {
                self.apply_atomlike_radical(props, atom)?;
            }
            if props.isotope_mass.is_some() {
                self.apply_atomlike_isotope(props, atom, extended_isotopes)?;
            }
            if props.hydrogen_count.is_some() {
                self.apply_atomlike_hydrogen_count(props, atom)?;
            }
            if props.ring_bond_count.is_some() {
                self.apply_ring_bond_count(props, atom)?;
            }
            if props.substitution_count.is_some() {
                self.apply_substitution_count(props, atom)?;
            }
            if props.unsaturated.is_some() {
                self.apply_unsaturated_atom(props, atom)?;
            }
            if props.link_atom.is_some() {
                self.apply_link_atom(props, atom)?;
            }
            if props.atom_list_elements.is_some() {
                self.apply_atom_list(props, atom)?;
            }
            if props.attachment_point.is_some() {
                self.apply_attachment_point(props, atom)?;
            }
            if props.attachment_order.is_some() {
                self.apply_attachment_order(props, atom)?;
            }
            if props.rgroup_label.is_some() {
                self.apply_rgroup_label(props, atom)?;
            }
        }
        for (rgroup_label, props) in &self.rgroup_properties {
            self.apply_rgroup_logic(*rgroup_label, props, molecule)?;
        }
        self.validate_sgroup_data()?;
        self.apply_sgroup(molecule)?;
        self.apply_atomlike_zero_order_bonds(molecule)?;
        Ok(())
    }

    fn apply_atomlike_alias(&self, props: &AtomProperties, atom: &mut AtomLike) -> Result<()> {
        let alias = props.alias.as_ref().unwrap();
        if let Some(existing_alias) = atom.properties.get("molFileAlias") {
            if existing_alias != alias {
                return Err(ValidationError::InvalidComponent(format!(
                    "Atom alias conflict: existing '{}' vs new '{}'",
                    existing_alias, alias
                ))
                .into());
            }
        }
        atom.properties
            .insert("molFileAlias".to_string(), alias.clone());
        Ok(())
    }

    fn apply_atom_alias(
        &self,
        atom_index: usize,
        props: &AtomProperties,
        molecule: &mut Molecule,
    ) -> Result<()> {
        let alias = props.alias.as_ref().unwrap();
        let atom = molecule
            .atom_mut(atom_index)
            .ok_or(DataError::MissingAtomIndex(atom_index))?;
        if let Some(existing_alias) = atom.properties.get("molFileAlias") {
            if existing_alias != alias {
                return Err(ValidationError::InvalidComponent(format!(
                    "Atom alias conflict for atom {}: existing '{}' vs new '{}'",
                    atom_index, existing_alias, alias
                ))
                .into());
            }
        }
        atom.properties
            .insert("molFileAlias".to_string(), alias.clone());
        Ok(())
    }

    fn apply_atomlike_value(&self, props: &AtomProperties, atom: &mut AtomLike) -> Result<()> {
        let value = props.value.as_ref().unwrap();
        if let Some(existing_value) = atom.properties.get("molFileValue") {
            if existing_value != value {
                return Err(ValidationError::InvalidComponent(format!(
                    "Atom value conflict: existing '{}' vs new '{}'",
                    existing_value, value
                ))
                .into());
            }
        }
        atom.properties
            .insert("molFileValue".to_string(), value.clone());
        Ok(())
    }

    fn apply_atom_value(
        &self,
        atom_index: usize,
        props: &AtomProperties,
        molecule: &mut Molecule,
    ) -> Result<()> {
        let value = props.value.as_ref().unwrap();
        let atom = molecule
            .atom_mut(atom_index)
            .ok_or(DataError::MissingAtomIndex(atom_index))?;
        if let Some(existing_value) = atom.properties.get("molFileValue") {
            if existing_value != value {
                return Err(ValidationError::InvalidComponent(format!(
                    "Atom value conflict for atom {}: existing '{}' vs new '{}'",
                    atom_index, existing_value, value
                ))
                .into());
            }
        }
        atom.properties
            .insert("molFileValue".to_string(), value.clone());
        Ok(())
    }

    fn apply_atomlike_charge(&self, props: &AtomProperties, atom: &mut AtomLike) -> Result<()> {
        let charge = props.charge.unwrap();
        atom.charge = charge;
        atom.radical = None;
        Ok(())
    }

    fn apply_atom_charge(
        &self,
        atom_index: usize,
        props: &AtomProperties,
        molecule: &mut Molecule,
    ) -> Result<()> {
        let charge = props.charge.unwrap();
        let atom = molecule
            .atom_mut(atom_index)
            .ok_or(DataError::MissingAtomIndex(atom_index))?;
        atom.charge = charge;
        atom.radical = None;
        Ok(())
    }

    fn apply_atomlike_radical(&self, props: &AtomProperties, atom: &mut AtomLike) -> Result<()> {
        atom.radical = props.radical;
        atom.charge = 0;
        Ok(())
    }

    fn apply_atom_radical(
        &self,
        atom_index: usize,
        props: &AtomProperties,
        molecule: &mut Molecule,
    ) -> Result<()> {
        let atom = molecule
            .atom_mut(atom_index)
            .ok_or(DataError::MissingAtomIndex(atom_index))?;
        atom.radical = props.radical;
        atom.charge = 0;
        Ok(())
    }

    fn apply_atomlike_isotope(
        &self,
        props: &AtomProperties,
        atom: &mut AtomLike,
        extended_isotopes: bool,
    ) -> Result<()> {
        let element = match &atom.symbol {
            AtomSymbol::Element(element) => Ok::<_, umol::Error>(*element),
            AtomSymbol::NamedIsotope(isotope) => Ok::<_, umol::Error>(isotope.element()),
            _ => Err(ValidationError::InvalidComponent(format!(
                "Cannot set isotope for atom: {:?}",
                atom.symbol
            ))
            .into()),
        }?;
        let mass = convert_atom_isotope_mass_number(
            element,
            props.isotope_mass.unwrap(),
            extended_isotopes,
        )?;
        if let Some(existing) = atom.isotope_mass {
            if let Some(mass) = mass {
                if existing != mass {
                    return Err(ValidationError::InvalidComponent(format!(
                        "Isotope conflict: existing '{}' vs new '{}'",
                        existing, mass
                    ))
                    .into());
                }
            }
        }
        atom.isotope_mass = mass;
        Ok(())
    }

    fn apply_atom_isotope(
        &self,
        atom_index: usize,
        props: &AtomProperties,
        molecule: &mut Molecule,
        extended_isotopes: bool,
    ) -> Result<()> {
        let atom = molecule
            .atom_mut(atom_index)
            .ok_or(DataError::MissingAtomIndex(atom_index))?;
        let mass = convert_atom_isotope_mass_number(
            atom.element,
            props.isotope_mass.unwrap(),
            extended_isotopes,
        )?;
        if let Some(existing) = atom.isotope_mass {
            if let Some(mass) = mass {
                if existing != mass {
                    return Err(ValidationError::InvalidComponent(format!(
                        "Isotope conflict for atom {}: existing '{}' vs new '{}'",
                        atom_index, existing, mass
                    ))
                    .into());
                }
            }
        }
        atom.isotope_mass = mass;
        Ok(())
    }

    fn apply_atomlike_hydrogen_count(
        &self,
        props: &AtomProperties,
        atom: &mut AtomLike,
    ) -> Result<()> {
        let hydrogen_count = props.hydrogen_count.unwrap();
        atom.hydrogen_count = Some(hydrogen_count);
        Ok(())
    }

    fn apply_atom_hydrogen_count(
        &self,
        atom_index: usize,
        props: &AtomProperties,
        molecule: &mut Molecule,
    ) -> Result<()> {
        let hydrogen_count = props.hydrogen_count.unwrap();
        let atom = molecule
            .atom_mut(atom_index)
            .ok_or(DataError::MissingAtomIndex(atom_index))?;
        atom.hydrogen_count = Some(hydrogen_count);
        Ok(())
    }

    fn apply_ring_bond_count(&self, props: &AtomProperties, atom: &mut AtomLike) -> Result<()> {
        if let Some(existing) = atom.ring_bond_count {
            if Some(existing) != props.ring_bond_count {
                return Err(ValidationError::InvalidComponent(format!(
                    "Ring bond count conflict: existing {:?} vs new {:?}",
                    existing, props.ring_bond_count
                ))
                .into());
            }
        }
        atom.ring_bond_count = props.ring_bond_count;
        Ok(())
    }

    fn apply_substitution_count(&self, props: &AtomProperties, atom: &mut AtomLike) -> Result<()> {
        if let Some(existing) = atom.substitution_count {
            if Some(existing) != props.substitution_count {
                return Err(ValidationError::InvalidComponent(format!(
                    "Substitution count conflict: existing {:?} vs new {:?}",
                    existing, props.substitution_count
                ))
                .into());
            }
        }
        atom.substitution_count = props.substitution_count;
        Ok(())
    }

    fn apply_unsaturated_atom(&self, props: &AtomProperties, atom: &mut AtomLike) -> Result<()> {
        atom.unsaturated = props.unsaturated;
        Ok(())
    }

    fn apply_link_atom(&self, props: &AtomProperties, atom: &mut AtomLike) -> Result<()> {
        if atom.link_atom.is_some() {
            return Err(ValidationError::InvalidComponent("Link atom conflict".to_string()).into());
        }
        atom.link_atom = props.link_atom;
        Ok(())
    }

    fn apply_atom_list(&self, props: &AtomProperties, atom: &mut AtomLike) -> Result<()> {
        if !matches!(
            atom.symbol,
            AtomSymbol::Element(_) | AtomSymbol::AtomList(_)
        ) {
            return Err(ValidationError::InvalidComponent(format!(
                "Atom list can only be applied to an element or atom list, not {:?}",
                atom.symbol
            ))
            .into());
        }
        atom.symbol = AtomSymbol::AtomList(AtomList {
            elements: props.atom_list_elements.clone().unwrap(),
            exclusion: props.atom_list_exclusion.unwrap(),
        });
        Ok(())
    }

    fn apply_attachment_point(&self, props: &AtomProperties, atom: &mut AtomLike) -> Result<()> {
        if atom.attachment_point.is_some() {
            return Err(
                ValidationError::InvalidComponent("Attachment point conflict".to_string()).into(),
            );
        }
        atom.attachment_point = props.attachment_point;
        Ok(())
    }

    fn apply_attachment_order(&self, props: &AtomProperties, atom: &mut AtomLike) -> Result<()> {
        if atom.attachment_order.is_some() {
            return Err(
                ValidationError::InvalidComponent("Attachment order conflict".to_string()).into(),
            );
        }
        atom.attachment_order = props.attachment_order.clone();
        Ok(())
    }

    fn apply_rgroup_label(&self, props: &AtomProperties, atom: &mut AtomLike) -> Result<()> {
        if let Some(new_label) = props.rgroup_label {
            match &mut atom.symbol {
                AtomSymbol::Element(_) => {
                    // Convert element to R-group with the label
                    atom.symbol = AtomSymbol::RGroup(RGroup::new(Some(new_label)));
                }
                AtomSymbol::RGroup(rgroup) => {
                    // Verify existing label matches or is None
                    if let Some(existing_label) = rgroup.label {
                        if existing_label != new_label {
                            return Err(ValidationError::InvalidComponent(format!(
                                "R-group label conflict: existing '{}' vs new '{}'",
                                existing_label, new_label
                            ))
                            .into());
                        }
                    } else {
                        rgroup.label = Some(new_label);
                    }
                }
                _ => {
                    return Err(ValidationError::InvalidComponent(format!(
                        "R-group label can only be applied to an element or R-group, not {:?}",
                        atom.symbol
                    ))
                    .into());
                }
            }
        }
        Ok(())
    }

    fn apply_rgroup_logic(
        &self,
        rgroup_label: usize,
        props: &RGroupProperties,
        molecule: &mut MoleculeLike,
    ) -> Result<()> {
        for i in 0..molecule.atom_count() {
            if let Some(atom) = molecule.atom_mut(i) {
                if let AtomSymbol::RGroup(rgroup) = &mut atom.symbol {
                    if rgroup.label == Some(rgroup_label as u32) {
                        rgroup.dependent_label = props.dependent_label;
                        if let Some(rgroup_or_h) = props.rgroup_or_h {
                            rgroup.rgroup_or_h = rgroup_or_h;
                        }
                        if let Some(occurrence) = &props.occurrence {
                            rgroup.occurrence = occurrence.clone();
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn apply_atomlike_zero_order_bonds(&self, molecule: &mut MoleculeLike) -> Result<()> {
        for (bond_index, props) in &self.bond_properties {
            if let Some(bond) = molecule.bond_mut(*bond_index) {
                bond.bond_type =
                    convert_bondlike_type_code(props.order_override.unwrap(), true, true)?;
            } else {
                return Err(DataError::InvalidFragment(format!(
                    "Zero-order bond for undefined bond {}",
                    bond_index
                ))
                .into());
            }
        }
        Ok(())
    }

    fn apply_atom_zero_order_bond(
        &self,
        bond_index: usize,
        props: &BondProperties,
        molecule: &mut Molecule,
    ) -> Result<()> {
        let bond = molecule
            .bond_mut(bond_index)
            .ok_or(DataError::MissingBondIndex(bond_index))?;
        bond.bond_type = convert_bondlike_type_code(props.order_override.unwrap(), true, true)?;
        Ok(())
    }

    fn apply_sgroup(&self, molecule: &mut MoleculeLike) -> Result<()> {
        for (sgroup_index, props) in &self.sgroup_properties {
            let sgroup_type = props.group_type.ok_or_else(|| {
                DataError::InvalidFragment(format!("S-group {} has no type", sgroup_index))
            })?;
            let mut sgroup = SGroup::new(sgroup_type);
            sgroup.label = props.label.or(Some(*sgroup_index as u32));
            sgroup.group_subtype = props.group_subtype;
            sgroup.connectivity = props.connectivity;
            if let Some(expansion) = props.expansion {
                sgroup.expansion = expansion;
            }
            if let Some(atom_indices) = &props.atom_indices {
                sgroup.atom_indices = atom_indices.clone();
            }
            if let Some(bond_indices) = &props.bond_indices {
                sgroup.bond_indices = bond_indices.clone();
            }
            if let Some(parent_atom_indices) = &props.parent_atom_indices {
                sgroup.parent_atom_indices = Some(parent_atom_indices.clone());
            }
            if let Some(subscript) = &props.subscript {
                sgroup.subscript = Some(subscript.clone());
            }
            if let Some(multiplier) = &props.multiplier {
                sgroup.multiplier = Some(*multiplier);
            }
            if let Some(correspondence) = &props.correspondence {
                sgroup.correspondence = Some(correspondence.clone());
            }
            if let Some(bracket_coords) = &props.bracket_coords {
                sgroup.bracket_coords = Some(*bracket_coords);
            }
            if let Some(connecting_bond) = &props.connecting_bond {
                sgroup.connecting_bond = Some(*connecting_bond);
            }
            if let Some(hierarchy_parent) = &props.hierarchy_parent {
                sgroup.hierarchy_parent = Some(*hierarchy_parent);
            }
            if let Some(component_number) = &props.component_number {
                sgroup.component_number = Some(*component_number);
            }
            if !props.data.is_empty() {
                sgroup.data = props.data.clone();
            }
            if let Some(display) = &props.display {
                sgroup.display = Some(*display);
            }
            molecule.add_sgroup(*sgroup_index, sgroup);
        }
        Ok(())
    }

    fn validate_sgroup_data(&mut self) -> Result<()> {
        if self.context.current_sgroup_index.is_some() {
            let sgroup_index = self.context.current_sgroup_index.unwrap();
            self.finalize_sgroup_data(sgroup_index, None)?;
            self.context.current_sgroup_index = None;
            self.context.current_data_field = None;
        }
        Ok(())
    }

    fn finalize_sgroup_data(
        &mut self,
        sgroup_index: usize,
        sed_content: Option<String>,
    ) -> Result<()> {
        let data_field = self.context.current_data_field.as_ref().ok_or_else(|| {
            DataError::InvalidFragment(
                "No data field context when finalizing S-group data".to_string(),
            )
        })?;

        let props = self
            .sgroup_properties
            .get_mut(&sgroup_index)
            .ok_or_else(|| {
                DataError::InvalidFragment(format!(
                    "S-group data content for undefined S-group {}",
                    sgroup_index
                ))
            })?;

        let has_sed_content = sed_content.is_some();

        let buffer = self.context.current_data_content.as_mut();
        let mut content = if let Some(buffer) = buffer {
            if buffer.is_empty() {
                String::new()
            } else {
                buffer.join("")
            }
        } else {
            String::new()
        };

        if let Some(sed_data) = sed_content {
            content.push_str(&sed_data);
        }

        if content.len() > 200 {
            content.truncate(200);
        }
        let content = content.trim_end().to_string();

        if let Some(data) = props.data.get_mut(data_field) {
            if !content.is_empty() || has_sed_content {
                if let Some(existing_content) = &mut data.data_content {
                    existing_content.push(content);
                } else {
                    data.data_content = Some(vec![content]);
                }
            }
        }

        self.context.current_data_content = None;

        // TODO: Check attachment points are valid
        // From RDKit C++ code
        //     const RWMol *mol, std::pair<const int, SubstanceGroup> &sgroup) {
        //     bool res = true;
        //     int nAtoms = static_cast<int>(mol->getNumAtoms());
        //     std::vector<SubstanceGroup::AttachPoint> &attachPoints umol-models-graph/benches/parsing_bench.rs=
        //         sgroup.second.getAttachPoints();
        //     for (auto &attachPoint : attachPoints) {
        //         if (attachPoint.lvIdx == nAtoms) {
        //         const std::vector<unsigned int> &bonds = sgroup.second.getBonds();
        //         if (bonds.size() == 1) {
        //             const auto bond = mol->getBondWithIdx(bonds.front());
        //             if (bond->getBeginAtomIdx() == attachPoint.aIdx ||
        //                 bond->getEndAtomIdx() == attachPoint.aIdx) {
        //             attachPoint.lvIdx = bond->getOtherAtomIdx(attachPoint.aIdx);
        //             }
        //         }
        //         }
        //         if (attachPoint.lvIdx == nAtoms) {
        //         BOOST_LOG(rdWarningLog)
        //             << "Could not infer missing lvIdx on malformed SAP line for SGroup "
        //             << sgroup.first << std::endl;
        //         res = false;
        //         }
        //     }
        //     return res;

        Ok(())
    }
}

#[cfg(test)]
mod tests;
