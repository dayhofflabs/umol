//! Accumulator for molecular properties from a CTAB file

use std::collections::BTreeMap;

use crate::io::ctfile::error::SemanticError;
use umol_data::Element;

type Result<T> = std::result::Result<T, SemanticError>;

use super::context::Context;
use super::convert::{
    convert_atom_isotope_mass_number, convert_attachment_point_code, convert_extended_bond_type_code,
    convert_radical_type_code, convert_ring_bond_count_code, convert_substitution_count_code,
    convert_unsaturated_atom_code,
};
use crate::io::ctab::config::CtabParseFlags;
use crate::io::ctab::parser::properties::{PropertyEntries, SGroupDataEntry};
use crate::table_ir::{
    AtomList, AtomSymbol, AttachmentPointType, CtfileData, ExtendedMolecule,
    LinkAtom, RGroup, RGroupOccurrence, RingBondCount, SGroup, SGroupBracketCoords,
    SGroupConnectingBond, SGroupConnectivity, SGroupData, SGroupDataDisplay, SGroupMultiplier,
    SGroupSubtype, SGroupType, SubstitutionCount, UnsaturatedAtom,
};

// Accumulator for properties of a single atom
#[derive(Debug, Default)]
pub(super) struct AtomProperties {
    pub alias: Option<String>,
    pub value: Option<String>,
    pub charge: Option<i8>,
    pub unpaired_e: Option<u8>,
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
pub(super) struct BondProperties {
    pub order_override: Option<u8>,
}

// Accumulator for properties of a single R-Group
#[derive(Debug, Default)]
pub(super) struct RGroupProperties {
    pub dependent_label: Option<u32>,
    pub rgroup_or_h: Option<bool>,
    pub occurrence: Option<Vec<RGroupOccurrence>>,
}

// Accumulator for properties of a single S-Group
#[derive(Debug, Default)]
pub(super) struct SGroupProperties {
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
    pub(crate) fn new(sgroup_type: SGroupType) -> Self {
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
pub(super) struct MoleculeProperties {
    context: Context,
    pub atom_properties: BTreeMap<usize, AtomProperties>,
    pub bond_properties: BTreeMap<usize, BondProperties>,
    pub rgroup_properties: BTreeMap<usize, RGroupProperties>,
    pub sgroup_properties: BTreeMap<usize, SGroupProperties>,
}

impl MoleculeProperties {
    pub(crate) fn new() -> Self {
        Self {
            context: Context::new(),
            atom_properties: BTreeMap::new(),
            bond_properties: BTreeMap::new(),
            rgroup_properties: BTreeMap::new(),
            sgroup_properties: BTreeMap::new(),
        }
    }

    pub(crate) fn add_entry(&mut self, entry: PropertyEntries, flags: CtabParseFlags) -> Result<()> {
        let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
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
                    props.unpaired_e = convert_radical_type_code(entry.radical_type)?;
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
                        return Err(SemanticError::DuplicateProperty(
                            format!("link atom for atom {}", entry.atom_index)
                        ).into());
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
                        return Err(SemanticError::DuplicateProperty(
                            format!("attachment point for atom {}", entry.atom_index)
                        ).into());
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
                        return Err(SemanticError::DuplicateProperty(
                            format!("S-group type for index {}", entry.sgroup_index)
                        ).into());
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
                        .ok_or_else(|| SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "subtype",
                        })?;
                    props.group_subtype = Some(entry.sgroup_subtype);
                }
            }
            PropertyEntries::SGroupLabelEntries(entries) => {
                for entry in entries {
                    let props = self
                        .sgroup_properties
                        .get_mut(&entry.sgroup_index)
                        .ok_or_else(|| SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "label",
                        })?;
                    props.label = Some(entry.label);
                }
            }
            PropertyEntries::SGroupConnectivityEntries(entries) => {
                for entry in entries {
                    let props = self
                        .sgroup_properties
                        .get_mut(&entry.sgroup_index)
                        .ok_or_else(|| SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "connectivity",
                        })?;
                    props.connectivity = Some(entry.connectivity);
                }
            }
            PropertyEntries::SGroupExpansionEntries(entries) => {
                for entry in entries {
                    let props = self
                        .sgroup_properties
                        .get_mut(&entry.sgroup_index)
                        .ok_or_else(|| SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "expansion",
                        })?;
                    props.expansion = Some(true);
                }
            }
            PropertyEntries::SGroupAtomListEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "atom list",
                        }
                    })?;
                props.atom_indices = Some(entry.atom_indices);
            }
            PropertyEntries::SGroupBondListEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "bond list",
                        }
                    })?;
                props.bond_indices = Some(entry.bond_indices);
            }
            PropertyEntries::SGroupParentAtomEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "parent atom list",
                        }
                    })?;
                props.parent_atom_indices = Some(entry.atom_indices);
            }
            PropertyEntries::SGroupSubscriptEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "subscript",
                        }
                    })?;
                let sgroup_type = self
                    .context
                    .sgroup_types
                    .get(&entry.sgroup_index)
                    .ok_or_else(|| {
                        SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "subscript",
                        }
                    })?;

                if !sgroup_accepts_subscript(*sgroup_type)
                    && !sgroup_accepts_multiplier(*sgroup_type)
                {
                    return Err(SemanticError::SGroupTypeConstraint {
                        sgroup_type: *sgroup_type,
                        message: "subscript and multiplier not allowed",
                    }.into());
                } else if sgroup_accepts_subscript(*sgroup_type) {
                    if entry.subscript.is_none() {
                        return Err(SemanticError::SGroupTypeConstraint {
                            sgroup_type: *sgroup_type,
                            message: "subscript required",
                        }.into());
                    }
                    props.subscript = entry.subscript;
                } else if sgroup_accepts_multiplier(*sgroup_type) {
                    if entry.multiplier.is_none() {
                        return Err(SemanticError::SGroupTypeConstraint {
                            sgroup_type: *sgroup_type,
                            message: "multiplier required",
                        }.into());
                    }
                    props.multiplier = entry.multiplier;
                } else {
                    return Err(SemanticError::SGroupTypeConstraint {
                        sgroup_type: *sgroup_type,
                        message: "cannot have both subscript and multiplier",
                    }.into());
                }
            }
            PropertyEntries::SGroupCorrespondenceEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "correspondence",
                        }
                    })?;
                props.correspondence = Some(entry.bond_indices);
            }
            PropertyEntries::SGroupDisplayInfoEntry(entry) => {
                let props = self
                    .sgroup_properties
                    .get_mut(&entry.sgroup_index)
                    .ok_or_else(|| {
                        SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "display info",
                        }
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
                        SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "connecting bond",
                        }
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
                        SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "data description",
                        }
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
                        SemanticError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "data description",
                        }
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
                        return Err(SemanticError::MissingSGroupDataContext);
                    }

                    let context_sgroup = self.context.current_sgroup_index.unwrap();
                    if sgroup_index != context_sgroup {
                        return Err(SemanticError::SGroupIndexMismatch {
                            expected: context_sgroup,
                            actual: sgroup_index,
                        }.into());
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
                            return Err(SemanticError::SGroupIndexMismatch {
                                expected: context_sgroup,
                                actual: sgroup_index,
                            }.into());
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
                        return Err(SemanticError::MissingSGroupDataContext);
                    }

                    let context_sgroup = self.context.current_sgroup_index.unwrap();
                    if sgroup_index != context_sgroup {
                        return Err(SemanticError::SGroupIndexMismatch {
                            expected: context_sgroup,
                            actual: sgroup_index,
                        }.into());
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
                            SemanticError::UndefinedSGroup {
                                index: entry.sgroup_index,
                                property: "hierarchy",
                            }
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
                            SemanticError::UndefinedSGroup {
                                index: entry.sgroup_index,
                                property: "component",
                            }
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
            SemanticError::MissingSGroupDataContext
        })?;

        let props = self
            .sgroup_properties
            .get_mut(&sgroup_index)
            .ok_or_else(|| {
                SemanticError::UndefinedSGroup {
                    index: sgroup_index,
                    property: "data content",
                }
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
        Ok(())
    }

    /// Apply all properties to ExtendedMolecule (TableIR type)
    pub(crate) fn update_extended_molecule(
        &mut self,
        molecule: &mut ExtendedMolecule,
        flags: CtabParseFlags,
    ) -> Result<()> {
        let extended_isotopes = flags.contains(CtabParseFlags::EXTENDED_ISOTOPES);

        // Apply atom properties
        for (&atom_idx, props) in &self.atom_properties {
            let Some(atom) = molecule.atoms.get_mut(atom_idx) else {
                return Err(SemanticError::IndexOutOfBounds(atom_idx));
            };

            // Apply alias
            if let Some(ref alias) = props.alias {
                atom.properties
                    .insert("molFileAlias".to_string(), alias.clone());
            }

            // Apply value
            if let Some(ref value) = props.value {
                atom.properties
                    .insert("molFileValue".to_string(), value.clone());
            }

            // Apply charge
            if let Some(charge) = props.charge {
                atom.charge = Some(charge);
                atom.unpaired_e = None; // charge overrides radical from atom block
            }

            // Apply radical (unpaired electrons)
            if let Some(unpaired_e) = props.unpaired_e {
                atom.unpaired_e = Some(unpaired_e);
                atom.charge = None; // radical overrides charge from atom block
            }

            // Apply isotope
            if let Some(isotope) = props.isotope_mass {
                // Validate isotope for elements, apply directly for named isotopes
                match atom.symbol {
                    AtomSymbol::Element(element) => {
                        let validated =
                            convert_atom_isotope_mass_number(element, isotope, extended_isotopes)?;
                        atom.isotope_mass = validated;
                    }
                    AtomSymbol::NamedIsotope(named) => {
                        // For named isotopes (D, T), validate against the element
                        let validated = convert_atom_isotope_mass_number(
                            named.element(),
                            isotope,
                            extended_isotopes,
                        )?;
                        atom.isotope_mass = validated;
                    }
                    _ => {
                        // For other symbol types (queries, etc.), just set the isotope
                        atom.isotope = Some(isotope);
                    }
                }
            }

            // Apply hydrogen count
            if let Some(h_count) = props.hydrogen_count {
                atom.hydrogens = Some(h_count);
            }

            // Apply ring bond count
            if let Some(ref rbc) = props.ring_bond_count {
                atom.ring_bond_count = Some(rbc.clone());
            }

            // Apply substitution count
            if let Some(ref sub) = props.substitution_count {
                atom.substitution_count = Some(sub.clone());
            }

            // Apply unsaturated
            if props.unsaturated.is_some() {
                atom.unsaturated = Some(UnsaturatedAtom);
            }

            // Apply link atom
            if let Some(ref link) = props.link_atom {
                atom.link_atom = Some(link.clone());
            }

            // Apply atom list
            if let Some(ref elements) = props.atom_list_elements {
                let exclusion = props.atom_list_exclusion.unwrap_or(false);
                atom.symbol = AtomSymbol::AtomList(AtomList {
                    elements: elements.clone(),
                    exclusion,
                });
            }

            // Apply attachment point
            if let Some(ref ap) = props.attachment_point {
                atom.attachment_point = Some(ap.clone());
            }

            // Apply attachment order
            if let Some(ref ao) = props.attachment_order {
                atom.attachment_order = Some(ao.clone());
            }

            // Apply R-group label
            if let Some(rgroup_label) = props.rgroup_label {
                // Set the symbol to RGroup with a minimal RGroup struct
                // The full RGroup details are handled separately via rgroup_properties
                atom.symbol = AtomSymbol::RGroup(RGroup {
                    label: Some(rgroup_label),
                    dependent_label: None,
                    rgroup_or_h: false,
                    occurrence: vec![],
                });
            }
        }

        // Apply bond properties (zero order bonds)
        for (&bond_idx, props) in &self.bond_properties {
            if let Some(order_code) = props.order_override {
                let extended = flags.contains(CtabParseFlags::EXTENDED_RANGE);
                let queries = flags.contains(CtabParseFlags::QUERIES);
                if let Some(bond) = molecule.bonds.get_mut(bond_idx) {
                    bond.order = convert_extended_bond_type_code(order_code, extended, queries)?;
                }
            }
        }

        // Initialize ctfile_data if we have any CTFile-specific data
        let has_ctfile_data = !self.rgroup_properties.is_empty()
            || !self.sgroup_properties.is_empty();
        if has_ctfile_data {
            if molecule.ctfile_data.is_none() {
                molecule.ctfile_data = Some(CtfileData::default());
            }
        }

        // Apply R-group logic
        if let Some(ref mut ctfile_data) = molecule.ctfile_data {
            for (&rgroup_label, props) in &self.rgroup_properties {
                let rgroup = RGroup {
                    label: Some(rgroup_label as u32),
                    dependent_label: props.dependent_label,
                    rgroup_or_h: props.rgroup_or_h.unwrap_or(false),
                    occurrence: props.occurrence.clone().unwrap_or_default(),
                };
                ctfile_data.rgroups.insert(rgroup_label, rgroup);
            }
        }

        // Validate and apply S-groups
        self.validate_sgroup_data()?;
        self.apply_sgroup(molecule)?;

        Ok(())
    }

    fn apply_sgroup(&mut self, molecule: &mut ExtendedMolecule) -> Result<()> {
        for (sgroup_index, props) in std::mem::take(&mut self.sgroup_properties) {
            let group_type = props.group_type.ok_or_else(|| {
                SemanticError::SGroupMissingType(sgroup_index)
            })?;
            let mut sgroup = SGroup::new(group_type);
            sgroup.label = props.label.or(Some(sgroup_index as u32));
            sgroup.group_subtype = props.group_subtype;
            sgroup.connectivity = props.connectivity;
            if let Some(expansion) = props.expansion {
                sgroup.expansion = expansion;
            }
            if let Some(atom_indices) = props.atom_indices {
                sgroup.atom_indices = atom_indices;
            }
            if let Some(bond_indices) = props.bond_indices {
                sgroup.bond_indices = bond_indices;
            }
            if let Some(parent_atom_indices) = props.parent_atom_indices {
                sgroup.parent_atom_indices = Some(parent_atom_indices);
            }
            if let Some(subscript) = props.subscript {
                sgroup.subscript = Some(subscript);
            }
            if let Some(correspondence) = props.correspondence {
                sgroup.correspondence = Some(correspondence);
            }
            if let Some(connecting_bond) = props.connecting_bond {
                sgroup.connecting_bond = Some(connecting_bond);
            }
            if let Some(hierarchy_parent) = props.hierarchy_parent {
                sgroup.hierarchy_parent = Some(hierarchy_parent);
            }
            if let Some(component_number) = props.component_number {
                sgroup.component_number = Some(component_number);
            }
            if !props.data.is_empty() {
                sgroup.data = props.data;
            }
            if let Some(multiplier) = props.multiplier {
                sgroup.multiplier = Some(multiplier);
            }
            if let Some(bracket_coords) = props.bracket_coords {
                sgroup.bracket_coords = Some(bracket_coords);
            }
            if let Some(display) = props.display {
                sgroup.display = Some(display);
            }
            
            // Store in ctfile_data
            if molecule.ctfile_data.is_none() {
                molecule.ctfile_data = Some(CtfileData::default());
            }
            molecule.ctfile_data.as_mut().unwrap().sgroups.insert(sgroup_index, sgroup);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
