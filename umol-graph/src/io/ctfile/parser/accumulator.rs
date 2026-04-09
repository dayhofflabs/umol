//! Accumulator for molecular properties from a CTAB file

use std::collections::BTreeMap;
use std::mem;

use umol_data::Element;

use super::context::Context;
use super::convert::{
    convert_atom_isotope_mass_number, convert_attachment_point_code, convert_radical_type_code,
    convert_ring_bond_count_code, convert_substitution_count_code, convert_unsaturated_atom_code,
};
use crate::io::ctfile::config::CtabParseFlags;
use crate::io::ctfile::error::ParseError;
use crate::io::ctfile::parser::properties::{PropertyEntries, SGroupDataEntry};
use crate::table_ir::atom::ImplicitHydrogens;
use crate::table_ir::{
    AtomList, AtomSymbol, AttachmentPointType, BondOrder, CtfileData, ExtendedMolecule,
    LegacyGroupAbbreviation, LinkAtom, Molecule, RGroup, RGroupOccurrence, RingBondCount, SGroup,
    SGroupBracketCoords, SGroupConnectingBond, SGroupConnectivity, SGroupData, SGroupDataDisplay,
    SGroupMultiplier, SGroupSubtype, SGroupType, StereoInterpretation, SubstitutionCount,
    UnsaturatedAtom,
};

/// Accumulator for global molecule properties
#[derive(Debug, Default)]
pub(super) struct MoleculeProperties {
    pub chiral_flag: Option<bool>,
}

// Accumulator for properties of a single atom
#[derive(Debug, Default)]
pub(super) struct AtomProperties {
    pub label: Option<String>,
    pub value: Option<String>,
    pub charge: Option<i8>,
    pub unpaired_electrons: Option<u8>,
    pub isotope_mass: Option<u32>,
    pub hydrogen_count: Option<u8>,
    pub pattern: Option<String>,
    pub ring_bond_count: Option<RingBondCount>,
    pub substitution_count: Option<SubstitutionCount>,
    pub unsaturated: Option<UnsaturatedAtom>,
    pub link_atom: Option<LinkAtom>,
    pub attachment_point: Option<AttachmentPointType>,
    pub attachment_order: Option<Vec<(u32, u8)>>,
    pub rgroup_label: Option<u32>,
    pub atom_list_elements: Option<Vec<Element>>,
    pub atom_list_exclusion: Option<bool>,
}

// Accumulator for properties of a single bond
#[derive(Debug, Default)]
pub(super) struct BondProperties {
    pub order_override: Option<BondOrder>,
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
    pub atom_indices: Option<Vec<u32>>,
    pub bond_indices: Option<Vec<u32>>,
    pub parent_atom_indices: Option<Vec<u32>>,
    pub multiplier: Option<SGroupMultiplier>,
    pub subscript: Option<String>,
    pub correspondence: Option<Vec<u32>>,
    pub bracket_coords: Option<SGroupBracketCoords>,
    pub connecting_bond: Option<SGroupConnectingBond>,
    pub hierarchy_parent: Option<u32>,
    pub component_number: Option<u32>,
    pub data: Option<SGroupData>,
    pub display: Option<SGroupDataDisplay>,
}

impl SGroupProperties {
    pub(crate) fn new(sgroup_type: SGroupType) -> Self {
        Self {
            group_type: Some(sgroup_type),
            group_subtype: None,
            label: None,
            connectivity: None,
            expansion: None,
            atom_indices: None,
            bond_indices: None,
            parent_atom_indices: None,
            multiplier: None,
            subscript: None,
            correspondence: None,
            bracket_coords: None,
            connecting_bond: None,
            hierarchy_parent: None,
            component_number: None,
            data: None,
            display: None,
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
pub(super) struct PropertyAccumulator {
    context: Context,
    pub molecule_properties: Vec<MoleculeProperties>,
    pub atom_properties: BTreeMap<u32, AtomProperties>,
    pub bond_properties: BTreeMap<u32, BondProperties>,
    pub rgroup_properties: BTreeMap<u32, RGroupProperties>,
    pub sgroup_properties: BTreeMap<u32, SGroupProperties>,
    pub legacy_group_abbreviations: Vec<LegacyGroupAbbreviation>,
}

impl PropertyAccumulator {
    pub(crate) fn new() -> Self {
        Self {
            context: Context::new(),
            molecule_properties: Vec::new(),
            atom_properties: BTreeMap::new(),
            bond_properties: BTreeMap::new(),
            rgroup_properties: BTreeMap::new(),
            sgroup_properties: BTreeMap::new(),
            legacy_group_abbreviations: Vec::new(),
        }
    }

    pub(crate) fn add_entry(
        &mut self,
        entry: PropertyEntries,
        flags: CtabParseFlags,
    ) -> Result<(), ParseError> {
        let extended_range = flags.contains(CtabParseFlags::EXTENDED_RANGE);
        match entry {
            PropertyEntries::MoleculeChiralFlagEntry(e) => {
                self.molecule_properties.push(MoleculeProperties {
                    chiral_flag: Some(e.chiral_flag),
                });
            }
            PropertyEntries::AtomAliasEntry(e) => {
                let props = self.atom_properties.entry(e.atom_index).or_default();
                props.label = Some(e.alias);
            }
            PropertyEntries::LegacyGroupAbbreviationEntry(e) => {
                self.legacy_group_abbreviations
                    .push(LegacyGroupAbbreviation {
                        atom_index1: e.atom_index1,
                        atom_index2: e.atom_index2,
                        label: e.label,
                    });
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
                    props.unpaired_electrons = convert_radical_type_code(entry.radical_type)?;
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
                        return Err(ParseError::DuplicateProperty(format!(
                            "link atom for atom {}",
                            entry.atom_index
                        )));
                    }
                    props.link_atom = Some(LinkAtom {
                        min_repeat: 0,
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
                        return Err(ParseError::DuplicateProperty(format!(
                            "attachment point for atom {}",
                            entry.atom_index
                        )));
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
                let props = self.rgroup_properties.entry(entry.label).or_default();
                props.dependent_label = entry.dependent_label;
                props.rgroup_or_h = Some(entry.rgroup_or_h);
                props.occurrence = Some(entry.occurrence);
            }
            PropertyEntries::SGroupTypeEntries(entries) => {
                for entry in entries {
                    if self.sgroup_properties.contains_key(&entry.sgroup_index) {
                        return Err(ParseError::DuplicateProperty(format!(
                            "S-group type for index {}",
                            entry.sgroup_index
                        )));
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
                    let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                        ParseError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "subtype",
                        },
                    )?;
                    props.group_subtype = Some(entry.sgroup_subtype);
                }
            }
            PropertyEntries::SGroupLabelEntries(entries) => {
                for entry in entries {
                    let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                        ParseError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "label",
                        },
                    )?;
                    props.label = Some(entry.label);
                }
            }
            PropertyEntries::SGroupConnectivityEntries(entries) => {
                for entry in entries {
                    let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                        ParseError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "connectivity",
                        },
                    )?;
                    props.connectivity = Some(entry.connectivity);
                }
            }
            PropertyEntries::SGroupExpansionEntries(entries) => {
                for entry in entries {
                    let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                        ParseError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "expansion",
                        },
                    )?;
                    props.expansion = Some(true);
                }
            }
            PropertyEntries::SGroupAtomListEntry(entry) => {
                let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                    ParseError::UndefinedSGroup {
                        index: entry.sgroup_index,
                        property: "atom list",
                    },
                )?;
                props.atom_indices = Some(entry.atom_indices);
            }
            PropertyEntries::SGroupBondListEntry(entry) => {
                let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                    ParseError::UndefinedSGroup {
                        index: entry.sgroup_index,
                        property: "bond list",
                    },
                )?;
                props.bond_indices = Some(entry.bond_indices);
            }
            PropertyEntries::SGroupParentAtomEntry(entry) => {
                let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                    ParseError::UndefinedSGroup {
                        index: entry.sgroup_index,
                        property: "parent atom",
                    },
                )?;
                props.parent_atom_indices = Some(entry.atom_indices);
            }
            PropertyEntries::SGroupSubscriptEntry(entry) => {
                let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                    ParseError::UndefinedSGroup {
                        index: entry.sgroup_index,
                        property: "subscript",
                    },
                )?;
                let sgroup_type = self.context.sgroup_types.get(&entry.sgroup_index).ok_or(
                    ParseError::UndefinedSGroup {
                        index: entry.sgroup_index,
                        property: "subscript",
                    },
                )?;

                if !sgroup_accepts_subscript(*sgroup_type)
                    && !sgroup_accepts_multiplier(*sgroup_type)
                {
                    return Err(ParseError::SGroupTypeConstraint {
                        sgroup_type: *sgroup_type,
                        message: "subscript and multiplier not allowed",
                    });
                } else if sgroup_accepts_subscript(*sgroup_type) {
                    if entry.subscript.is_none() {
                        return Err(ParseError::SGroupTypeConstraint {
                            sgroup_type: *sgroup_type,
                            message: "subscript required",
                        });
                    }
                    props.subscript = entry.subscript;
                } else if sgroup_accepts_multiplier(*sgroup_type) {
                    if entry.multiplier.is_none() {
                        return Err(ParseError::SGroupTypeConstraint {
                            sgroup_type: *sgroup_type,
                            message: "multiplier required",
                        });
                    }
                    props.multiplier = entry.multiplier;
                } else {
                    return Err(ParseError::SGroupTypeConstraint {
                        sgroup_type: *sgroup_type,
                        message: "cannot have both subscript and multiplier",
                    });
                }
            }
            PropertyEntries::SGroupCorrespondenceEntry(entry) => {
                let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                    ParseError::UndefinedSGroup {
                        index: entry.sgroup_index,
                        property: "correspondence",
                    },
                )?;
                props.correspondence = Some(entry.bond_indices);
            }
            PropertyEntries::SGroupDisplayInfoEntry(entry) => {
                let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                    ParseError::UndefinedSGroup {
                        index: entry.sgroup_index,
                        property: "display info",
                    },
                )?;
                props.bracket_coords = Some(SGroupBracketCoords {
                    bracket1: (entry.bracket_coords[0], entry.bracket_coords[1]),
                    bracket2: (entry.bracket_coords[2], entry.bracket_coords[3]),
                    bracket3: None,
                    bracket4: None,
                });
            }
            PropertyEntries::SGroupConnectingBondEntry(entry) => {
                let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                    ParseError::UndefinedSGroup {
                        index: entry.sgroup_index,
                        property: "connecting bond",
                    },
                )?;
                props.connecting_bond = Some(SGroupConnectingBond {
                    bond_index: entry.bond_index,
                    bond_vector: entry.bond_vector,
                });
            }
            PropertyEntries::SGroupDataDescriptionEntry(entry) => {
                let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                    ParseError::UndefinedSGroup {
                        index: entry.sgroup_index,
                        property: "data description",
                    },
                )?;
                let field_name = entry.field_name.clone();
                props.data = Some(SGroupData {
                    field_type: entry.field_type,
                    field_name: entry.field_name,
                    field_units: entry.field_units,
                    query_identifier: entry.query_identifier,
                    data_query_operator: entry.data_query_operator,
                    data_content: None,
                });

                self.context.current_data_sgroup_index = Some(entry.sgroup_index);
                self.context.current_data_field = Some(field_name);
            }
            PropertyEntries::SGroupDataDisplayEntry(entry) => {
                let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                    ParseError::UndefinedSGroup {
                        index: entry.sgroup_index,
                        property: "data description",
                    },
                )?;
                props.display = Some(SGroupDataDisplay {
                    coords: entry.coords,
                    display_type: entry.display_type,
                    display_placement: entry.display_placement,
                    display_units: entry.display_units,
                    display_chars: entry.display_chars,
                });
            }
            // TODO: Verify that multiple records can be attached to a single data SGroup.
            PropertyEntries::SGroupDataEntry(entry) => {
                // NOTE: SGroup Data Description (SDT) field must appear before the SCD/SED fields
                // SDT sets self.context.current_data_sgroup_index and self.context.current_data_field
                let sgroup_index = match entry {
                    SGroupDataEntry::Continuation {
                        sgroup_index: sg, ..
                    } => sg,
                    SGroupDataEntry::EndBlank { sgroup_index: sg } => sg,
                    SGroupDataEntry::EndWithData {
                        sgroup_index: sg, ..
                    } => sg,
                };
                if let Some(context_sgroup) = self.context.current_data_sgroup_index {
                    debug_assert!(
                        self.context.current_data_field.is_some(),
                        "current data field must be set together with current data sgroup index"
                    );

                    if sgroup_index != context_sgroup {
                        return Err(ParseError::SGroupIndexMismatch {
                            expected: context_sgroup,
                            actual: sgroup_index,
                        });
                    }
                } else {
                    let location = match entry {
                        SGroupDataEntry::Continuation { .. } => "continuation",
                        SGroupDataEntry::EndBlank { .. } => "end blank",
                        SGroupDataEntry::EndWithData { .. } => "end with data",
                    };
                    return Err(ParseError::MissingSGroupDataContext {
                        index: sgroup_index,
                        location,
                    });
                }

                match entry {
                    SGroupDataEntry::Continuation {
                        sgroup_index: _sg,
                        data_content,
                    } => {
                        if let Some(context_content) = self.context.current_data_content.as_mut() {
                            context_content.push(data_content);
                        } else {
                            self.context.current_data_content = Some(vec![data_content]);
                        }
                    }
                    SGroupDataEntry::EndWithData {
                        sgroup_index,
                        data_content,
                    } => {
                        if let Some(context_content) = self.context.current_data_content.as_mut() {
                            context_content.push(data_content);
                        } else {
                            self.context.current_data_content = Some(vec![data_content]);
                        }
                        self.finalize_sgroup_data(sgroup_index)?;
                    }
                    SGroupDataEntry::EndBlank { sgroup_index } => {
                        if self.context.current_data_content.is_none() {
                            return Err(ParseError::MissingSGroupDataContext {
                                index: sgroup_index,
                                location: "end",
                            });
                        }
                        self.finalize_sgroup_data(sgroup_index)?;
                    }
                }
            }
            PropertyEntries::SGroupHierarchyEntries(entries) => {
                for entry in entries {
                    if !self
                        .sgroup_properties
                        .contains_key(&entry.parent_sgroup_index)
                    {
                        return Err(ParseError::UndefinedSGroup {
                            index: entry.parent_sgroup_index,
                            property: "hierarchy parent",
                        });
                    }
                    let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                        ParseError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "hierarchy",
                        },
                    )?;
                    props.hierarchy_parent = Some(entry.parent_sgroup_index);
                }
            }
            PropertyEntries::SGroupComponentEntries(entries) => {
                for entry in entries {
                    let props = self.sgroup_properties.get_mut(&entry.sgroup_index).ok_or(
                        ParseError::UndefinedSGroup {
                            index: entry.sgroup_index,
                            property: "component",
                        },
                    )?;
                    props.component_number = Some(entry.component_number);
                }
            }
            PropertyEntries::BondOrderOverrideEntries(entries) => {
                for entry in entries {
                    self.bond_properties.insert(
                        entry.bond_index,
                        BondProperties {
                            order_override: Some(entry.bond_order),
                        },
                    );
                }
            }
            PropertyEntries::AtomChargeOverrideEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    props.charge = Some(entry.charge);
                }
            }
            PropertyEntries::AtomHydrogenCountEntries(entries) => {
                for entry in entries {
                    let props = self.atom_properties.entry(entry.atom_index).or_default();
                    if entry.hydrogen_count.is_some() {
                        props.hydrogen_count = entry.hydrogen_count;
                    }
                }
            }
            PropertyEntries::ChemSketchLabelEntry(e) => {
                let props = self.atom_properties.entry(e.atom_index).or_default();
                props.label = Some(e.label);
            }
            PropertyEntries::MarvinSmartsPatternEntry(e) => {
                let props = self.atom_properties.entry(e.atom_index).or_default();
                props.pattern = Some(e.smarts_pattern);
            }
        }

        Ok(())
    }

    fn validate_sgroup_data(&mut self) -> Result<(), ParseError> {
        if self.context.current_data_sgroup_index.is_some()
            || self.context.current_data_field.is_some()
        {
            return Err(ParseError::MissingSgroupDataEnd {
                index: self.context.current_data_sgroup_index.unwrap(),
            });
        }
        Ok(())
    }

    fn finalize_sgroup_data(&mut self, sgroup_index: u32) -> Result<(), ParseError> {
        if let Some(context_sgroup) = self.context.current_data_sgroup_index {
            debug_assert!(
                self.context.current_data_field.is_some(),
                "current data field must be set together with current data sgroup index"
            );
            debug_assert!(
                self.context.current_data_content.is_some(),
                "current data context must be set together with current data sgroup index"
            );
            if sgroup_index != context_sgroup {
                return Err(ParseError::SGroupIndexMismatch {
                    expected: context_sgroup,
                    actual: sgroup_index,
                });
            }
        } else {
            return Err(ParseError::MissingSGroupDataContext {
                index: sgroup_index,
                location: "finalization",
            });
        }

        let props =
            self.sgroup_properties
                .get_mut(&sgroup_index)
                .ok_or(ParseError::UndefinedSGroup {
                    index: sgroup_index,
                    property: "data content",
                })?;

        let mut content = self.context.current_data_content.as_mut().unwrap().join("");
        if content.len() > 200 {
            content.truncate(200);
        }
        let content = content.trim_end().to_string();

        if props.data.is_none() {
            return Err(ParseError::MissingSGroupDataContext {
                index: sgroup_index,
                location: "finalization",
            });
        }

        let data = props.data.as_mut().unwrap();
        if let Some(data_content) = data.data_content.as_mut() {
            data_content.push(content);
        } else {
            data.data_content = Some(vec![content]);
        }

        self.context.current_data_sgroup_index = None;
        self.context.current_data_field = None;
        self.context.current_data_content = None;
        Ok(())
    }

    /// Apply all properties to Molecule
    pub(crate) fn update_molecule(
        &mut self,
        molecule: &mut Molecule,
        flags: CtabParseFlags,
    ) -> Result<(), ParseError> {
        let extended_isotopes = flags.contains(CtabParseFlags::EXTENDED_ISOTOPES);

        // Apply molecule properties (chiral flag)
        if let Some(chiral_flag) = self.molecule_properties.first().and_then(|p| p.chiral_flag) {
            molecule
                .properties
                .insert("chiral_flag".to_string(), chiral_flag.to_string());
            if chiral_flag {
                molecule.stereo_interpretation = Some(StereoInterpretation::Absolute);
            }
        }

        // Apply atom properties (only basic ones compatible with Atom)
        for (&atom_idx, props) in &self.atom_properties {
            let Some(atom) = molecule.atoms.get_mut(atom_idx as usize) else {
                return Err(ParseError::IndexOutOfBounds(atom_idx));
            };

            // Apply alias
            if let Some(ref label) = props.label {
                atom.label = Some(label.clone());
            }

            // Apply value
            if let Some(ref value) = props.value {
                atom.value = Some(value.clone());
            }

            // Apply charge
            if let Some(charge) = props.charge {
                atom.charge = Some(charge);
                atom.unpaired_electrons = None; // charge overrides radical from atom block
                atom.multiplicity = None;
            }

            // Apply radical (unpaired electrons)
            if let Some(unpaired_electrons) = props.unpaired_electrons {
                atom.unpaired_electrons = Some(unpaired_electrons);
                atom.multiplicity = None;
                atom.charge = None; // radical overrides charge from atom block
            }

            // Apply isotope
            if let Some(isotope) = props.isotope_mass {
                let mass =
                    convert_atom_isotope_mass_number(atom.element, isotope, extended_isotopes)?;
                atom.isotope_mass = mass;
            }

            // Apply hydrogen count
            if let Some(h_count) = props.hydrogen_count {
                atom.implicit_hydrogens = Some(ImplicitHydrogens::Hydrogens(h_count));
            }
        }

        // Apply bond properties (bond order override)
        for (&bond_idx, props) in &self.bond_properties {
            if let Some(bo) = props.order_override {
                let Some(bond) = molecule.bonds.get_mut(bond_idx as usize) else {
                    return Err(ParseError::IndexOutOfBounds(bond_idx));
                };
                bond.order = bo;
            }
        }

        Ok(())
    }

    pub(crate) fn update_extended_molecule(
        &mut self,
        molecule: &mut ExtendedMolecule,
        flags: CtabParseFlags,
    ) -> Result<(), ParseError> {
        let extended_isotopes = flags.contains(CtabParseFlags::EXTENDED_ISOTOPES);

        // Apply molecule properties (chiral flag)
        if let Some(chiral_flag) = self.molecule_properties.first().and_then(|p| p.chiral_flag) {
            molecule
                .properties
                .insert("chiral_flag".to_string(), chiral_flag.to_string());
            if chiral_flag {
                molecule.stereo_interpretation = Some(StereoInterpretation::Absolute);
            }
        }

        // Apply atom properties
        for (&atom_idx, props) in &self.atom_properties {
            let Some(atom) = molecule.atoms.get_mut(atom_idx as usize) else {
                return Err(ParseError::IndexOutOfBounds(atom_idx));
            };

            // Apply alias
            if let Some(ref label) = props.label {
                atom.label = Some(label.clone());
            }

            // Apply value
            if let Some(ref value) = props.value {
                atom.value = Some(value.clone());
            }

            // Apply charge
            if let Some(charge) = props.charge {
                atom.charge = Some(charge);
                atom.unpaired_electrons = None; // charge overrides radical from atom block
                atom.multiplicity = None;
            }

            // Apply radical (unpaired electrons)
            if let Some(unpaired_electrons) = props.unpaired_electrons {
                atom.unpaired_electrons = Some(unpaired_electrons);
                atom.multiplicity = None;
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
                        atom.isotope_mass = Some(isotope);
                    }
                }
            }

            // Apply hydrogen count
            if let Some(h_count) = props.hydrogen_count {
                atom.implicit_hydrogens = Some(ImplicitHydrogens::Hydrogens(h_count));
            }

            // Apply SMARTS pattern
            if let Some(ref pattern) = props.pattern {
                if atom.pattern.is_some() {
                    return Err(ParseError::DuplicateProperty(format!(
                        "SMARTS pattern conflict: existing value for atom {}",
                        atom_idx
                    )));
                }
                atom.pattern = Some(pattern.clone());
            }

            // Apply ring bond count
            if let Some(rbc) = props.ring_bond_count {
                if atom.ring_bond_count.is_some() {
                    return Err(ParseError::DuplicateProperty(format!(
                        "Ring bond count conflict: existing value for atom {}",
                        atom_idx
                    )));
                }
                atom.ring_bond_count = Some(rbc);
            }

            // Apply substitution count
            if let Some(sub) = props.substitution_count {
                if atom.substitution_count.is_some() {
                    return Err(ParseError::DuplicateProperty(format!(
                        "Substitution count conflict: existing value for atom {}",
                        atom_idx
                    )));
                }
                atom.substitution_count = Some(sub);
            }

            // Apply unsaturated
            if props.unsaturated.is_some() {
                if atom.unsaturated.is_some() {
                    return Err(ParseError::DuplicateProperty(format!(
                        "Unsaturated conflict: existing value for atom {}",
                        atom_idx
                    )));
                }
                atom.unsaturated = Some(UnsaturatedAtom);
            }

            // Apply link atom
            if let Some(link) = props.link_atom {
                if atom.link_atom.is_some() {
                    return Err(ParseError::DuplicateProperty(format!(
                        "Link atom conflict: existing value for atom {}",
                        atom_idx
                    )));
                }
                atom.link_atom = Some(link);
            }

            // Apply atom list
            if let Some(ref elements) = props.atom_list_elements {
                let exclusion = props.atom_list_exclusion.unwrap_or(false);
                // Allow overwriting Element, NamedIsotope, or empty AtomList (L placeholder)
                let can_set = match &atom.symbol {
                    AtomSymbol::Element(_) | AtomSymbol::NamedIsotope(_) => true,
                    AtomSymbol::AtomList(list) => list.elements.is_empty(),
                    _ => false,
                };
                if !can_set {
                    return Err(ParseError::DuplicateProperty(format!(
                        "Atom list conflict: existing symbol for atom {}",
                        atom_idx
                    )));
                }
                atom.symbol = AtomSymbol::AtomList(AtomList {
                    elements: elements.clone(),
                    exclusion,
                });
            }

            // Apply attachment point
            if let Some(ap) = props.attachment_point {
                atom.attachment_point = Some(ap);
            }

            // Apply attachment order
            if let Some(ref ao) = props.attachment_order {
                atom.attachment_order = Some(ao.clone());
            }

            // Apply R-group label
            if let Some(rgroup_label) = props.rgroup_label {
                // Check for conflicts: can't set RGroup label on AtomList
                if matches!(atom.symbol, AtomSymbol::AtomList(_)) {
                    return Err(ParseError::DuplicateProperty(format!(
                        "RGroup label conflict: atom {} already has AtomList",
                        atom_idx
                    )));
                }
                // Check for conflicts: can't change existing RGroup label (but can overwrite None)
                if let AtomSymbol::RGroup(existing_rgroup) = &atom.symbol {
                    if let Some(existing_label) = existing_rgroup.label {
                        if existing_label != rgroup_label {
                            return Err(ParseError::DuplicateProperty(
                                format!("RGroup label conflict: atom {} already has RGroup label {}, cannot set {}", atom_idx, existing_label, rgroup_label),
                            ));
                        }
                    }
                }
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
            if let Some(bo) = props.order_override {
                let Some(bond) = molecule.bonds.get_mut(bond_idx as usize) else {
                    return Err(ParseError::IndexOutOfBounds(bond_idx));
                };
                bond.order = bo;
            }
        }

        // Initialize ctfile_data if we have any CTFile-specific data
        if (!self.rgroup_properties.is_empty()
            || !self.sgroup_properties.is_empty()
            || !self.legacy_group_abbreviations.is_empty())
            && molecule.ctfile_data.is_none()
        {
            molecule.ctfile_data = Some(CtfileData::default());
        }

        // Apply R-group logic and legacy group abbreviations
        if let Some(ref mut ctfile_data) = molecule.ctfile_data {
            for (&rgroup_label, props) in &self.rgroup_properties {
                let rgroup = RGroup {
                    label: Some(rgroup_label),
                    dependent_label: props.dependent_label,
                    rgroup_or_h: props.rgroup_or_h.unwrap_or(false),
                    occurrence: props.occurrence.clone().unwrap_or_default(),
                };
                ctfile_data.rgroups.insert(rgroup_label, rgroup);
            }
            // TODO: Add verifiication for atom indices
            ctfile_data
                .legacy_group_abbreviations
                .append(&mut self.legacy_group_abbreviations);
        }

        // Validate and apply S-groups
        self.validate_sgroup_data()?;
        self.apply_sgroup(molecule)?;

        Ok(())
    }

    fn apply_sgroup(&mut self, molecule: &mut ExtendedMolecule) -> Result<(), ParseError> {
        for (sgroup_index, props) in mem::take(&mut self.sgroup_properties) {
            debug_assert!(
                molecule.ctfile_data.is_some(),
                "ctfile_data should be present"
            );

            let group_type = props
                .group_type
                .ok_or(ParseError::SGroupMissingType(sgroup_index))?;
            let mut sgroup = SGroup::new(group_type);
            sgroup.label = props.label.or(Some(sgroup_index));
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
            if let Some(data) = props.data {
                sgroup.data = Some(data);
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
            molecule
                .ctfile_data
                .as_mut()
                .unwrap()
                .sgroups
                .insert(sgroup_index, sgroup);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
