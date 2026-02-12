//! CTFile-specific data for roundtripping.
//!
//! Contains data structures that are specific to CTFile formats (MOL, SDF).
//! This data is preserved for exact roundtripping, scheduled to be replaced
//! by semantically defined structures.

use std::collections::{BTreeMap, HashMap};

use super::rgroup::RGroup;
use super::sgroup::SGroup;

/// Legacy group abbreviation
#[derive(Debug, Clone, PartialEq)]
pub struct LegacyGroupAbbreviation {
    pub atom_index1: u32, // Atoms on this side are abbreviated
    pub atom_index2: u32, // Attachment point to main structure
    pub label: String,
}

/// CTFile-specific data container for ExtendedMolecule
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CtfileData {
    pub sgroups: BTreeMap<u32, SGroup>,
    pub rgroups: BTreeMap<u32, RGroup>,
    pub legacy_group_abbreviations: Vec<LegacyGroupAbbreviation>,
}

impl CtfileData {
    pub fn update_atoms_bonds(
        &self,
        atom_index_map: &HashMap<u32, u32>,
        bond_index_map: &HashMap<u32, u32>,
    ) -> Option<Self> {
        let mut sgroups = BTreeMap::new();
        for (label, sgroup) in &self.sgroups {
            if let Some(remapped) = sgroup.remap_indices(atom_index_map, bond_index_map) {
                sgroups.insert(*label, remapped);
            }
        }

        let legacy_group_abbreviations = self
            .legacy_group_abbreviations
            .iter()
            .filter_map(|entry| {
                let atom_index1 = atom_index_map.get(&entry.atom_index1).copied()?;
                let atom_index2 = atom_index_map.get(&entry.atom_index2).copied()?;
                Some(LegacyGroupAbbreviation {
                    atom_index1,
                    atom_index2,
                    label: entry.label.clone(),
                })
            })
            .collect::<Vec<_>>();

        let data = Self {
            sgroups,
            rgroups: self.rgroups.clone(),
            legacy_group_abbreviations,
        };
        if data.sgroups.is_empty()
            && data.rgroups.is_empty()
            && data.legacy_group_abbreviations.is_empty()
        {
            None
        } else {
            Some(data)
        }
    }
}
