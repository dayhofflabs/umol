//! CXSMILES annotation data for roundtripping.
//!
//! Contains format-specific data that doesn't have clean semantic representation
//! but is needed for faithful roundtripping of CXSMILES.

use std::collections::{BTreeMap, HashMap};

use super::atom::{BicycloStereo, BicycloStereoData};
use super::rgroup::RGroup;
use super::sgroup::SGroup;
use crate::table_ir::atom::Chirality;

/// Local parity entry (@: / @@:). Chiral center with ordered substituents.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalParityCenter {
    pub center: u32,
    pub substituents: Vec<u32>,
    pub chirality: Chirality,
}

/// A set of stereocenters with a common interpretation mode
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoSet {
    pub atoms: Vec<u32>,
    pub relation: StereoSetRelation,
}

/// How to interpret a group of stereocenters
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StereoSetRelation {
    /// Centers are correlated (e.g., racemate-like)
    Correlated,
    /// Centers are independent (e.g., mixture)
    Independent,
}

/// CXSMILES annotation data
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CxAnnotationData {
    /// Enhanced stereo groups: index -> set of atoms with their mode
    pub stereo_groups: BTreeMap<u32, StereoSet>,

    /// Component groupings (atom indices per component)
    /// Used when explicit grouping differs from graph connectivity
    pub components: Option<Vec<Vec<u32>>>,

    /// S-groups from CXSMILES Sg/SgD/SgH tags
    pub sgroups: BTreeMap<u32, SGroup>,

    /// R-groups from CXSMILES LOG: tag (label, occurrence, rgroup_or_h)
    pub rgroups: BTreeMap<u32, RGroup>,

    /// R-group member structures from RG: tag (label -> SMILES strings)
    pub rgroup_members: BTreeMap<u32, Vec<String>>,

    /// Local parity from @: / @@: (chiral center, ordered substituents, chirality)
    pub local_parity: Option<Vec<LocalParityCenter>>,

    /// Bicyclic stereo from THB: / TLB: / TEB:
    pub bicyclo_stereo: Option<Vec<BicycloStereo>>,
}

impl CxAnnotationData {
    pub fn update_atoms_bonds(
        &self,
        atom_index_map: &HashMap<u32, u32>,
        bond_index_map: &HashMap<u32, u32>,
    ) -> Option<Self> {
        let mut stereo_groups = BTreeMap::new();
        for (idx, set) in &self.stereo_groups {
            let atoms = set
                .atoms
                .iter()
                .map(|a| atom_index_map.get(a).copied())
                .collect::<Option<Vec<u32>>>();
            if let Some(atoms) = atoms {
                stereo_groups.insert(
                    *idx,
                    StereoSet {
                        atoms,
                        relation: set.relation,
                    },
                );
            }
        }

        let components = self.components.as_ref().and_then(|groups| {
            let remapped = groups
                .iter()
                .filter_map(|g| {
                    g.iter()
                        .map(|a| atom_index_map.get(a).copied())
                        .collect::<Option<Vec<u32>>>()
                })
                .collect::<Vec<Vec<u32>>>();
            if remapped.is_empty() {
                None
            } else {
                Some(remapped)
            }
        });

        let mut sgroups = BTreeMap::new();
        for (label, sgroup) in &self.sgroups {
            if let Some(remapped) = sgroup.remap_indices(atom_index_map, bond_index_map) {
                sgroups.insert(*label, remapped);
            }
        }

        let local_parity = self
            .local_parity
            .as_ref()
            .map(|entries| {
                entries
                    .iter()
                    .filter_map(|entry| {
                        let center = atom_index_map.get(&entry.center).copied()?;
                        let substituents = entry
                            .substituents
                            .iter()
                            .map(|s| atom_index_map.get(s).copied())
                            .collect::<Option<Vec<u32>>>()?;
                        Some(LocalParityCenter {
                            center,
                            substituents,
                            chirality: entry.chirality,
                        })
                    })
                    .collect::<Vec<LocalParityCenter>>()
            })
            .filter(|entries| !entries.is_empty());

        let bicyclo_stereo = self
            .bicyclo_stereo
            .as_ref()
            .map(|entries| {
                entries
                    .iter()
                    .filter_map(|entry| update_bicyclo_atoms(entry, atom_index_map))
                    .collect::<Vec<BicycloStereo>>()
            })
            .filter(|entries| !entries.is_empty());

        let data = Self {
            stereo_groups,
            components,
            sgroups,
            rgroups: self.rgroups.clone(),
            rgroup_members: self.rgroup_members.clone(),
            local_parity,
            bicyclo_stereo,
        };
        if data.stereo_groups.is_empty()
            && data.components.is_none()
            && data.sgroups.is_empty()
            && data.rgroups.is_empty()
            && data.rgroup_members.is_empty()
            && data.local_parity.is_none()
            && data.bicyclo_stereo.is_none()
        {
            None
        } else {
            Some(data)
        }
    }
}

fn update_bicyclo_atoms(
    entry: &BicycloStereo,
    atom_index_map: &HashMap<u32, u32>,
) -> Option<BicycloStereo> {
    let update_atoms = |data: &BicycloStereoData| -> Option<BicycloStereoData> {
        let ligand_atom = atom_index_map.get(&data.ligand_atom).copied()?;
        let connection_atom = atom_index_map.get(&data.connection_atom).copied()?;
        let lower_bridge_atoms = data
            .lower_bridge_atoms
            .iter()
            .map(|a| atom_index_map.get(a).copied())
            .collect::<Option<Vec<u32>>>()?;
        let higher_bridge_atoms = data
            .higher_bridge_atoms
            .iter()
            .map(|a| atom_index_map.get(a).copied())
            .collect::<Option<Vec<u32>>>()?;
        Some(BicycloStereoData {
            ligand_atom,
            connection_atom,
            lower_bridge_atoms,
            higher_bridge_atoms,
        })
    };

    match entry {
        BicycloStereo::TowardsHigherBridge(data) => {
            Some(BicycloStereo::TowardsHigherBridge(update_atoms(data)?))
        }
        BicycloStereo::TowardsLowerBridge(data) => {
            Some(BicycloStereo::TowardsLowerBridge(update_atoms(data)?))
        }
        BicycloStereo::TowardsEitherBridge(data) => {
            Some(BicycloStereo::TowardsEitherBridge(update_atoms(data)?))
        }
    }
}
