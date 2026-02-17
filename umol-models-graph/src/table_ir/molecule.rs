//! Molecule types for TableIR.
//!
//! `Molecule` is TableIR molecule representation for molecules with fixed composition
//! `ExtendedMolecule` is temporary container for generalized molecules (query molecules, molecule
//!   libraries, substitutions, polymers, etc.). Currently used to hold extensions from CTFile formats
//!   (Query features, SGroups, RGroups, etc.). It is supposed to be split into multiple semantically
//!   defined structures.

use std::collections::{BTreeMap, HashMap, HashSet};

use indexmap::IndexMap;
use itertools::Itertools;
use petgraph::unionfind::UnionFind;
use umol_data::Element;

use super::atom::{Atom, AtomSymbol, ExtendedAtom};
use super::bond::{Bond, ExtendedBond};
use super::ctfile_data::CtfileData;
use super::cx_data::CxAnnotationData;
use super::error::JoinError;
use super::multicenter::MulticenterBond;
use super::rgroup::RGroup;
use super::sgroup::SGroup;
use super::source::SourceFormat;
use super::stereo::StereoInterpretation;
use super::topology::Ring;
use super::utils::{element_symbol_key, format_sum_formula};
use crate::position::Point3D;

/// Basic molecule IR
#[derive(Clone, Debug, PartialEq)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub rings: Vec<Ring>,
    pub positions: Option<Vec<Point3D>>,
    pub multicenter_bonds: Vec<MulticenterBond>,
    pub stereo_interpretation: Option<StereoInterpretation>,
    pub comments: Vec<String>,
    pub properties: IndexMap<String, String>,
    pub source_format: SourceFormat,
}

impl Molecule {
    pub fn empty() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            rings: Vec::new(),
            positions: None,
            multicenter_bonds: Vec::new(),
            stereo_interpretation: None,
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::UNKNOWN,
        }
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    pub fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.len()
    }

    /// Count of molecular properties (from SDF/MOL properties or CXSMILES annotations)
    pub fn property_count(&self) -> usize {
        self.properties.len()
    }

    /// Component labels per atom
    fn component_labels(&self) -> Vec<u32> {
        let mut union_find = UnionFind::new(self.atom_count());
        for bond in &self.bonds {
            let (first, second) = bond.atoms.as_tuple();
            union_find.union(first, second);
        }
        for multicenter_bond in &self.multicenter_bonds {
            let indices: Vec<u32> = multicenter_bond
                .contributions()
                .iter()
                .flat_map(|c| c.atoms().iter().copied())
                .collect();
            for (first, second) in indices.into_iter().tuple_combinations() {
                union_find.union(first, second);
            }
        }
        union_find.into_labeling()
    }

    /// Count of connected components in the molecule
    pub fn component_count(&self) -> u32 {
        self.component_labels().iter().collect::<HashSet<_>>().len() as u32
    }

    /// Atom indices in each connected component
    pub fn component_atom_indices(&self) -> Vec<Vec<u32>> {
        let mut components = BTreeMap::new();
        for (atom_index, component_label) in self.component_labels().into_iter().enumerate() {
            components
                .entry(component_label)
                .or_insert_with(Vec::new)
                .push(atom_index as u32);
        }
        components.into_values().collect()
    }

    /// Split molecule into connected components with remapped atom indices.
    pub fn split_components(&self) -> Vec<Molecule> {
        let mut components = Vec::new();
        for component_atoms in self.component_atom_indices() {
            let index_map: HashMap<u32, u32> = component_atoms
                .iter()
                .copied()
                .enumerate()
                .map(|(new_idx, old_idx)| (old_idx, new_idx as u32))
                .collect();

            let atoms = component_atoms
                .iter()
                .map(|&old_idx| self.atoms[old_idx as usize].clone())
                .collect();

            let bonds = self
                .bonds
                .iter()
                .filter_map(|bond| {
                    let (a, b) = bond.atoms.as_tuple();
                    let new_a = *index_map.get(&a)?;
                    let new_b = *index_map.get(&b)?;
                    Some(bond.update_atoms(new_a, new_b))
                })
                .collect();

            let rings =
                self.rings
                    .iter()
                    .filter_map(|ring| {
                        let new_start = match ring.start_atom {
                            Some(old_idx) => Some(*index_map.get(&old_idx)?),
                            None => None,
                        };
                        let new_end = match ring.end_atom {
                            Some(old_idx) => Some(*index_map.get(&old_idx)?),
                            None => None,
                        };
                        if new_start.is_none() && new_end.is_none() {
                            return None;
                        }
                        Some(ring.update_atoms(
                            new_start.unwrap_or_default(),
                            new_end.unwrap_or_default(),
                        ))
                    })
                    .collect();

            let positions = self.positions.as_ref().map(|positions| {
                component_atoms
                    .iter()
                    .filter_map(|&old_idx| positions.get(old_idx as usize).copied())
                    .collect()
            });

            let multicenter_bonds = self
                .multicenter_bonds
                .iter()
                .filter_map(|multicenter| multicenter.update_atoms(&index_map))
                .collect();

            components.push(Molecule {
                atoms,
                bonds,
                rings,
                positions,
                multicenter_bonds,
                stereo_interpretation: self.stereo_interpretation,
                comments: self.comments.clone(),
                properties: self.properties.clone(),
                source_format: self.source_format,
            });
        }
        components
    }

    /// Combine a list of molecules into one molecule.
    pub fn join_components(components: &[Molecule]) -> Molecule {
        if components.is_empty() {
            return Molecule::empty();
        }

        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let mut rings = Vec::new();
        let mut multicenter_bonds = Vec::new();
        let mut comments = Vec::new();
        let mut properties = IndexMap::new();
        let mut all_positions = true;
        let mut positions_acc = Vec::new();
        let source_format = components[0].source_format;
        let stereo_interpretation = components[0].stereo_interpretation;

        for mol in components {
            let atom_offset = atoms.len() as u32;
            atoms.extend(mol.atoms.iter().cloned());

            for bond in &mol.bonds {
                let (a, b) = bond.atoms.as_tuple();
                bonds.push(bond.update_atoms(a + atom_offset, b + atom_offset));
            }

            for ring in &mol.rings {
                let new_start = ring.start_atom.map(|idx| idx + atom_offset);
                let new_end = ring.end_atom.map(|idx| idx + atom_offset);
                if new_start.is_none() && new_end.is_none() {
                    continue;
                }
                rings.push(
                    ring.update_atoms(new_start.unwrap_or_default(), new_end.unwrap_or_default()),
                );
            }

            let index_map: HashMap<u32, u32> = (0..mol.atom_count() as u32)
                .map(|idx| (idx, idx + atom_offset))
                .collect();
            multicenter_bonds.extend(
                mol.multicenter_bonds
                    .iter()
                    .filter_map(|multicenter| multicenter.update_atoms(&index_map)),
            );

            if let Some(component_positions) = mol.positions.as_ref() {
                positions_acc.extend(component_positions.iter().copied());
            } else {
                all_positions = false;
            }

            comments.extend(mol.comments.iter().cloned());
            for (key, value) in &mol.properties {
                properties.insert(key.clone(), value.clone());
            }
        }

        Molecule {
            atoms,
            bonds,
            rings,
            positions: if all_positions {
                Some(positions_acc)
            } else {
                None
            },
            multicenter_bonds,
            stereo_interpretation,
            comments,
            properties,
            source_format,
        }
    }

    /// Get sum formula in Hill notation (C first, H second, then alphabetically)
    pub fn sum_formula(&self) -> String {
        let mut atom_counts: BTreeMap<[u8; 2], (Element, usize)> = BTreeMap::new();
        let mut c_count = 0usize;
        let mut h_count = 0usize;
        let mut charge = 0i32;

        for atom in &self.atoms {
            let element = atom.element;
            match element {
                Element::C => c_count += 1,
                Element::H => h_count += 1,
                e => {
                    let key = element_symbol_key(e);
                    atom_counts.entry(key).or_insert((e, 0)).1 += 1;
                }
            }
            if let Some(ch) = atom.charge {
                charge += ch as i32;
            }
        }

        format_sum_formula(c_count, h_count, atom_counts, charge)
    }
}

/// Extended molecule IR - includes MDL extensions (SGroups, RGroups, etc.)
/// This is a flat structure using ExtendedAtom and ExtendedBond.
#[derive(Clone, Debug, PartialEq)]
pub struct ExtendedMolecule {
    pub atoms: Vec<ExtendedAtom>,
    pub bonds: Vec<ExtendedBond>,
    pub rings: Vec<Ring>,
    pub positions: Option<Vec<Point3D>>,
    pub multicenter_bonds: Vec<MulticenterBond>,

    pub stereo_interpretation: Option<StereoInterpretation>,
    pub comments: Vec<String>,
    pub properties: IndexMap<String, String>,

    pub ctfile_data: Option<CtfileData>,
    pub cx_data: Option<CxAnnotationData>,

    pub source_format: SourceFormat,
}

impl ExtendedMolecule {
    pub fn empty() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            rings: Vec::new(),
            positions: None,
            multicenter_bonds: Vec::new(),
            stereo_interpretation: None,
            comments: Vec::new(),
            properties: IndexMap::new(),
            ctfile_data: None,
            cx_data: None,
            source_format: SourceFormat::UNKNOWN,
        }
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    pub fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.len()
    }

    /// Get sum formula in Hill notation (C first, H second, then alphabetically)
    /// Extended atoms (wildcards, atom lists, R-groups, pseudoatoms, lone pairs) are appended
    /// after elements using their symbolic representation.
    pub fn sum_formula(&self) -> String {
        let mut atom_counts: BTreeMap<[u8; 2], (Element, usize)> = BTreeMap::new();
        let mut c_count = 0usize;
        let mut h_count = 0usize;
        let mut charge = 0i32;

        // Count extended atom types (in enum variant order)
        let mut wildcard_count = 0usize;
        let mut atomlist_count = 0usize;
        let mut rgroup_count = 0usize;
        let mut pseudoatom_count = 0usize;
        let mut lonepair_count = 0usize;

        for atom in &self.atoms {
            match &atom.symbol {
                AtomSymbol::Element(e) => match e {
                    Element::C => c_count += 1,
                    Element::H => h_count += 1,
                    e => {
                        let key = element_symbol_key(*e);
                        atom_counts.entry(key).or_insert((*e, 0)).1 += 1;
                    }
                },
                AtomSymbol::NamedIsotope(i) => match i.element() {
                    Element::C => c_count += 1,
                    Element::H => h_count += 1,
                    e => {
                        let key = element_symbol_key(e);
                        atom_counts.entry(key).or_insert((e, 0)).1 += 1;
                    }
                },
                AtomSymbol::WildcardAtom(_) => wildcard_count += 1,
                AtomSymbol::AtomList(_) => atomlist_count += 1,
                AtomSymbol::RGroup(_) => rgroup_count += 1,
                AtomSymbol::Pseudoatom(_) => pseudoatom_count += 1,
                AtomSymbol::LonePair => lonepair_count += 1,
            }
            if let Some(ch) = atom.charge {
                charge += ch as i32;
            }
        }

        let mut result = format_sum_formula(c_count, h_count, atom_counts, charge);

        // Append extended atom counts (symbol, count with index 1 elided)
        let extended = [
            ("*", wildcard_count),
            ("[L]", atomlist_count),
            ("R", rgroup_count),
            ("[Ps]", pseudoatom_count),
            ("[LP]", lonepair_count),
        ];
        for (symbol, count) in extended {
            if count > 1 {
                result.push_str(&format!("{}{}", symbol, count));
            } else if count == 1 {
                result.push_str(symbol);
            }
        }

        result
    }

    /// Extract basic molecule (converts ExtendedAtom/ExtendedBond to basic types)
    pub fn to_molecule(&self) -> Result<Molecule, super::error::ConversionError> {
        Ok(Molecule {
            atoms: self
                .atoms
                .iter()
                .map(|a| Atom::try_from(a.clone()))
                .collect::<Result<Vec<_>, _>>()?,
            bonds: self
                .bonds
                .iter()
                .map(|b| Bond::try_from(b.clone()))
                .collect::<Result<Vec<_>, _>>()?,
            rings: self.rings.clone(),
            positions: self.positions.clone(),
            multicenter_bonds: self.multicenter_bonds.clone(),
            comments: self.comments.clone(),
            stereo_interpretation: self.stereo_interpretation,
            properties: self.properties.clone(),
            source_format: self.source_format,
        })
    }

    /// View over SGroups from both ctfile_data and cx_data.
    /// Iteration yields ctfile sgroups first, then cx sgroups.
    /// Lookup checks ctfile first, then cx.
    pub fn sgroups(&self) -> SgroupsView<'_> {
        SgroupsView {
            ctfile: self.ctfile_data.as_ref().map(|d| &d.sgroups),
            cx: self.cx_data.as_ref().map(|c| &c.sgroups),
        }
    }

    /// Get mutable reference to SGroups
    /// Uses ctfile_data if present, else cx_data (creating it if needed)
    pub fn sgroups_mut(&mut self) -> &mut BTreeMap<u32, SGroup> {
        if self.ctfile_data.is_some() {
            return &mut self.ctfile_data.as_mut().unwrap().sgroups;
        }
        if self.cx_data.is_none() {
            self.cx_data = Some(CxAnnotationData::default());
        }
        &mut self.cx_data.as_mut().unwrap().sgroups
    }

    /// View over RGroups from both ctfile_data and cx_data.
    /// Iteration yields ctfile rgroups first, then cx rgroups.
    pub fn rgroups(&self) -> RgroupsView<'_> {
        RgroupsView {
            ctfile: self.ctfile_data.as_ref().map(|d| &d.rgroups),
            cx: self.cx_data.as_ref().map(|c| &c.rgroups),
        }
    }

    /// Get mutable reference to RGroups.
    /// Uses ctfile_data if present, else cx_data (creating it if needed)
    pub fn rgroups_mut(&mut self) -> &mut BTreeMap<u32, RGroup> {
        if self.ctfile_data.is_some() {
            return &mut self.ctfile_data.as_mut().unwrap().rgroups;
        }
        if self.cx_data.is_none() {
            self.cx_data = Some(CxAnnotationData::default());
        }
        &mut self.cx_data.as_mut().unwrap().rgroups
    }

    /// Count of SDF/MOL properties
    pub fn property_count(&self) -> usize {
        self.properties.len()
    }

    /// Count of extended atoms (non-Element/NamedIsotope: wildcards, atom lists, rgroups, etc.)
    pub fn extended_atom_count(&self) -> usize {
        self.atoms.iter().filter(|a| a.symbol.is_extended()).count()
    }

    /// Count of extended bonds (query or extended bond orders)
    pub fn extended_bond_count(&self) -> usize {
        self.bonds
            .iter()
            .filter(|b| b.has_extended_features())
            .count()
    }

    /// Count of RGroups defined in CTFile data
    pub fn rgroup_count(&self) -> usize {
        self.rgroups().len()
    }

    /// Count of SGroups defined in CTFile data
    pub fn sgroup_count(&self) -> usize {
        self.sgroups().len()
    }

    fn component_labels(&self) -> Vec<u32> {
        let mut union_find = UnionFind::new(self.atom_count());
        for bond in &self.bonds {
            let (first, second) = bond.atoms.as_tuple();
            union_find.union(first, second);
        }
        for multicenter_bond in &self.multicenter_bonds {
            let indices: Vec<u32> = multicenter_bond
                .contributions()
                .iter()
                .flat_map(|c| c.atoms().iter().copied())
                .collect();
            for (first, second) in indices.into_iter().tuple_combinations() {
                union_find.union(first, second);
            }
        }
        union_find.into_labeling()
    }

    pub fn component_atom_indices(&self) -> Vec<Vec<u32>> {
        let mut components = BTreeMap::new();
        for (atom_index, component_label) in self.component_labels().into_iter().enumerate() {
            components
                .entry(component_label)
                .or_insert_with(Vec::new)
                .push(atom_index as u32);
        }
        components.into_values().collect()
    }

    pub fn split_components(&self) -> Vec<ExtendedMolecule> {
        let mut components = Vec::new();
        for component_atoms in self.component_atom_indices() {
            let atom_index_map: HashMap<u32, u32> = component_atoms
                .iter()
                .copied()
                .enumerate()
                .map(|(new_idx, old_idx)| (old_idx, new_idx as u32))
                .collect();

            let atoms = component_atoms
                .iter()
                .map(|&old_idx| self.atoms[old_idx as usize].clone())
                .collect();

            let mut bonds = Vec::new();
            let mut bond_index_map: HashMap<u32, u32> = HashMap::new();
            for (old_bond_idx, bond) in self.bonds.iter().enumerate() {
                let (a, b) = bond.atoms.as_tuple();
                if let (Some(new_a), Some(new_b)) = (
                    atom_index_map.get(&a).copied(),
                    atom_index_map.get(&b).copied(),
                ) {
                    let new_bond_idx = bonds.len() as u32;
                    bonds.push(bond.update_atoms(new_a, new_b));
                    bond_index_map.insert(old_bond_idx as u32, new_bond_idx);
                }
            }

            let rings =
                self.rings
                    .iter()
                    .filter_map(|ring| {
                        let new_start = match ring.start_atom {
                            Some(old_idx) => Some(*atom_index_map.get(&old_idx)?),
                            None => None,
                        };
                        let new_end = match ring.end_atom {
                            Some(old_idx) => Some(*atom_index_map.get(&old_idx)?),
                            None => None,
                        };
                        if new_start.is_none() && new_end.is_none() {
                            return None;
                        }
                        Some(ring.update_atoms(
                            new_start.unwrap_or_default(),
                            new_end.unwrap_or_default(),
                        ))
                    })
                    .collect();

            let positions = self.positions.as_ref().map(|positions| {
                component_atoms
                    .iter()
                    .filter_map(|&old_idx| positions.get(old_idx as usize).copied())
                    .collect()
            });

            let multicenter_bonds = self
                .multicenter_bonds
                .iter()
                .filter_map(|multicenter| multicenter.update_atoms(&atom_index_map))
                .collect();

            let ctfile_data = self
                .ctfile_data
                .as_ref()
                .and_then(|d| d.update_atoms_bonds(&atom_index_map, &bond_index_map));
            let cx_data = self
                .cx_data
                .as_ref()
                .and_then(|d| d.update_atoms_bonds(&atom_index_map, &bond_index_map));

            components.push(ExtendedMolecule {
                atoms,
                bonds,
                rings,
                positions,
                multicenter_bonds,
                stereo_interpretation: self.stereo_interpretation,
                comments: self.comments.clone(),
                properties: self.properties.clone(),
                ctfile_data,
                cx_data,
                source_format: self.source_format,
            });
        }
        components
    }

    pub fn join_components(components: &[ExtendedMolecule]) -> ExtendedMolecule {
        Self::join_components_inner(components, false).unwrap()
    }

    pub fn try_join_components(
        components: &[ExtendedMolecule],
    ) -> Result<ExtendedMolecule, JoinError> {
        Self::join_components_inner(components, true)
    }

    fn join_components_inner(
        components: &[ExtendedMolecule],
        fail_on_collision: bool,
    ) -> Result<ExtendedMolecule, JoinError> {
        if components.is_empty() {
            return Ok(ExtendedMolecule::empty());
        }

        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let mut rings = Vec::new();
        let mut multicenter_bonds = Vec::new();
        let mut comments = Vec::new();
        let mut properties = IndexMap::new();
        let mut all_positions = true;
        let mut positions_acc = Vec::new();
        let source_format = components[0].source_format;
        let stereo_interpretation = components[0].stereo_interpretation;

        let mut merged_ctfile_data = CtfileData::default();
        let mut has_ctfile_data = false;
        let mut merged_cx_data = CxAnnotationData::default();
        let mut has_cx_data = false;

        for mol in components {
            let atom_offset = atoms.len() as u32;
            atoms.extend(mol.atoms.iter().cloned());

            let mut bond_index_map: HashMap<u32, u32> = HashMap::new();
            for (old_bond_idx, bond) in mol.bonds.iter().enumerate() {
                let (a, b) = bond.atoms.as_tuple();
                bond_index_map.insert(old_bond_idx as u32, bonds.len() as u32);
                bonds.push(bond.update_atoms(a + atom_offset, b + atom_offset));
            }

            for ring in &mol.rings {
                let new_start = ring.start_atom.map(|idx| idx + atom_offset);
                let new_end = ring.end_atom.map(|idx| idx + atom_offset);
                if new_start.is_none() && new_end.is_none() {
                    continue;
                }
                rings.push(
                    ring.update_atoms(new_start.unwrap_or_default(), new_end.unwrap_or_default()),
                );
            }

            let atom_index_map: HashMap<u32, u32> = (0..mol.atom_count() as u32)
                .map(|idx| (idx, idx + atom_offset))
                .collect();
            multicenter_bonds.extend(
                mol.multicenter_bonds
                    .iter()
                    .filter_map(|multicenter| multicenter.update_atoms(&atom_index_map)),
            );

            if let Some(component_positions) = mol.positions.as_ref() {
                positions_acc.extend(component_positions.iter().copied());
            } else {
                all_positions = false;
            }

            comments.extend(mol.comments.iter().cloned());
            for (key, value) in &mol.properties {
                properties.insert(key.clone(), value.clone());
            }

            if let Some(data) = mol
                .ctfile_data
                .as_ref()
                .and_then(|d| d.update_atoms_bonds(&atom_index_map, &bond_index_map))
            {
                has_ctfile_data = true;
                merge_ctfile_data(&mut merged_ctfile_data, data, fail_on_collision)?;
            }

            if let Some(data) = mol
                .cx_data
                .as_ref()
                .and_then(|d| d.update_atoms_bonds(&atom_index_map, &bond_index_map))
            {
                has_cx_data = true;
                merge_cx_data(&mut merged_cx_data, data, fail_on_collision)?;
            }
        }

        Ok(ExtendedMolecule {
            atoms,
            bonds,
            rings,
            positions: if all_positions {
                Some(positions_acc)
            } else {
                None
            },
            multicenter_bonds,
            stereo_interpretation,
            comments,
            properties,
            ctfile_data: if has_ctfile_data {
                Some(merged_ctfile_data)
            } else {
                None
            },
            cx_data: if has_cx_data {
                Some(merged_cx_data)
            } else {
                None
            },
            source_format,
        })
    }
}

fn merge_ctfile_data(
    target: &mut CtfileData,
    source: CtfileData,
    fail_on_collision: bool,
) -> Result<(), JoinError> {
    for (key, value) in source.sgroups {
        if fail_on_collision && target.sgroups.contains_key(&key) {
            return Err(JoinError::CtfileSgroupCollision { label: key });
        }
        target.sgroups.insert(key, value);
    }
    for (key, value) in source.rgroups {
        if fail_on_collision && target.rgroups.contains_key(&key) {
            return Err(JoinError::CtfileRgroupCollision { label: key });
        }
        target.rgroups.insert(key, value);
    }
    target
        .legacy_group_abbreviations
        .extend(source.legacy_group_abbreviations);
    Ok(())
}

fn merge_cx_data(
    target: &mut CxAnnotationData,
    source: CxAnnotationData,
    fail_on_collision: bool,
) -> Result<(), JoinError> {
    for (key, value) in source.stereo_groups {
        target.stereo_groups.insert(key, value);
    }
    if let Some(groups) = source.components {
        target
            .components
            .get_or_insert_with(Vec::new)
            .extend(groups);
    }
    for (key, value) in source.sgroups {
        if fail_on_collision && target.sgroups.contains_key(&key) {
            return Err(JoinError::CxSgroupCollision { label: key });
        }
        target.sgroups.insert(key, value);
    }
    for (key, value) in source.rgroups {
        if fail_on_collision && target.rgroups.contains_key(&key) {
            return Err(JoinError::CxRgroupCollision { label: key });
        }
        target.rgroups.insert(key, value);
    }
    for (key, value) in source.rgroup_members {
        target.rgroup_members.insert(key, value);
    }
    if let Some(entries) = source.local_parity {
        target
            .local_parity
            .get_or_insert_with(Vec::new)
            .extend(entries);
    }
    if let Some(entries) = source.bicyclo_stereo {
        target
            .bicyclo_stereo
            .get_or_insert_with(Vec::new)
            .extend(entries);
    }
    Ok(())
}

impl From<Molecule> for ExtendedMolecule {
    fn from(mol: Molecule) -> Self {
        Self {
            atoms: mol.atoms.into_iter().map(ExtendedAtom::from).collect(),
            bonds: mol.bonds.into_iter().map(ExtendedBond::from).collect(),
            rings: mol.rings,
            positions: mol.positions,
            multicenter_bonds: mol.multicenter_bonds,
            stereo_interpretation: mol.stereo_interpretation,
            comments: mol.comments,
            properties: mol.properties,
            ctfile_data: None,
            cx_data: None,
            source_format: mol.source_format,
        }
    }
}

/// View over SGroups from ctfile_data and cx_data, chained for iteration.
#[derive(Clone, Debug)]
pub struct SgroupsView<'a> {
    ctfile: Option<&'a BTreeMap<u32, SGroup>>,
    cx: Option<&'a BTreeMap<u32, SGroup>>,
}

impl<'a> SgroupsView<'a> {
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn len(&self) -> usize {
        self.ctfile.map(|m| m.len()).unwrap_or(0) + self.cx.map(|m| m.len()).unwrap_or(0)
    }

    pub fn get(&self, key: &u32) -> Option<&'a SGroup> {
        self.ctfile
            .and_then(|m| m.get(key))
            .or_else(|| self.cx.and_then(|m| m.get(key)))
    }

    pub fn contains_key(&self, key: &u32) -> bool {
        self.ctfile.map_or(false, |m| m.contains_key(key))
            || self.cx.map_or(false, |m| m.contains_key(key))
    }

    pub fn iter(&self) -> impl Iterator<Item = (u32, &'a SGroup)> {
        self.ctfile
            .into_iter()
            .flat_map(|m| m.iter())
            .chain(self.cx.into_iter().flat_map(|m| m.iter()))
            .map(|(k, v)| (*k, v))
    }
}

/// View over RGroups from ctfile_data and cx_data, chained for iteration.
#[derive(Clone, Debug)]
pub struct RgroupsView<'a> {
    ctfile: Option<&'a BTreeMap<u32, RGroup>>,
    cx: Option<&'a BTreeMap<u32, RGroup>>,
}

impl<'a> RgroupsView<'a> {
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn len(&self) -> usize {
        self.ctfile.map(|m| m.len()).unwrap_or(0) + self.cx.map(|m| m.len()).unwrap_or(0)
    }

    pub fn get(&self, key: &u32) -> Option<&'a RGroup> {
        self.ctfile
            .and_then(|m| m.get(key))
            .or_else(|| self.cx.and_then(|m| m.get(key)))
    }

    pub fn contains_key(&self, key: &u32) -> bool {
        self.ctfile.map_or(false, |m| m.contains_key(key))
            || self.cx.map_or(false, |m| m.contains_key(key))
    }

    pub fn iter(&self) -> impl Iterator<Item = (u32, &'a RGroup)> {
        self.ctfile
            .into_iter()
            .flat_map(|m| m.iter())
            .chain(self.cx.into_iter().flat_map(|m| m.iter()))
            .map(|(k, v)| (*k, v))
    }
}

#[cfg(test)]
mod tests;
