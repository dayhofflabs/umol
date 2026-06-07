//! Molecule types for TableIR.
//!
//! `Molecule` is TableIR molecule representation for molecules with fixed composition
//! `ExtendedMolecule` is temporary container for generalized molecules (query molecules, molecule
//!   libraries, substitutions, polymers, etc.). Currently used to hold extensions from CTFile formats
//!   (Query features, SGroups, RGroups, etc.). It is supposed to be split into multiple semantically
//!   defined structures.

use std::collections::BTreeMap;

use indexmap::IndexMap;
use umol_shared::element::Element;

use super::atom::{Atom, AtomSymbol, ExtendedAtom};
use super::bond::{Bond, ExtendedBond};
use super::ctfile_data::CtfileData;
use super::cx_data::CxAnnotationData;
use super::error::ConversionError;
use super::multicenter::MulticenterBond;
use super::rgroup::RGroup;
use super::sgroup::SGroup;
use super::source::SourceFormat;
use super::stereo::{ChiralityFrame, ConfigurationScope};
use super::utils::{element_symbol_key, format_sum_formula};
use crate::position::Point3D;

/// Basic molecule IR
#[derive(Clone, Debug, PartialEq)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub positions: Option<Vec<Point3D>>,
    pub multicenter_bonds: Vec<MulticenterBond>,
    pub configuration_scope: Option<ConfigurationScope>,
    pub chirality_frame: Option<ChiralityFrame>,
    pub comments: Vec<String>,
    pub properties: IndexMap<String, String>,
    pub source_format: SourceFormat,
}

impl Molecule {
    pub fn empty() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            positions: None,
            multicenter_bonds: Vec::new(),
            configuration_scope: None,
            chirality_frame: None,
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
    pub positions: Option<Vec<Point3D>>,
    pub multicenter_bonds: Vec<MulticenterBond>,
    pub configuration_scope: Option<ConfigurationScope>,
    pub chirality_frame: Option<ChiralityFrame>,
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
            positions: None,
            multicenter_bonds: Vec::new(),
            configuration_scope: None,
            chirality_frame: None,
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
        if let Some(data) = self.ctfile_data.as_mut() {
            &mut data.sgroups
        } else {
            self.cx_data = Some(CxAnnotationData::default());
            &mut self.cx_data.as_mut().unwrap().sgroups
        }
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
        if let Some(data) = self.ctfile_data.as_mut() {
            &mut data.rgroups
        } else {
            self.cx_data = Some(CxAnnotationData::default());
            &mut self.cx_data.as_mut().unwrap().rgroups
        }
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
}

impl From<Molecule> for ExtendedMolecule {
    fn from(mol: Molecule) -> Self {
        Self {
            atoms: mol.atoms.into_iter().map(ExtendedAtom::from).collect(),
            bonds: mol.bonds.into_iter().map(ExtendedBond::from).collect(),
            positions: mol.positions,
            multicenter_bonds: mol.multicenter_bonds,
            configuration_scope: mol.configuration_scope,
            chirality_frame: mol.chirality_frame,
            comments: mol.comments,
            properties: mol.properties,
            ctfile_data: None,
            cx_data: None,
            source_format: mol.source_format,
        }
    }
}

impl TryFrom<ExtendedMolecule> for Molecule {
    type Error = ConversionError;

    fn try_from(extended: ExtendedMolecule) -> Result<Self, Self::Error> {
        Ok(Molecule {
            atoms: extended
                .atoms
                .iter()
                .map(|a| Atom::try_from(a.clone()))
                .collect::<Result<Vec<_>, _>>()?,
            bonds: extended
                .bonds
                .iter()
                .map(|b| Bond::try_from(b.clone()))
                .collect::<Result<Vec<_>, _>>()?,
            positions: extended.positions.clone(),
            multicenter_bonds: extended.multicenter_bonds.clone(),
            comments: extended.comments.clone(),
            configuration_scope: extended.configuration_scope,
            chirality_frame: extended.chirality_frame,
            properties: extended.properties.clone(),
            source_format: extended.source_format,
        })
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
        self.ctfile.is_some_and(|m| m.contains_key(key))
            || self.cx.is_some_and(|m| m.contains_key(key))
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
        self.ctfile.is_some_and(|m| m.contains_key(key))
            || self.cx.is_some_and(|m| m.contains_key(key))
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
