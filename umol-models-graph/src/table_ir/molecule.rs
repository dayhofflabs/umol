//! Molecule IR for TableIR.
//!
//! `Molecule` is TableIR molecule representation for molecules with fixed composition
//! `ExtendedMolecule` is temporary container for generalized molecules (query molecules, molecule
//!   libraries, substitutions, polymers, etc.). Currently used to hold extensions from CTFile formats
//!   (Query features, SGroups, RGroups, etc.). It is supposed to be split into multiple semantically
//!   defined structures.

use std::collections::BTreeMap;
use std::sync::LazyLock;

use indexmap::IndexMap;
use umol_data::Element;

use super::atom::{Atom, AtomSymbol, ExtendedAtom};
use super::bond::{Bond, ExtendedBond};
use super::ctfile_data::CtfileData;
use super::rgroup::RGroup;
use super::sgroup::SGroup;
use super::source::SourceFormat;
use super::topology::{Fragment, Link, Ring};
use super::utils::{element_symbol_key, format_sum_formula};
use crate::position::Point3D;

/// Basic molecule IR
#[derive(Clone, Debug, PartialEq)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub rings: Vec<Ring>,
    pub positions: Option<Vec<Point3D>>,
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

    /// Count of SDF/MOL properties
    pub fn property_count(&self) -> usize {
        self.properties.len()
    }
}

/// Extended molecule IR - includes MDL extensions (SGroups, RGroups, etc.)
/// This is a flat structure using ExtendedAtom and ExtendedBond.
#[derive(Clone, Debug, PartialEq)]
pub struct ExtendedMolecule {
    // Core structure with extended atom/bond types
    pub atoms: Vec<ExtendedAtom>,
    pub bonds: Vec<ExtendedBond>,
    pub rings: Vec<Ring>,
    pub positions: Option<Vec<Point3D>>,

    // Additional structure
    pub fragments: Vec<Fragment>,
    pub links: Vec<Link>,
    pub electrons: Option<u32>,

    // Properties and metadata
    pub comments: Vec<String>,
    pub properties: IndexMap<String, String>,

    // Format-specific data for roundtripping (CTFile formats)
    pub ctfile_data: Option<CtfileData>,

    pub source_format: SourceFormat,
}

impl ExtendedMolecule {
    pub fn empty() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            rings: Vec::new(),
            positions: None,
            fragments: Vec::new(),
            links: Vec::new(),
            electrons: None,
            comments: Vec::new(),
            properties: IndexMap::new(),
            ctfile_data: None,
            source_format: SourceFormat::UNKNOWN,
        }
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Get sum formula in Hill notation (C first, H second, then alphabetically)
    pub fn sum_formula(&self) -> String {
        let mut atom_counts: BTreeMap<[u8; 2], (Element, usize)> = BTreeMap::new();
        let mut c_count = 0usize;
        let mut h_count = 0usize;
        let mut charge = 0i32;

        for atom in &self.atoms {
            let element = match &atom.symbol {
                AtomSymbol::Element(e) => Some(*e),
                AtomSymbol::NamedIsotope(i) => Some(i.element()),
                _ => None,
            };
            if let Some(element) = element {
                match element {
                    Element::C => c_count += 1,
                    Element::H => h_count += 1,
                    e => {
                        let key = element_symbol_key(e);
                        atom_counts.entry(key).or_insert((e, 0)).1 += 1;
                    }
                }
            }
            if let Some(ch) = atom.charge {
                charge += ch as i32;
            }
        }

        format_sum_formula(c_count, h_count, atom_counts, charge)
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
            comments: self.comments.clone(),
            properties: self.properties.clone(),
            source_format: self.source_format,
        })
    }

    /// Get reference to SGroups (CTFile-specific)
    /// Returns empty map if ctfile_data is not present
    pub fn sgroups(&self) -> &std::collections::BTreeMap<usize, SGroup> {
        static EMPTY: LazyLock<BTreeMap<usize, SGroup>> = LazyLock::new(BTreeMap::new);

        self.ctfile_data
            .as_ref()
            .map(|d| &d.sgroups)
            .unwrap_or(&*EMPTY)
    }

    /// Get mutable reference to SGroups (CTFile-specific)
    /// Initializes ctfile_data if not present
    pub fn sgroups_mut(&mut self) -> &mut std::collections::BTreeMap<usize, SGroup> {
        if self.ctfile_data.is_none() {
            self.ctfile_data = Some(CtfileData::default());
        }
        &mut self.ctfile_data.as_mut().unwrap().sgroups
    }

    /// Get reference to RGroups (CTFile-specific)
    /// Returns empty map if ctfile_data is not present
    pub fn rgroups(&self) -> &std::collections::BTreeMap<usize, RGroup> {
        static EMPTY: LazyLock<BTreeMap<usize, RGroup>> = LazyLock::new(BTreeMap::new);

        self.ctfile_data
            .as_ref()
            .map(|d| &d.rgroups)
            .unwrap_or(&*EMPTY)
    }

    /// Get mutable reference to RGroups (CTFile-specific)
    /// Initializes ctfile_data if not present
    pub fn rgroups_mut(&mut self) -> &mut std::collections::BTreeMap<usize, RGroup> {
        if self.ctfile_data.is_none() {
            self.ctfile_data = Some(CtfileData::default());
        }
        &mut self.ctfile_data.as_mut().unwrap().rgroups
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
        self.bonds.iter().filter(|b| b.has_extended_features()).count()
    }

    /// Count of RGroups defined in CTFile data
    pub fn rgroup_count(&self) -> usize {
        self.rgroups().len()
    }

    /// Count of SGroups defined in CTFile data
    pub fn sgroup_count(&self) -> usize {
        self.sgroups().len()
    }
}

impl From<Molecule> for ExtendedMolecule {
    fn from(mol: Molecule) -> Self {
        Self {
            atoms: mol.atoms.into_iter().map(ExtendedAtom::from).collect(),
            bonds: mol.bonds.into_iter().map(ExtendedBond::from).collect(),
            rings: mol.rings,
            positions: mol.positions,
            fragments: Vec::new(),
            links: Vec::new(),
            electrons: None,
            comments: mol.comments,
            properties: mol.properties,
            ctfile_data: None,
            source_format: mol.source_format,
        }
    }
}

#[cfg(test)]
mod tests;
