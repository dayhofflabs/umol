//! Molecule IR for TableIR.
//!
//! `Molecule` is TableIR molecule representation for molecules with fixed composition
//! `ExtendedMolecule` is temporary container for generalized molecules (query molecules, molecule
//!   libraries, substitutions, polymers, etc.). Currently used to hold extensions from CTFile formats
//!   (Query features, SGroups, RGroups, etc.). It is supposed to be split into multiple semantically
//!   defined structures.

use std::collections::BTreeMap;
use std::sync::LazyLock;

use umol_data::Element;

use super::atom::{Atom, AtomSymbol, ExtendedAtom};
use super::bond::{Bond, ExtendedBond};
use super::ctfile_data::CtfileData;
use super::property::Property;
use super::rgroup::RGroup;
use super::sgroup::SGroup;
use super::source::SourceFormat;
use super::topology::{Fragment, Link, Ring};
use super::utils::{element_symbol_key, format_sum_formula};

/// Basic molecule IR
#[derive(Clone, Debug, PartialEq)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub rings: Vec<Ring>,

    pub source_format: SourceFormat,
}

impl Molecule {
    pub fn empty() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            rings: Vec::new(),
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
}

/// Extended molecule IR - includes MDL extensions (SGroups, RGroups, etc.)
/// This is a flat structure using ExtendedAtom and ExtendedBond.
#[derive(Clone, Debug, PartialEq)]
pub struct ExtendedMolecule {
    // Core structure with extended atom/bond types
    pub atoms: Vec<ExtendedAtom>,
    pub bonds: Vec<ExtendedBond>,
    pub rings: Vec<Ring>,
    pub source_format: SourceFormat,

    // Additional structure
    pub fragments: Vec<Fragment>,
    pub links: Vec<Link>,
    pub electrons: Option<u32>,

    // Properties and metadata
    pub properties: Vec<Property>,
    pub comments: Vec<String>,

    // Format-specific data for roundtripping (CTFile formats)
    pub ctfile_data: Option<CtfileData>,
}

impl ExtendedMolecule {
    pub fn empty() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            rings: Vec::new(),
            source_format: SourceFormat::UNKNOWN,
            fragments: Vec::new(),
            links: Vec::new(),
            electrons: None,
            properties: Vec::new(),
            comments: Vec::new(),
            ctfile_data: None,
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
}

impl From<Molecule> for ExtendedMolecule {
    fn from(mol: Molecule) -> Self {
        Self {
            atoms: mol.atoms.into_iter().map(ExtendedAtom::from).collect(),
            bonds: mol.bonds.into_iter().map(ExtendedBond::from).collect(),
            rings: mol.rings,
            source_format: mol.source_format,
            fragments: Vec::new(),
            links: Vec::new(),
            electrons: None,
            properties: Vec::new(),
            comments: Vec::new(),
            ctfile_data: None,
        }
    }
}

#[cfg(test)]
mod tests;
