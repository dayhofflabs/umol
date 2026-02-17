//! Reaction IR for TableIR.
//!
//! `Reaction` is TableIR reaction representation for reactions with fixed composition
//! Delaying the implementation for `ExtendedReaction` as long its semantics is unclear.
//! (ordinary reactions between extended molecules (variable composition, polymers, etc.)
//! and reaction templates are distinct use cases).

use std::collections::BTreeMap;

use indexmap::IndexMap;

use super::molecule::{ExtendedMolecule, Molecule};
use super::source::SourceFormat;

/// Basic reaction IR
#[derive(Clone, Debug, PartialEq)]
pub struct Reaction {
    pub reactants: Molecule,
    pub products: Molecule,
    pub agents: Molecule,
    pub atom_mapping: BTreeMap<u32, (Vec<u32>, Vec<u32>)>,

    pub comments: Vec<String>,
    pub properties: IndexMap<String, String>,
    pub source_format: SourceFormat,
}

impl Reaction {
    pub fn empty() -> Self {
        Self {
            reactants: Molecule::empty(),
            products: Molecule::empty(),
            agents: Molecule::empty(),
            atom_mapping: BTreeMap::new(),
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::UNKNOWN,
        }
    }

    pub fn from_molecules(reactants: Molecule, products: Molecule, agents: Molecule) -> Self {
        Self {
            reactants,
            products,
            agents,
            atom_mapping: BTreeMap::new(),
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::UNKNOWN,
        }
    }
}

/// Extended reaction IR
#[derive(Clone, Debug, PartialEq)]
pub struct ExtendedReaction {
    pub reactants: ExtendedMolecule,
    pub products: ExtendedMolecule,
    pub agents: ExtendedMolecule,
    pub atom_mapping: BTreeMap<u32, (Vec<u32>, Vec<u32>)>,

    pub comments: Vec<String>,
    pub properties: IndexMap<String, String>,
    pub source_format: SourceFormat,
}

impl ExtendedReaction {
    pub fn empty() -> Self {
        Self {
            reactants: ExtendedMolecule::empty(),
            products: ExtendedMolecule::empty(),
            agents: ExtendedMolecule::empty(),
            atom_mapping: BTreeMap::new(),
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::UNKNOWN,
        }
    }

    pub fn from_extended_molecules(
        reactants: ExtendedMolecule,
        products: ExtendedMolecule,
        agents: ExtendedMolecule,
    ) -> Self {
        Self {
            reactants,
            products,
            agents,
            atom_mapping: BTreeMap::new(),
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::UNKNOWN,
        }
    }
}
